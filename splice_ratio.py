import parsers as ps


class SpliceRatioCounter:
    def __init__(self, bam_paths='', reads_orientation='forward', min_qual=10, max_mm_indels=2, min_overlap=10,
                 test=False):
        """
        responsible of computing the count table of every junction (junction, intron, invalid counts)
        :param bam_paths: blank space separated list of bam files paths (if no index, one is built with samtools)
        :param reads_orientation: must be 'reverse' if reads (or first reads if PE) are on the opposite strand of gene
        :param min_qual: minimal mapping quality of a read to be considered
        :param max_mm_indels: maximal number of non-matches and indels of a read to be considered in an intron
        :param min_overlap: minimal nucleotide overlap of a read in an intron to be considered in it
        :param test: for testing
        """
        if not test:
            assert bam_paths
        self.min_overlap, self.max_mm_indels = min_overlap, max_mm_indels
        self.bams, self.splices_dic, self.start_dic, self.last_read_len = [], {}, {}, {}
        self.arrangements = ['junction', 'intron', 'invalid', 'surrounding']
        for bam_path in bam_paths.split(' '):
            if bam_path:
                self.bams.append(ps.Bam(bam_path, reads_orientation))
        self.n_bams = len(self.bams)
        for n_bam in range(self.n_bams):
            splice_sites = self.bams[n_bam].get_splice_sites(min_qual)
            self._add_splice_sites(splice_sites, n_bam)
        self.min_qual = min_qual

    def get_splices_coverage(self):
        """
        get the number of reads from within all the splice sites
        likely to be originating from the retained intron
        """
        for site_str in self.splices_dic.keys():
            site = _parse_site(site_str)
            for n_bam in range(self.n_bams):
                for read in self.bams[n_bam].fetch(site['chrom'], site['start'], site['stop']):
                    read_arrangement = self._read_arrangement(site_str, read)
                    self.splices_dic[site_str][n_bam][read_arrangement] += 1

    def _add_splice_sites(self, splice_sites, n_bam):
        """
        populates self.splices_dic and self.start_dic
        :param splice_sites: bam splice sites dict, from Bam().get_splice_sites()
        :param n_bam: bam index
        """
        for site, junction_count in splice_sites.items():
            if site not in self.splices_dic:
                self.splices_dic[site] = []
                for bam in range(self.n_bams):
                    empty_dic = dict(zip(self.arrangements, [0] * 4))
                    self.splices_dic[site].append(empty_dic)
            self.splices_dic[site][n_bam]['junction'] = junction_count

    def _read_arrangement(self, site_str, read):
        """
        determines the position of the read relative to the site:
        'surrounding': the splice site is completely inside a read non-match (another splice site)
        'intron': the read is mostly inside the splice site (min_overlap)
        'invalid': the read overlaps an exon or splices outside the site
        :param site_str: 'chr1_18315_29387_+'
        :param read: pysam.AlignedRead
        """
        site = _parse_site(site_str)
        reader = ps.Read(read.cigar, read.reference_start)
        if self._read_surrounds_site(reader, site):
            return 'surrounding'
        if read.mapq >= self.min_qual and \
                        site['strand'] == self.bams[0].determine_strand(read) and \
                        len(reader.mm_indels) <= self.max_mm_indels and \
                        self._overlap(read, site) >= self.min_overlap:
            self.start_dic[site_str].append(read.reference_start)
            self.last_read_len[site_str] = read.query_length
            return 'intron'
        return 'invalid'

    @staticmethod
    def _read_surrounds_site(reader, site):
        """
        checks if the splice site is completely inside a read non-match (another splice site)
        """
        for read_site in reader.splice_sites:
            if read_site[0] <= site['start'] and read_site[1] >= site['stop']:
                return True
        return False

    @staticmethod
    def _overlap(read, site):
        """
        computes the number of bases of overlap between a read and a site
        """
        read_len = read.query_length
        overlap = read_len
        if site['start'] > read.reference_start:
            overlap = read_len - (site['start'] - read.reference_start)
        if site['stop'] < read.reference_start + read_len:
            overlap = min(site['stop'] - read.reference_start, overlap)
        return max(overlap, 0)


class SpliceRatioFilter:
    def __init__(self, splice_ratio_counts=None, sample_conditions=None, min_valid=20, min_covered=0.8,
                 max_invalid=0.05, n_bins=5, max_unequal=5, test=False):
        """
        responsible of filtering junction counts to keep only high confidence intron (or pre-mRNA) reads
        :param splice_ratio_counts: SpliceRatioCounter object, on which get_splices_coverage was run
        :param sample_conditions: sample condition factor eg: ['A', 'A', 'B', 'B']
        :param min_valid: minimal number of valid reads (both junction, intron) present in the totals of each condition
        :param min_covered: minimal % of the intron which is covered with reads, samples of same condition aggregated
        :param max_invalid: maximal number of invalid reads present in the totals of each condition
        :param n_bins: number of bins in which the valid reads covering the intron are divided
        :param max_unequal: maximal value of chi_square_test / mean(bins)
        :param test: for testing
        """
        if not test:
            assert splice_ratio_counts and sample_conditions
        self.condition = _parse_factor(sample_conditions)
        self.splice_ratio_counts = splice_ratio_counts
        self.min_int_jun = min_valid
        self.max_invalid = max_invalid
        self.max_unequal = max_unequal
        self.max_uncovered = (1 - min_covered) / 2.0
        self.n_bins = n_bins
        self.masked_sites, self.licit_sites = {}, {}

    def filter_sites(self):
        """
        sites from splices_dic are added to licit_sites if reads are likely to be
        coming from pre-mRNA or intron retention (...)
        otherwise they are added to masked_sites
        """
        for site in self.splice_ratio_counts.splices_dic.keys():
            counts = self.splice_ratio_counts.splices_dic[site]
            starts = self.splice_ratio_counts.start_dic[site]
            last_read_len = self.splice_ratio_counts.last_read_len[site]
            invalids = self._few_valid_or_many_invalid(counts)
            covered = self._all_covered(site, starts, last_read_len)
            unequal = self._unequal_cov(starts)
            if invalids or not covered or unequal:
                self.masked_sites[site] = [invalids, not covered, unequal]
            else:
                self.licit_sites[site] = counts

    def table(self, case):
        """
        returns a table (string) with the information about sites specified in case:
        'chr1_34098_44019_+ 10 14 23 9
        chr12_21385_89144_- 21 49 20 0'
        :param case: must be one of: 'junction', 'intron', 'invalid', 'surrounding', 'masked'
        :rtype str
        """
        string = ''
        if case == 'masked':
            for site, checks in self.masked_sites.items():
                string += site + ' ' + ' '.join([str(x) for x in checks]) + '\n'
        elif case in self.splice_ratio_counts.arrangements:
            for site, site_dic in self.licit_sites.items():
                string += site + ' ' + ' '.join(str(x[case]) for x in site_dic) + '\n'
        return string

    def _few_valid_or_many_invalid(self, counts):
        """
        returns True if there are too few valid or too many invalid reads in site
        """
        for condition_i in self.condition.keys():
            intron_count, junction_count, invalid_count = 0, 0, 0
            for sample in self.condition[condition_i]:
                junction_count += counts[sample]['junction']
                intron_count += counts[sample]['intron']
                invalid_count += counts[sample]['invalid']
            if min(junction_count, intron_count) < self.min_int_jun or invalid_count > self.max_invalid:
                return True
        return False

    def _unequal_cov(self, starts):
        """
        returns True if coverage of the site is too unequal (defined by max_unequal)
        """
        assert len(starts) >= self.min_int_jun
        bin_counts = self._fill_bins(starts)
        mean_bin_count = sum(bin_counts) / len(bin_counts)
        chi_test = _chi_test(bin_counts, mean_bin_count)
        if chi_test / mean_bin_count > self.max_unequal:
            return True
        return False

    def _all_covered(self, site_str, starts, last_read_len):
        """
        returns True if the coverage of the site is completely covered (up to max_uncovered%)
        """
        site_dic = _parse_site(site_str)
        end_last = starts[len(starts) - 1] + last_read_len
        start_first = starts[0]
        site_len = site_dic['stop'] - site_dic['start']
        max_uncovered_bases = site_len * self.max_uncovered
        if start_first > site_dic['start'] + max_uncovered_bases or \
                end_last < site_dic['stop'] - max_uncovered_bases:
            return False
        return True

    def _fill_bins(self, starts):
        """
        divides a list of numbers into deciles
        """
        starts = sorted(starts)
        bin_width = (starts[len(starts) - 1] - starts[0]) / float(self.n_bins)
        if bin_width == 0:
            return [len(starts)] + [0] * (self.n_bins - 1)
        bins, bin_count = [], 0
        bin_pos = starts[0] + bin_width
        for start in starts:
            while start > bin_pos:
                bins.append(bin_count)
                bin_pos += bin_width
                bin_count = 0
            bin_count += 1
        bins.append(bin_count)
        return bins


def _parse_site(site_str):
    """
    :param site_str: splice site string 'chr1_1024_1560_+'
    :return: dict: chrom, start, stop, strand
    """
    site_dic = dict(zip(['chrom', 'start', 'stop', 'strand'], site_str.split('_')))
    site_dic['start'] = int(site_dic['start'])
    site_dic['stop'] = int(site_dic['stop'])
    assert site_dic['stop'] > site_dic['start']
    return site_dic


def _parse_factor(factor):
    if not factor:
        return None
    samples_dic = {}
    for i, k in enumerate(factor):
        if k not in samples_dic:
            samples_dic[k] = []
        samples_dic[k].append(i)
    assert len(samples_dic) > 1
    return samples_dic


def _chi_test(obs, exp):
    return sum([((x - exp) ** 2) / exp for x in obs])
