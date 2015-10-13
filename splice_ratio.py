import parsers as ps


class SpliceRatioCounter:
    def __init__(self, bam_paths='', reads_orientation='forward', min_qual=10, max_mm_indels=2, min_overlap=10):
        """
        responsible of computing the count table of every junction (junction, intron, invalid counts)
        :param bam_path: must be indexed, otherwise an index will be built, provided samtools is available
        :param reads_orientation: must be 'reverse' if reads (or first reads if PE) are on the opposite strand of gene
        :param min_qual: minimal mapping quality of a read to be considered
        :param max_mm_indels: maximal number of mismatches (non-matches) and indels of a read to be considered in an intron
        :param min_overlap: minimal nucleotide overlap of a read in an intron to be considered in it
        """
        self.min_overlap, self.max_mm_indels = min_overlap, max_mm_indels
        self.bams, self.splices_dic, self.start_dic, self.last_read_len = [], {}, {}, {}
        for bam_path in bam_paths:
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
        for site, junction_count in splice_sites:
            if site not in self.splices_dic:
                self.splices_dic[site] = []
                empty_dic = dict(zip(['junction', 'intron', 'invalid', 'surrounding'], [0] * 4))
                for bam in range(self.n_bams):
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
    def __init__(self, sample_conditions, splice_ratio_counts, min_counts=20, max_invalid=0.05, max_unequal=10, n_bins=10, max_uncovered=0.1):
        """
        responsible of filtering junction counts to keep only high confidence intron (or pre-mRNA) reads
        :param sample_conditions:
        :param splice_ratio_counts:
        :param min_counts:
        :param max_invalid:
        :param max_unequal:
        :param n_bins:
        :param max_uncovered:
        :return:
        """
        self.condition = _parse_factor(sample_conditions)
        self.splice_ratio_counts = splice_ratio_counts
        self.min_counts = min_counts
        self.max_invalid = max_invalid
        self.max_unequal = max_unequal
        self.max_uncovered = max_uncovered
        self.n_bins = n_bins

    def filter_sites(self):
        """
        removes sites from self.splices_dic if reads are unlikely to be
        coming from pre-mRNA or intron retention (...)
        """
        for site in self.splice_ratio_counts.splices_dic.keys():
            counts = self.splice_ratio_counts.splices_dic[site]
            starts = self.splice_ratio_counts.start_dic[site]
            last_read_len = self.splice_ratio_counts.last_read_len[site]
            if self._few_valid_or_many_invalid(counts) or \
                    self._not_all_covered(site, starts, last_read_len) or \
                    self._unequal_cov(starts):
                del self.splice_ratio_counts.splices_dic[site]

    def _few_valid_or_many_invalid(self, counts):
        """
        returns True if there are too few valid or too many invalid reads in site
        """
        for condition_i in self.condition.keys():
            cond_sum_i, cond_sum_j, cond_sum_n = 0, 0, 0
            for sample in self.condition[condition_i]:
                cond_sum_j += counts[sample]['junction']
                cond_sum_i += counts[sample]['intron']
                cond_sum_n += counts[sample]['invalid']
            if min(cond_sum_j, cond_sum_i) < self.min_counts or cond_sum_n > self.max_invalid:
                return True
        return False

    def _unequal_cov(self, starts):
        """
        returns True if coverage of the site is too unequal (defined by max_unequal)
        """
        bin_counts = self._fill_bins(starts)
        mean_bin_count = sum(bin_counts)/len(bin_counts)
        chi_test = self._chi_test(bin_counts, mean_bin_count)
        if chi_test / mean_bin_count > self.max_unequal:
            return True
        return False

    def _chi_test(self, obs, exp):
        return sum([((x - exp)**2) / exp for x in obs])

    def _fill_bins(self, starts):
        first_start = starts[0]
        len_starts = starts[len(starts) - 1] - first_start
        bin_width = len_starts / self.n_bins
        bins, bin_count = [], 0
        for start in starts:
            dist_from_first = start - first_start
            if dist_from_first < bin_width:
                bin_count += 1
            else:
                bins.append(bin_count)
                first_start += bin_width
                bin_count = 1
        return bins

    def _not_all_covered(self, site, starts, last_read_len):
        """
        returns True if coverage of the site is not completely covered (up to max_uncovered%)
        """
        site_dic = _parse_site(site)
        end_last = starts[len(starts) - 1] + last_read_len
        start_first = starts[0]
        site_len = site_dic['stop'] - site_dic['start']
        max_uncovered_bases = site_len * self.max_uncovered
        if  start_first > site_dic['start'] + max_uncovered_bases or \
            end_last < site_dic['stop'] - max_uncovered_bases:
            return True
        return False

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
    samples_dic = {}
    for i,k in enumerate(factor):
        if k not in samples_dic:
            samples_dic[k] = []
        samples_dic[k].append(i)
    assert len(samples_dic) > 1
    # add exception!
    return samples_dic
