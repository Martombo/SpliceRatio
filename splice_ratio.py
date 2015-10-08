import parsers as ps


class SpliceRatio:
    def __init__(self, bam_paths='', reads_orientation='forward', min_qual=10):
        """
        :param bam_path: must be indexed, otherwise an index will be built, provided samtools is available
        :param reads_orientation: must be 'reverse' if reads (or first reads if PE) are on the opposite strand of gene
        """
        self.bams = []
        self.n_bams = len(self.bams)
        self.splices_dic = {}
        for bam_path in bam_paths:
            self.bams.append(ps.Bam(bam_path, reads_orientation))
        for n_bam in range(self.n_bams):
            splice_sites = self.bams[n_bam].get_splice_sites(min_qual)
            self._add_splice_sites(splice_sites, n_bam)
        self.min_qual = min_qual

    def get_splices_coverage(self, max_mm=0, min_overlap=10):
        """
        get the number of reads from within all the splice sites
        likely to be originating from the retained intron
        :param min_qual: default for TopHat: only uniquely mapped reads
        :param max_mm: max tolerated mismatches, indels on read
        :param min_overlap: minimal bases overlap between read and site
        """
        self.min_overlap, self.max_mm = min_overlap, max_mm
        for site_str in self.splices_dic.keys():
            site = self._parse_site(site_str)
            for n_bam in range(self.n_bams):
                (count, na) = self._splice_counts(site, self.bams[n_bam])
                self.splices_dic[site_str][n_bam]['intron'] = count
                self.splices_dic[site_str][n_bam]['invalid'] = na

    def filter_sites(self, max_sites_within=0):
        for site in self.splices_dic.keys():
            if self._too_many_invalid(site):
                del self.splices_dic[site]
            elif self._unequal_cov(site):
                del self.splices_dic[site]

    def _add_splice_sites(self, splice_sites, n_bam):
        for k,v in splice_sites:
            if k not in self.splices_dic:
                self.splices_dic[k] = []
                for bam in range(self.n_bams):
                    self.splices_dic[k].append({'junction':0})
            self.splices_dic[k][n_bam]['junction'] = v

    def _parse_site(self, site_str):
        site_dic = dict(zip(('chrom', 'start', 'stop', 'strand'), site_str.split('_')))
        site_dic['start'] = int(site_dic['start'])
        site_dic['stop'] = int(site_dic['stop'])
        return site_dic

    def _splice_counts(self, site, bam):
        """
        returns the number of reads from within one given site
        likely to be originating from the retained intron
        and the number of invalid reads
        :param site: dict with chrom, start, stop, strand
        """
        n_reads, n_na = 0, 0
        for read in bam.fetch(site['chrom'], site['start'], site['stop']):
            flag = self._inspect_read(read, site, bam)
            if flag <= 1:
                n_reads += flag
            else:
                n_na += 1
        return (n_reads, n_na)

    def _inspect_read(self, read, site, bam):
        read_splicer = ps.Read(read.cigar, read.reference_start)
        if read.mapq >= self.min_qual and \
                        site['strand'] == bam.determine_strand(read) and \
                        read.get_mismatches() <= self.max_mm and \
                        self._overlap(read, site) >= self.min_overlap:
            return 1
        else:
            return 0

    def _read_in_site(self, read, site):
        """
        checks if a read is fully inside region and a good enough match
        :param read: read to be checked (pysam.AlignedRead)
        :param start: start of the splice site
        :param stop: end of splice site
        :return: flag: 0 not ,1 is in site 2: Na, other splice site within site thus biased coverage
        """
        mm = 0
        for c_part in read.cigar:
            if c_part[0] != 0:
                mm += c_part[1]
            if c_part[0] == 3 and c_part[1] > 2:
                if self._site_in_read_splice(read, site):
                    return 0
                return 2
        overlap = self._overlap(read, site)
        if overlap < self.min_overlap or mm > self.max_mm:
            return 0
        return 1

    def _overlap(self, read, site):
        """
        computes the number of bases of overlap between a read and a site
        """
        read_len = read.query_length
        overlap = read_len
        if site['start'] > read.reference_start:
            overlap = read_len - (site['start'] - read.reference_start)
        if site['stop'] < read.reference_start + read_len:
            overlap = min(site['stop'] - read.reference_start, overlap)
        return overlap

    def _site_in_read_splice(self, read, site):
        """
        checks if the splice site is completely inside a read non-match (another splice site)
        """
        read_splicer = ps.Read(read.cigar, read.reference_start)
        read_splice_sites = read_splicer.get_splice_sites()
        for read_site in read_splice_sites:
            if read_site[0] <= site['start'] and read_site[1] >= site['stop']:
                return True
        return False
