import parsers as prs
import pysam as ps


class SpliceRatio:

    def __init__(self, bam_path='', reads_orientation='forward'):
        self.bam = prs.Bam(bam_path, reads_orientation=reads_orientation)
        self.bam_file = ps.AlignmentFile(bam_path,'rb')
        self.splices_dic = {}

    def get_splicing_coverage(self, min_qual = 40, max_mm=0, max_sites_within = 0, min_overlap = 10):
        """
        get the number of reads from within splice site
        :param chrom: str chromosome name
        :param start: int start
        :param stop: int stop
        :param min_qual: default TopHat: only uniquely mapped reads
        :param max_mm: max tolerated mismatches on read ( including indel
        :param min_overlap number of nucleotides inside splice junctions
        """
        for splice_site in self.splices_dic.keys():
            (chrom, start, stop, strand) = splice_site.split('_')
            start = int(start)
            stop = int(stop)
            (count, na) = self._count_coverage(chrom, start, stop, strand, min_qual, min_overlap)
            if na > max_sites_within:
                self.splices_dic[splice_site].append('Na')
            else:
                self.splices_dic[splice_site].append(count)

    def _count_coverage(self, chrom, start, stop, strand,  min_qual, min_overlap):
        """
        returns the number of reads from within given region
        reads is only counted if fully inside and with less than max_mm mismatches
        """
        fetch = self.bam_file.fetch(chrom, start, stop)
        n_reads = 0
        n_na = 0
        for read in fetch:
            if read.mapq  >= min_qual and strand == self.bam._determine_strand(read) and not self._is_around_site(read, start, stop):
                count = self._read_in_intron(read, start, stop, min_overlap)
                if  count <= 1:
                    n_reads += count
                else:
                    n_na += 1
        return (n_reads, n_na)

    def _read_in_intron(self, read, start, stop, min_overlap):
        """
        checks if a read is fully inside region and a good enough match
        :param read: read of interest, pysam object
        :param start: satrt of region
        :param stop: end of region
        :param max_mm: tolerate number of mismatches
        :return: is_in_intron 0 not ,1 is in site 3: Na, other splice site within site thus biased coverage
        """
        is_in_intron = 1
        leng = 0
        mm = 0
        has_intron = False

        for c_part in read.cigar:
            if c_part[0] != 1: # insertion = 1 no extension on reference
                leng += c_part[1]
            if c_part[0] != 0:
                mm += c_part[1]
            if c_part[0] == 3:
                has_intron = True
        if start > read.reference_start: # read befor start
            overlap = leng - (start - read.reference_start)
        elif stop < read.reference_start + leng : # read longer than junction
            overlap = stop - read.reference_start
        else:
            overlap = leng
        if overlap < min_overlap:
            is_in_intron = 0
        if has_intron:
            is_in_intron = 2
        return is_in_intron

    def _is_around_site(self, read, start, stop):
        read_splicer = Read(read.cigar, read.reference_start)
        read_splice_sites = read_splicer.get_splice_sites()
        for site in read_splice_sites:
            if site[0] <= start and site[1] >= stop:
                return True
        return False

