import subprocess as sp
import pysam as ps


class Bam():
    """utility functions to parse bam"""

    def __init__(self, path, reads_orientation='forward'):
        """
        sam data can be provided for testing
        splices_dic contains all splice sites (key) mapped to a list with the count of the site and the intron coverage
        :param reads_orientation: either 'forward', 'reverse' or 'mixed'
        """
        assert path
        assert path[-4:] == '.bam'
        assert os.path.isfile(path)
        if not os.path.isfile(path + '.bai'):
            p_index = sp.Popen(['samtools', 'index', path])
            p_index.communicate()
        assert os.path.isfile(path + '.bai')
        assert reads_orientation in ['forward', 'reverse', 'mixed']
        self.path = path
        self.pysam = ps.AlignmentFile(path, 'rb')
        self.reads_orientation = reads_orientation
        self.splices_dic = {}

    def get_matching(self, chrom, pos, strand='', min_qual=40):
        """
        retrieves all reads that match the reference at pos
        """
        fetch = self.pysam.fetch(chrom, pos, pos + 1)
        matching_reads = []
        for read in fetch:
            if read.mapq >= min_qual:
                if strand and strand == self.determine_strand(read):
                    reader = Read(read.cigar, read.reference_start)
                    if pos in reader.matches:
                        matching_reads.append(read)
        return matching_reads

    def get_coverage(self, chrom, start, stop, min_qual=40):
        """
        get the number of reads in region
        reads is counted even if only 1 base overlaps region
        :param chrom: str chromosome name
        :param start: int start
        :param stop: int stop
        :param min_qual: default TopHat: only uniquely mapped reads
        """
        fetch = self.pysam.fetch(chrom, start, stop)
        n_reads = 0
        for read in fetch:
            if read.mapq >= min_qual:
                n_reads += 1
        return n_reads

    def get_splice_sites(self, min_qual=40):
        """
        get splice sites counts as dictionary
        :param min_qual: 40 default (TopHat only uniquely mapped reads)
        :return dict: key: junction location, value: count. eg: {"chr1_12038_13759_+": 56, ...}
        """
        for read in self.pysam.fetch():
            if read.mapq < min_qual:
                continue
            reader = Read(read.cigar, read.reference_start, only_splicing=True)
            if reader.splice_sites:
                self._add_read_sites(read, reader.splice_sites)
        return self.splices_dic

    def _add_read_sites(self, read, read_splice_sites):
        chrom = self.pysam.getrname(read.reference_id)
        strand = self.determine_strand(read)
        for read_splice_site in read_splice_sites:
            self._add_site(read_splice_site, chrom, strand)

    def _add_site(self, splice_site, chrom, strand):
        locus_string = '_'.join([str(x) for x in [chrom, splice_site[0], splice_site[1], strand]])
        if locus_string in self.splices_dic:
            self.splices_dic[locus_string] += 1
        else:
            self.splices_dic[locus_string] = 1

    def determine_strand(self, read):
        if self.reads_orientation == 'mixed':
            return 'NA'
        strand_bool = True
        if read.is_reverse:
            strand_bool = not strand_bool
        if self.reads_orientation == 'reverse':
            strand_bool = not strand_bool
        if read.is_read2:
            strand_bool = not strand_bool
        return '+' if strand_bool else '-'

    def fetch(self, chrom, start, stop):
        return self.pysam.fetch(chrom, start, stop)

    def delete(self):
        os.remove(self.path)
        os.remove(self.path + '.bai')

class Sam:
    """
    handles Sam data or files.
    raw data should be provided only for tests.
    """

    def __init__(self, sam_data='', sam_path=''):
        assert bool(sam_data) != bool(sam_path)
        if sam_data:
            k = 1
            while os.path.isfile('tmp%d.sam' % k): k += 1
            sam_path = 'tmp%d.sam' % k
            with open(sam_path, 'w') as fin:
                fin.write(sam_data)
        self.sam_path = sam_path

    def to_indexed_bam(self):
        bam_file = self._to_bam()
        sorted_bam_path = self._sort(bam_file.filename.decode())
        os.remove(bam_file.filename.decode())
        return sorted_bam_path

    def _to_bam(self):
        tmp_bam_path = self.sam_path + '.bam'
        p_view = sp.Popen(['samtools', 'view', '-Sb', '-o', tmp_bam_path, self.sam_path])
        p_view.communicate()
        assert os.path.isfile(tmp_bam_path)
        return ps.AlignmentFile(tmp_bam_path, 'rb')

    def _sort(self, bam_path):
        sorted_bam_path = bam_path[:len(bam_path)-8]
        ps.sort(bam_path, sorted_bam_path)
        assert os.path.isfile(sorted_bam_path + '.bam')
        return sorted_bam_path + '.bam'

    def delete(self):
        os.remove(self.sam_path)


class Read:
    """
    retrieves read info about mapping, mostly splicing (junction) gaps
    :param cigar: read cigar string
    :param start: read mapping start position
    """
    def __init__(self, cigar, start, only_splicing = False):
        self.before_splice = True
        self.pos = start
        self.cigar = cigar
        self.splice_sites, self.matches, self.mm_indels, self.nonmatches = [], [], [], []
        if only_splicing:
            self._get_splicing_info()
        else:
            self._get_mapping_info()

    def _get_splicing_info(self):
        """
        computes:
        splice_sites: the donor and acceptor sites (junction) of a read
        list of donor and acceptor position: [(152683, 153107), (153194, 153867)]
                 None if read overlaps no junction
        """
        for cigar_token in self.cigar:
            if cigar_token[0] == 0:
                self._match(cigar_token[1])
            elif cigar_token[0] == 1:
                pass
            elif cigar_token[0] == 2:
                self._move_pos(cigar_token[1])
            elif cigar_token[0] == 3:
                self._non_match(cigar_token[1])
            else:
                break

    def _get_mapping_info(self):
        """
        computes:
        1) splice_sites: the donor and acceptor sites (junction) of a read
        list of donor and acceptor position: [(152683, 153107), (153194, 153867)]
                 None if read overlaps no junction
        2) mm_indels: the list of base positions with mismatches or indels
        3) matches: the list of base positions with matches
        4) nonmatches: the list of base positions with nonmatches
        """
        for cigar_token in self.cigar:
            if cigar_token[0] == 0:
                for k in range(cigar_token[1]):
                    self.matches.append(self.pos + k)
                self._match(cigar_token[1])
            elif cigar_token[0] == 1:
                for k in range(cigar_token[1]):
                    self.mm_indels.append(self.pos + k)
            elif cigar_token[0] == 2:
                for k in range(cigar_token[1]):
                    self.mm_indels.append(self.pos + k)
                self._move_pos(cigar_token[1])
            elif cigar_token[0] == 3:
                for k in range(cigar_token[1]):
                    self.mm_indels.append(self.pos + k)
                    self.nonmatches.append(self.pos + k)
                self._non_match(cigar_token[1])
            else:
                break

    def _match(self, leng):
        if not self.before_splice:
            self.splice_sites.append((self.start_site,self.pos))
            self.before_splice = True
        self._move_pos(leng)

    def _non_match(self, leng):
        if self.before_splice:
            self.before_splice = False
            self.start_site = self.pos
        self._move_pos(leng)

    def _move_pos(self, leng):
        self.pos += leng
