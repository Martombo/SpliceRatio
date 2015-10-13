import unittest as ut
import unittest.mock as um
import parsers as ps
import splice_ratio as sr


class TestParsers(ut.TestCase):
    def test_parse_site(self):
        site_str = 'chr1_18253_25983_+'
        site_dic = sr._parse_site(site_str)
        self.assertEquals('chr1',site_dic['chrom'])
        self.assertEquals(18253,site_dic['start'])
        self.assertEquals(25983,site_dic['stop'])
        self.assertEquals('+',site_dic['strand'])

    def test_parse_factor(self):
        factor = ['A', 'A', 'B', 'B']
        samples_dic = sr._parse_factor(factor)
        self.assertEquals(2, len(samples_dic))
        self.assertEquals([0,1], samples_dic['A'])
        self.assertEquals([2,3], samples_dic['B'])


class TestOverlap(ut.TestCase):
    def test_overlap(self):
        read = um.Mock(query_length = 100, reference_start = 1000)
        site = {'start': 900, 'stop':1100}
        overlap = sr.SpliceRatioCounter._overlap(read, site)
        self.assertEquals(100, overlap)

    def test_overlap_partial(self):
        read = um.Mock(query_length = 100, reference_start = 1000)
        site = {'start': 1050, 'stop':1100}
        overlap = sr.SpliceRatioCounter._overlap(read, site)
        self.assertEquals(50, overlap)

    def test_overlap_null(self):
        read = um.Mock(query_length = 100, reference_start = 10000)
        site = {'start': 1050, 'stop':1100}
        overlap = sr.SpliceRatioCounter._overlap(read, site)
        self.assertEquals(0, overlap)


class TestReadSurroundsSite(ut.TestCase):
    def test_surrounds(self):
        read = um.Mock(splice_sites=[(100, 200)])
        site = {'start': 120, 'stop':180}
        self.assertTrue(sr.SpliceRatioCounter._read_surrounds_site(read, site))

    def test_doesnt_surrounds(self):
        read = um.Mock(splice_sites=[(120, 180)])
        site = {'start': 100, 'stop':200}
        self.assertFalse(sr.SpliceRatioCounter._read_surrounds_site(read, site))

    def test_doesnt_surrounds2(self):
        read = um.Mock(splice_sites=[(80, 180)])
        site = {'start': 100, 'stop':200}
        self.assertFalse(sr.SpliceRatioCounter._read_surrounds_site(read, site))

    def test_no_splice_sites(self):
        read = um.Mock(splice_sites=[])
        site = {'start': 100, 'stop':200}
        self.assertFalse(sr.SpliceRatioCounter._read_surrounds_site(read, site))


class TestReadArrangement(ut.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.read = um.Mock(mapq=20, cigar=[], reference_start=100, query_length=80, is_reverse=False, is_read2=False)
        cls.site_str = '1_100_200_+'
        r_orientation = 'forward'
        bam = um.Mock(ps.Bam, reads_orientation=r_orientation)
        bam.determine_strand = um.Mock(return_value='+')
        cls.counter = sr.SpliceRatioCounter(reads_orientation=r_orientation)
        cls.counter.bams = [bam]
        cls.counter._read_surrounds_site = um.MagicMock(return_value=False)
        cls.counter._overlap = um.MagicMock(return_value=100)
        cls.counter.start_dic = {cls.site_str: []}
        cls.counter.last_read_len = {cls.site_str: []}

    def test_intron(self):
        self.assertEquals(self.counter._read_arrangement(self.site_str, self.read), 'intron')
        self.counter.start_dic[self.site_str] = []

    def test_surrounding(self):
        self.counter._read_surrounds_site = um.MagicMock(return_value=True)
        self.assertEquals(self.counter._read_arrangement(self.site_str, self.read), 'surrounding')
        self.counter._read_surrounds_site = um.MagicMock(return_value=False)

    def test_invalid_mapq(self):
        self.read.mapq = 2
        self.assertEquals(self.counter._read_arrangement(self.site_str, self.read), 'invalid')
        self.read.mapq = 20

    def test_invalid_strand(self):
        self.site_str = '1_100_200_-'
        self.assertEquals(self.counter._read_arrangement(self.site_str, self.read), 'invalid')
        self.site_str = '1_100_200_+'

    def test_invalid_mm(self):
        self.read.cigar = [(3, 3)]
        self.assertEquals(self.counter._read_arrangement(self.site_str, self.read), 'invalid')
        self.read.cigar = []

    def test_invalid_overlap(self):
        self.counter._overlap = um.MagicMock(return_value=5)
        self.assertEquals(self.counter._read_arrangement(self.site_str, self.read), 'invalid')
        self.counter._overlap = um.MagicMock(return_value=15)

    def test_start_dic(self):
        self.counter._read_arrangement(self.site_str, self.read)
        self.counter._read_arrangement(self.site_str, self.read)
        self.assertEquals(2, len(self.counter.start_dic[self.site_str]))
        self.assertEquals(80, self.counter.last_read_len[self.site_str])
        self.assertEquals(100, self.counter.start_dic[self.site_str][0])



class TestHighLevel(ut.TestCase):
    @classmethod
    def setUpClass(cls):
        bam_fields = ['asd', '0', 'chr1', '100', '50', '100M', '*', '0', '0', 'A' * 100, 'A' * 100 + '\n']
        header = '\t'.join(['@SQ', 'SN:chr1', 'LN:1000000'])
        sam_data = '\n'.join([header, '\t'.join(bam_fields)])
        parser_sam = ps.Sam(sam_data)
        cls.parser_bam = ps.Bam(parser_sam.to_indexed_bam())
        parser_sam.delete()
        path = cls.parser_bam.path

    @classmethod
    def tearDownClass(cls):
        cls.parser_bam.delete()
