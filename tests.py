import unittest as ut
import unittest.mock as um
import parsers as ps
import splice_ratio as sr
from copy import copy
from random import randrange
from scipy.stats import chisquare


class TestParsers(ut.TestCase):
    def test_parse_site(self):
        site_str = 'chr1_18253_25983_+'
        site_dic = sr._parse_site(site_str)
        self.assertEquals('chr1',site_dic['chrom'])
        self.assertEquals(18253,site_dic['start'])
        self.assertEquals(25983,site_dic['stop'])
        self.assertEquals('+',site_dic['strand'])

    def test_parse_site_underscore(self):
        site_str = 'chr1_gl0000asd_18263_25973_-'
        site_dic = sr._parse_site(site_str)
        self.assertEquals('chr1_gl0000asd',site_dic['chrom'])
        self.assertEquals(18263,site_dic['start'])
        self.assertEquals(25973,site_dic['stop'])
        self.assertEquals('-',site_dic['strand'])

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
        cls.counter = sr.SpliceRatioCounter(reads_orientation=r_orientation, test=True)
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


class TestAddSpliceSites(ut.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.counter = sr.SpliceRatioCounter(test=True)
        cls.counter.n_bams = 2
        cls.counter._add_splice_sites({'chr1_10_20_+': 10}, 0)
        cls.dic = {'intron': 0, 'junction': 10, 'invalid': 0, 'surrounding': 0}
        cls.empty_dic = {'intron': 0, 'junction': 0, 'invalid': 0, 'surrounding': 0}

    def test_add(self):
        self.assertEquals(1, len(self.counter.splices_dic))
        self.assertEquals(self.dic, self.counter.splices_dic['chr1_10_20_+'][0])

    def test_add_another(self):
        self.counter._add_splice_sites({'chr2_20_40_+': 20}, 0)
        self.assertEquals(2, len(self.counter.splices_dic))
        dic = copy(self.dic)
        dic['junction'] = 20
        self.assertEquals(dic, self.counter.splices_dic['chr2_20_40_+'][0])

    def test_add_bam2(self):
        self.counter._add_splice_sites({'chr3_30_50_+': 10}, 1)
        self.assertEquals(self.dic, self.counter.splices_dic['chr3_30_50_+'][1])
        self.assertEquals(self.empty_dic, self.counter.splices_dic['chr3_30_50_+'][0])


class TestGetSplicesCoverage(ut.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.dic = {'intron': 0, 'junction': 1, 'invalid': 0, 'surrounding': 0}
        header = '\t'.join(['@SQ', 'SN:chr1', 'LN:1000000'])
        bam_fields = ['asd', '0', 'chr1', '100', '50', '100M', '*', '0', '0', 'A'*100,'A'*100+'\n']
        parser_sam1 = ps.Sam('\n'.join([header, '\t'.join(bam_fields)]))
        cls.parser_bam1 = ps.Bam(parser_sam1.to_indexed_bam())
        splice_fields = ['asd', '0', 'chr1', '100', '50', '50M100N50M', '*', '0', '0', 'A'*100, 'A'*100+'\n']
        parser_sam2 = ps.Sam('\n'.join([header, '\t'.join(splice_fields)]))
        cls.parser_bam2 = ps.Bam(parser_sam2.to_indexed_bam())
        parser_sam3 = ps.Sam('\n'.join([header, '\t'.join(splice_fields) + '\t'.join(bam_fields)]))
        cls.parser_bam3 = ps.Bam(parser_sam3.to_indexed_bam())
        parser_sam1.delete()
        parser_sam2.delete()
        parser_sam3.delete()

    def test_simple(self):
        counter = sr.SpliceRatioCounter(self.parser_bam1.path)
        self.assertEquals(1, len(counter.bams))
        self.assertEquals(0, len(counter.splices_dic))
        self.assertEquals(0, len(counter.last_read_len))
        self.assertEquals(0, len(counter.start_dic))

    def test_splice(self):
        counter = sr.SpliceRatioCounter(self.parser_bam2.path)
        self.assertEquals(1, len(counter.splices_dic))
        self.assertEquals(self.dic, counter.splices_dic['chr1_149_249_+'][0])

    def test_splice_intron(self):
        counter = sr.SpliceRatioCounter(self.parser_bam3.path)
        self.assertEquals(1, len(counter.splices_dic))
        dic = copy(self.dic)
        dic['intron'] = 1
        self.assertEquals(self.dic, counter.splices_dic['chr1_149_249_+'][0])


    @classmethod
    def tearDownClass(cls):
        cls.parser_bam1.delete()
        cls.parser_bam2.delete()
        cls.parser_bam3.delete()


class TestChiTest(ut.TestCase):
    def test_chi(self):
        vec = [1,2,5,10,2,10,25,32]
        exp = sum(vec)/len(vec)
        self.assertAlmostEqual(chisquare(vec)[0], sr._chi_test(vec, exp))

    def test_chi_random(self):
        vec = [randrange(100) for _ in range(10)]
        exp = sum(vec)/len(vec)
        self.assertAlmostEqual(chisquare(vec)[0], sr._chi_test(vec, exp))


class TestFillBins(ut.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.n_bins = 10
        cls.filter = sr.SpliceRatioFilter(n_bins=cls.n_bins, min_valid=2, test=True)

    def test_simple(self):
        starts = range(1, self.filter.n_bins + 1)
        bins = self.filter._fill_bins(starts)
        self.assertEquals(self.n_bins, len(bins))
        self.assertEquals(len(starts), sum(bins))
        for k in bins:
            self.assertEquals(1, k)

    def test_example_small(self):
        starts = [235,265,285,305,315,335]
        bins = self.filter._fill_bins(starts)
        self.assertEquals(self.n_bins, len(bins))
        self.assertEquals(len(starts), sum(bins))
        for k in bins:
            self.assertGreaterEqual(k, 0)

    def test_example_big(self):
        starts = [235,265,285,305,315,335] * 10
        starts = sorted(starts)
        bins = self.filter._fill_bins(starts)
        self.assertEquals(self.n_bins, len(bins))
        self.assertEquals(len(starts), sum(bins))
        for k in bins:
            self.assertGreaterEqual(k, 0)

    def test_random(self):
        starts = sorted([randrange(10) for _ in range(1000)])
        bins = self.filter._fill_bins(starts)
        self.assertEquals(self.n_bins, len(bins))
        self.assertEquals(len(starts), sum(bins))
        for k in bins:
            self.assertGreaterEqual(k, 0)
        self.assertLess(max(bins) / min(bins), 3)

    def test_peak(self):
        starts = [10]*10
        bins = self.filter._fill_bins(starts)
        self.assertEquals(self.n_bins, len(bins))
        self.assertEquals(len(starts), sum(bins))
        self.assertEquals(len(starts), bins[0])
        for bin in bins[1:]:
            self.assertEquals(0, bin)


class TestNotAllCovered(ut.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.filter = sr.SpliceRatioFilter(min_covered=0.8, test=True)

    def test_covered(self):
        site = 'chr1_100_200_-'
        starts = [100]
        last_read_len = 100
        self.assertTrue(self.filter._all_covered(site, starts, last_read_len))

    def test_not_covered(self):
        site = 'chr1_100_200_+'
        starts = [300]
        last_read_len = 0
        self.assertFalse(self.filter._all_covered(site, starts, last_read_len))

    def test_enough_covered(self):
        site = 'chr1_100_200_+'
        starts = [109]
        last_read_len = 100
        self.assertTrue(self.filter._all_covered(site, starts, last_read_len))

    def test_not_enough_covered(self):
        site = 'chr1_100_200_+'
        starts = [111]
        last_read_len = 100
        self.assertFalse(self.filter._all_covered(site, starts, last_read_len))


class TestUnequalCov(ut.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.n_bins = 5
        cls.filter = sr.SpliceRatioFilter(min_valid=0, max_unequal=5, test=True)

    def test_uniform(self):
        self.filter._fill_bins = um.Mock(return_value=range(5))
        self.assertFalse(self.filter._unequal_cov([]))

    def test_equal_small(self):
        self.filter._fill_bins = um.Mock(return_value=[5,45,8,15,11])
        self.assertFalse(self.filter._unequal_cov([]))

    def test_unequal_small(self):
        self.filter._fill_bins = um.Mock(return_value=[3,55,8,15,11])
        self.assertTrue(self.filter._unequal_cov([]))

    def test_equal_big(self):
        self.filter._fill_bins = um.Mock(return_value=[50,680,80,150,200])
        self.assertFalse(self.filter._unequal_cov([]))

    def test_unequal_big(self):
        self.filter._fill_bins = um.Mock(return_value=[50,730,80,150,200])
        self.assertTrue(self.filter._unequal_cov([]))

    def test_peaky_small(self):
        self.filter._fill_bins = um.Mock(return_value=[1,18,4,2,3])
        self.assertTrue(self.filter._unequal_cov([]))

    def test_peaky_big(self):
        self.filter._fill_bins = um.Mock(return_value=[10,20,90,20,10])
        self.assertTrue(self.filter._unequal_cov([]))


class TestFewValidOrManyInvalid(ut.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.filter = sr.SpliceRatioFilter(sample_conditions=['A','B'], max_invalid_ratio=0.1, min_valid=20, test=True)
        cls.counts_ok = {'intron': 30, 'junction': 30, 'invalid': 0, 'surrounding': 0}
        cls.counts_invalid = {'intron': 30, 'junction': 30, 'invalid': 30, 'surrounding': 0}
        cls.counts_intron = {'intron': 5, 'junction': 30, 'invalid': 0, 'surrounding': 0}
        cls.counts_junction = {'intron': 30, 'junction': 10, 'invalid': 0, 'surrounding': 0}

    def test_simple(self):
        self.assertFalse(self.filter._few_valid_or_many_invalid([self.counts_ok]*2))

    def test_many_invalid(self):
        self.assertTrue(self.filter._few_valid_or_many_invalid([self.counts_ok, self.counts_invalid]))

    def test_few_introns_1sample(self):
        self.assertFalse(self.filter._few_valid_or_many_invalid([self.counts_ok, self.counts_intron]))

    def test_few_introns_2samples(self):
        self.assertTrue(self.filter._few_valid_or_many_invalid([self.counts_intron, self.counts_intron]))

    def test_few_junction(self):
        self.assertTrue(self.filter._few_valid_or_many_invalid([self.counts_ok, self.counts_junction]))

    def test_scn3a(self):
        filter = sr.SpliceRatioFilter(sample_conditions=['K']*2+['w']*2, max_invalid_ratio=0.2, min_valid=20, test=True)
        counts = [{'intron': 26, 'surrounding': 58, 'invalid': 3, 'junction': 50},
                  {'intron': 13, 'surrounding': 51, 'invalid': 1, 'junction': 46},
                  {'intron': 98, 'surrounding': 27, 'invalid': 2, 'junction': 21},
                  {'intron': 95, 'surrounding': 18, 'invalid': 0, 'junction': 16}]
        self.assertFalse(filter._few_valid_or_many_invalid(counts))


class TestFilterSites(ut.TestCase):
    @classmethod
    def setUp(cls):
        counter = sr.SpliceRatioCounter(test=True)
        counter.splices_dic = {'site1': []}
        counter.start_dic = {'site1': []}
        counter.last_read_len = {'site1': 0}
        cls.filter = sr.SpliceRatioFilter(sample_conditions=['A','B'], test=True)
        cls.filter._few_valid_or_many_invalid = um.Mock(return_value=False)
        cls.filter._unequal_cov = um.Mock(return_value=False)
        cls.filter._all_covered = um.Mock(return_value=True)
        cls.filter.splice_ratio_counts = counter

    def test_simple(self):
        self.filter.filter_sites()
        self.assertEquals(1, len(self.filter.licit_sites))
        self.assertEquals(0, len(self.filter.masked_sites))

    def test_few_valid(self):
        self.filter._few_valid_or_many_invalid = um.Mock(return_value=True)
        self.filter.filter_sites()
        self.assertEquals(0, len(self.filter.licit_sites))
        self.assertEquals(1, len(self.filter.masked_sites))

    def test_two_sites(self):
        self.filter.filter_sites()
        self.filter._unequal_cov = um.Mock(return_value=True)
        self.filter.filter_sites()
        self.assertEquals(1, len(self.filter.licit_sites))
        self.assertEquals(1, len(self.filter.masked_sites))


class TestTable(ut.TestCase):
    @classmethod
    def setUpClass(cls):
        counter = sr.SpliceRatioCounter(test=True)
        counter.start_dic = {'site_ok': [], 'site_bad': []}
        counter.last_read_len = {'site_ok': [], 'site_bad': []}
        arrangements = ['junction', 'intron', 'invalid', 'surrounding']
        counts_ok1 = dict(zip(arrangements, [500, 600, 0, 0]))
        counts_ok2 = dict(zip(arrangements, [400, 300, 0, 0]))
        counts_bad1 = dict(zip(arrangements, [0, 0, 700, 0]))
        counts_bad2 = dict(zip(arrangements, [0, 0, 200, 0]))
        counter.splices_dic = {'site_ok': [counts_ok1, counts_ok2], 'site_bad': [counts_bad1, counts_bad2]}
        cls.filter = sr.SpliceRatioFilter(splice_ratio_counts=counter, sample_conditions=['A','B'], test=True)
        cls.filter._fill_bins = um.Mock(return_value=[10]*5)
        cls.filter._all_covered = um.Mock(return_value=True)
        cls.filter._unequal_cov = um.Mock(return_value=False)

    def test_junction(self):
        self.filter.filter_sites()
        self.assertEquals('site_ok 500 400\n', self.filter.table('junction'))

    def test_intron(self):
        self.filter.filter_sites()
        self.assertEquals('site_ok 600 300\n', self.filter.table('intron'))

    def test_invalid(self):
        self.filter.filter_sites()
        self.assertEquals('site_bad True False False\n', self.filter.table('masked'))
