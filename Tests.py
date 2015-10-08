import unittest as ut
import parsers as ps
import splice_ratio


class TestBam(ut.TestCase):

    @classmethod
    def setUpClass(cls):
        bam_fields = ['asd', '0', 'chr1', '100', '50', '100M', '*', '0', '0', 'A'*100, 'A'*100 + '\n']
        splic_fields = ['asd', '0', 'chr1', '100', '50', '50M50N50M', '*', '0', '0', 'A'*100, 'A'*100 + '\n']
        header = '\t'.join(['@SQ','SN:chr1','LN:1000000'])
        sam_data = '\n'.join([header, '\t'.join(bam_fields)])
        splic_data = '\n'.join([header, '\t'.join(splic_fields)])
        parser_sam = ps.Sam(sam_data)
        cls.parser_bam = ps.Bam(parser_sam.to_indexed_bam())
        parser_splic = ps.Sam(splic_data)
        cls.parser_splic_bam = ps.Bam(parser_splic.to_indexed_bam())
        parser_sam.delete()
        parser_splic.delete()

    @classmethod
    def tearDownClass(cls):
        cls.parser_bam.delete()
        cls.parser_splic_bam.delete()
