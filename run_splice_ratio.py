import splice_ratio as sr

counter = sr.SpliceRatioCounter(bam_paths='', reads_orientation='', min_qual=10, max_mm_indels=0, min_overlap=0.1)
counter.get_splices_coverage()
filter = sr.SpliceRatioFilter(['a','b'], counter)