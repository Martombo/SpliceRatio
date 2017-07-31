import splice_ratio as sr
import subprocess as sp

# string of bam files, space separated
bam_paths = 'A1.bam A2.bam A3.bam B1.bam B2.bam B3.bam'

# list of conditions, in the same order as bam_paths 
condition = ['KO']*3 + ['wt']*3

# counting the reads across and within splice junctions, with options (default):
# reads_orientation ('forward'): can be either 'forward', 'reverse', 'mixed'
# min_qual (10): removes low quality reads (multi-mappers according to TopHat and similar)
# min_overlap (10): min overlap between read and intron to count it as intronic
counter = sr.SpliceRatioCounter(bam_paths, reads_orientation='reverse')
counter.get_splices_coverage()

# removing invalid splice sites according to these main options (default):
# min_valid_sample (5): min valid counts per condition
# max_invalid_ratio (0.2): max invalid counts / total
# min_covered (0.8): min % of intron covered with reads, samples of same condition aggregated
# max_invalid (20): max invalid counts per condition
weeder = sr.SpliceRatioFilter(counter, condition, min_valid_condition=10, max_invalid_ratio=0.1)
weeder.filter_sites()
for k in ['intron','junction','masked']:
    counts = weeder.table(k)
    with open(k + '_counts', 'w') as fout:
        fout.write(bam_paths + '\n')
        fout.write(counts)

# differential splice ratio analysis
p_rscript = sp.Popen(['Rscript', 'deseq2_normFactors.R'])
p_rscript.communicate()
