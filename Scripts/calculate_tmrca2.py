from __future__ import division
import sys, os, random, itertools, shutil, cStringIO
#from numpy import array, arange
sys.path.append('/net/akey/vol1/home/bvernot/tishkoff/mapped_snps/')
sys.path.append('/net/akey/vol1/home/bvernot/tishkoff/metrics/')
sys.path.append('/net/akey/vol1/home/bvernot/tishkoff/filter_files/')
sys.path.append('/net/akey/vol1/home/bvernot/archaic_exome/experiments/fdr_simulated_basic/latest/')
#from region_stats import region_type_stats
from myBedTools3 import myBedTools
from BaseLookup import BaseLookup
import fileinput
from operator import itemgetter
import argparse
#from numpy.random import binomial
#import tables
import re
import sqlite3, time
from bitarray import bitarray
from collections import Counter, defaultdict
import get_pct_arc_per_ind_from_ms_file_new as pct_arc
from mydefaultdict import mydefaultdict
from test_parse_tree import read_tree4, find_node, tree_to_str, get_terminals, get_path_to_root, get_dist_btwn_nodes

import time
start_time = time.time()
debug_ms = False

import locale
locale.setlocale(locale.LC_ALL, 'en_US')

class InvertBinaryBedFileAction(argparse.Action):
    def __call__(self, parser, namespace, filename, option_string=None):
        bt = myBedTools(myBedTools.binarybedfilegenome, invert = True, initialize = False, ref_version = namespace.ref_version)
        bt.read_as_regions(filename, invert = True)
        setattr(namespace, self.dest, bt)
        pass
    pass

class VCFFileAction(argparse.Action):
    def __call__(self, parser, namespace, filename, option_string=None):
        vcffile = open(filename, 'r')
        vcf = defaultdict(lambda : ['N'])
        sys.stderr.write("Reading VCF file %s..\n" % filename)
        c = 0
        for line in vcffile:
            if line.strip().startswith('#'): continue
            [chrom, pos, _, ref, alt] = line.strip().split()[:5]
            if chrom not in vcf:
                vcf[chrom] = defaultdict(lambda : ['N'])
                pass
            if int(pos) in vcf[chrom]:
                print "error - duplicate position in VCF file?"
                print chrom, pos, ref, alt
                print line
                sys.exit(-1)
                pass
            vcf[chrom][int(pos)] = alt.split(',')
            c += 1
            pass
        sys.stderr.write(" with %d lines.\n" % c)
        # def vcf_class(object):
        #     vcf = vcf
        #     def get_base_one_based(self, chrom, pos):
        #         return self.vcf[chrom][pos]
        #     pass
        setattr(namespace, self.dest, vcf)
        pass
    pass

class BinaryBedFileAction(argparse.Action):
    def __call__(self, parser, namespace, filename, option_string=None):
        bt = myBedTools(myBedTools.binarybedfilegenome, initialize = False, ref_version = namespace.ref_version)
        bt.read_as_regions(filename)
        setattr(namespace, self.dest, bt)
        pass
    pass

class BedFileAction(argparse.Action):
    def __call__(self, parser, namespace, filename, option_string=None):
        bt = myBedTools(myBedTools.bedfile, initialize = True, ref_version = namespace.ref_version)
        bt.read_as_regions(filename, file_type = myBedTools.bedfile)
        setattr(namespace, self.dest, bt)
        pass
    pass

class MultBinaryBedFileAction(argparse.Action):
    def __call__(self, parser, namespace, filenames, option_string=None):
        def fn_outer(i):
            bt = myBedTools(myBedTools.binarybedfilegenome, initialize = False, ref_version = namespace.ref_version)
            bt.read_as_regions(i)
            return bt
        setattr(namespace, self.dest, [fn_outer(i) for i in filenames])
        pass
    pass

class BinarySeqFileAction(argparse.Action):
    def __call__(self, parser, namespace, filename, option_string=None):
        bt = myBedTools(myBedTools.binaryseqfilegenome, initialize = False, ref_version = namespace.ref_version)
        bt.read_as_regions(filename, file_type=myBedTools.binaryseqfilegenome)
        setattr(namespace, self.dest, bt)
        pass
    pass

class MultBinarySeqFileAction(argparse.Action):
    def __call__(self, parser, namespace, filenames, option_string=None):
        def fn_outer(i):
            bt = myBedTools(myBedTools.binaryseqfilegenome, initialize = False, ref_version = namespace.ref_version)
            bt.read_as_regions(i, file_type=myBedTools.binaryseqfilegenome)
            return bt
        setattr(namespace, self.dest, [fn_outer(i) for i in filenames])
        pass
    pass

class MergeBinaryBedFilesAction(argparse.Action):
    def __call__(self, parser, namespace, filenames, option_string=None):
        sys.stderr.write("Merging %d bbg files\n" % len(filenames))
        bt = myBedTools(myBedTools.binarybedfilegenome, initialize = False, ref_version = namespace.ref_version)
        for f in filenames:
            bt.read_as_regions(f)
            pass
        setattr(namespace, self.dest, bt)
        pass
    pass

class IntersectBinaryBedFilesAction(argparse.Action):
    def __call__(self, parser, namespace, filenames, option_string=None):
        sys.stderr.write("intersecting %d bbg files\n" % len(filenames))
        bt = myBedTools(myBedTools.binarybedfilegenome, initialize = False, ref_version = namespace.ref_version)
        bt.read_to_bases(myBedTools.binarybedfilegenome, filenames[0], myBedTools.set_to_one)
        for f in filenames[1:]:
            bt.read_to_bases(myBedTools.binarybedfilegenome, f, myBedTools.bitfn_and)
            pass
        setattr(namespace, self.dest, bt)
        pass
    pass


parser = argparse.ArgumentParser(description='Get counts, frequencies, etc.')
## ref has to be set before any bed file options
## I achieve this requirement by having it be part of parser, and not a subparser
parser.add_argument('-ref', '--ref-version', type = str, default = 'b37')


subparsers = parser.add_subparsers(dest='func')

# build parsers for various arguments
mult_vars_parser = argparse.ArgumentParser(add_help = False)
mult_vars_parser.add_argument('-v', '--variant_file', nargs='+', required=True, type=argparse.FileType('r'), dest='varfiles')

vars_parser = argparse.ArgumentParser(add_help = False)
vars_parser_group = vars_parser.add_mutually_exclusive_group(required=True)
vars_parser_group.add_argument('-v', '--variant-file', type=argparse.FileType('r'), dest='varfile')
vars_parser_group.add_argument('-ms', '--ms-file', type=argparse.FileType('r'))
vars_parser.add_argument('-ms-pops', required=False, default = None, nargs='+', type = int)
vars_parser.add_argument('-ms-am', '--ms-pops-allow-mismatch', required=False, action='store_false', dest='ms_pops_force_match')
vars_parser.add_argument('-ms-reglen', required=False, default = None, type = int)
vars_parser.add_argument('-ms-no-randomize', action='store_true')
vars_parser.add_argument('-ms-forward', '--forward-ms-file', action='store_true', help="print ms file as it's being read - this is useful for piping to another program (maybe another calc_tmrca?).  You will want to suppress any other output, though, natch.")
vars_parser.add_argument('-ms-random-seed', default = None, type = int)
vars_parser.add_argument('-forward-ms-chr-sstar', '--forward-ms-chr-sstar', action='store_true', help="This option is only used for the very strange circumstance of printing the S* value for every ms chromosome, so that they can be later interpreted by get_pct_arc_per_ind_from_ms_file_new.py to annotate results with the S* threshold required to produce those results.")
vars_parser.add_argument('-filt-retain-arc', action='store_true')
vars_parser.add_argument('-filt-retain-known-anc', action='store_true')
vars_parser.add_argument('-ms-outgroup', required=False, default = None, nargs=2, metavar=['outgroup_chrnum', 'outgroup_dist'])
vars_parser.add_argument('-ms-reglens', '--ms-region-length-distribution', type = argparse.FileType('r'), required = False, default = None, help = 'Select the length of the region from this file.  Adjust the probability that a given base will be sampled by len/window_length.')
vars_parser.add_argument('-ms-probn', '--base_sample_probability_normal', type=float, nargs=2, required = False, default = None, help = 'Select the length of the region from this normal distribution (given by mean and stdev).  Adjust the probability that a given base will be sampled by len/window_length.')
vars_parser.add_argument('-prob', '--base_sample_probability', type=float, required = False, default = 1, help = 'The probability that a given base will be sampled.  This should not be used with base_sample_probability_normal.')
vars_parser.add_argument('-1kg', '--varfile-is-1kg', action='store_true')
vars_parser.add_argument('-1kg-ids', '--use-1kg-id', action='store_true')
vars_parser.add_argument('-esp', '--varfile-is-esp', action='store_true')
vars_parser.add_argument('-tair10', '--varfile-is-tair10', action='store_true')
vars_parser.add_argument('-prog', '--report-progress', default = None, type = int, help='Report progress every x snps.')
vars_parser.add_argument('-sim-prefix', '--sim-prefix', default = '', help='Prefix for simulated chromosomes; only used with -ms.')
vars_parser.add_argument('-allow-fixed', '--allow-fixed', action='store_false', dest='remove_fixed_vars')

base_parser = argparse.ArgumentParser(add_help = False)
base_parser.add_argument('-d', '--debug', action = 'store_true')
base_parser.add_argument('-dd', '--debug-2', action = 'store_true')
base_parser.add_argument('-ddd', '--debug-3', action = 'store_true')
base_parser.add_argument('-dsnps', '--debug-read-snps', action = 'store_true')
base_parser.add_argument('-dsstar', '--debug-sstar', action = 'store_true')
base_parser.add_argument('-kpc', '--keep-partial-calls', action="store_false", dest='fpc')
base_parser.add_argument('-fm', '--filter-missingness', default = .8, type = float, help='-fm 1 means that no sites with missing data are considered; -fm 0 allows all sites.')
base_parser.add_argument('-minf', '--min-freq', default = 0, type = float)
base_parser.add_argument('-maxf', '--max-freq', default = 1, type = float)
base_parser.add_argument('-kinds', '--keep-inds', nargs='+', help='keep these individuals (give numbers, not names - one based)', type = int)
base_parser.add_argument('-kinds-id', '--keep-inds-by-id', nargs='+', help='keep these individuals (give names - can force this to be Individual IDs for 1kg using --use-1kg-id)')
base_parser.add_argument('-ainds', '--add-inds', nargs='+', help='add these individuals (give numbers, not names - one based)', type = int)
base_parser.add_argument('-ainds-id', '--add-inds-by-id', nargs='+', help='add these individuals (give names - can force this to be Individual IDs for 1kg using --use-1kg-id)')
#base_parser.add_argument('-indsf', '--filter-inds-file', type = open, help='keep only these individuals (in a file)')
base_parser.add_argument('-pops', '--keep-pops', nargs='+', help='keep only these populations.  If calculating a two-population statistic, use -pops1 and -pops2')
base_parser.add_argument('-rpops', '--remove-pops', nargs='+', help='remove these populations.  If calculating a two-population statistic, use -pops1 and -pops2')
base_parser.add_argument('-r', '--regions', action=BinaryBedFileAction, required=False, default=None, help = 'A bbg file that specifies which regions to consider.  Only snps in this region are loaded.')
base_parser.add_argument('-rb', '--regions-bed', action=BedFileAction, required=False, default=None, help = 'A bbg file that specifies which regions to consider.  Only snps in this region are loaded.', dest = 'regions')
#base_parser.add_argument('-ir', '--intersect-region', action=MergeBinaryBedFilesAction, nargs='+', required=False, default=None, help = 'A bbg file that is intersected with the --regions file to produce a new set of regions to consider.')
base_parser.add_argument('-ir', '--intersect-region', action=IntersectBinaryBedFilesAction, nargs='+', required=False, default=None, help = 'A bbg file that is intersected with the --regions file to produce a new set of regions to consider.')
base_parser.add_argument('-x', '--exclude-region', nargs='+', action=MergeBinaryBedFilesAction, required=False, default=None, help = 'bbg file(s) that specify which regions should be excluded from the analysis.  If more than one file is given, the files are merged.')
base_parser.add_argument('-o', '--output-file', type=argparse.FileType('w'), required = False, default = sys.stdout, help = 'output file')
base_parser.add_argument('-pad', '--pad-regions', required=False, type=int, help='Pad regions by this amount on both sides.  This happens before any other regions operations (intersect/exclude).')
base_parser.add_argument('-nc100', '--exclude-all-nocalls', nargs='?', const='/net/akey/vol1/home/bvernot/tishkoff/cg_coverage/all_inds/var-IND.tsv/var-IND.tsv.nocalls_norefs.merged.bed.bbg')
base_parser.add_argument('-wfrac', '--required-window-fraction', required=False, default = .5, type=float, help='Require this fraction of the window\'s sites to pass filter (this is applied to chimp mapped and to human sequence).')
base_parser.add_argument('-data-name', '--data-name', default = 'whole_file', required=False)
base_parser.add_argument('-db', '--data-base-name', nargs = 2, default = None)
base_parser.add_argument('-tag', '--tag', nargs = 2, action = 'append', default = None)
base_parser.add_argument('-lim', '--limit', type = int, default = None)
base_parser.add_argument('-ilim', '--individual-limit', type = int, default = None)
base_parser.add_argument('-head', '--print-header', action = 'store_true')
base_parser.add_argument('-revisions', '--print-sstar-data-for-revisions-on-neand-paper', action = 'store_true')



selection_parser = argparse.ArgumentParser(add_help = False)
selection_parser.add_argument('-neu', '--neutral-regions', nargs='+', action = MergeBinaryBedFilesAction, required = True, help = 'One or more bbg files that denote putatively neutral regions.  If more than one file is given, all regions are merged.')
selection_parser.add_argument('-neux', '--neutral-regions-exclude', nargs='+', action = MergeBinaryBedFilesAction, required = False, help = 'One or more bbg files that should be excluded from the putatively neutral regions.  If more than one file is given, all files are merged.', default = None)
selection_parser.add_argument('-sel', '--selected-regions', nargs='+', action = MergeBinaryBedFilesAction, required = True, help = 'One or more bbg files that denote putatively selected regions.  If more than one file is given, all regions are merged')

multiple_population_parser = argparse.ArgumentParser(add_help = False)
multiple_population_parser.add_argument('-mpops', '--multiple-populations', nargs='+', required = False, help='list of populations, separated by commas')

# two_population_parser = argparse.ArgumentParser(add_help = False)
# two_population_parser.add_argument('-pops1', '--population1', nargs='+', required = True, help='In a two-population statistic, use these populations for population 1.')
# two_population_parser.add_argument('-pops2', '--population2', nargs='+', required = True, help='In a two-population statistic, use these populations for population 2.')

extra_population_parser = argparse.ArgumentParser(add_help = False)
extra_population_parser.add_argument('-tp', '--target-pops', nargs='+', required = True, help='In a statistic where a sub-population is targeted for analysis, use these populations.')

second_region_parser = argparse.ArgumentParser(add_help = False)
second_region_parser.add_argument('-sr', '--second-region', action=BinaryBedFileAction, required=True, help = 'A bbg file that specifies a second region to consider - often this is a subset of the first.  e.g., for region-ratio.')

chimp_snps_parser = argparse.ArgumentParser(add_help = False)
chimp_snps_parser.add_argument('-cs', '--chimp-snps', type=argparse.FileType('r'), required = False, help = 'chimp snps file')
# default = '/net/akey/vol1/home/bvernot/tishkoff/primate_sequences/alignments/hg19.panTro2.synNet.axt.5.snps', 

chimp_unmapped_parser = argparse.ArgumentParser(add_help = False)
chimp_unmapped_parser.add_argument('-cu', '--chimp-unmapped', action = InvertBinaryBedFileAction, required = False, help = 'chimp unmapped bbg file', dest = 'chimp_mapped')
# default = '/net/akey/vol1/home/bvernot/tishkoff/primate_sequences/alignments/hg19.panTro2.synNet.axt.5.unmapped.zerobased.bbg',

compare_calls_parser = argparse.ArgumentParser(add_help = False)
#compare_calls_parser.add_argument('-comp', '--compare-calls', action = BinarySeqFileAction, required = False, help = 'Calls for comparison with model prediction.')
compare_calls_parser.add_argument('-comp', '--compare-calls', action = MultBinarySeqFileAction, nargs = '+', required = False, help = 'Multiple BSG files.')
compare_calls_parser.add_argument('-compvcf', '--compare-calls-vcf', action = VCFFileAction, required = False, help = 'Multiple VCF files (I think only a single file works now..).')

mult_bsg_parser = argparse.ArgumentParser(add_help = False)
mult_bsg_parser.add_argument('-mbsg', '--bsg-files', action = MultBinarySeqFileAction, nargs = '+', required = False, help = 'Multiple BSG files.')

subset_parser = argparse.ArgumentParser(add_help = False)
subset_parser.add_argument('-subset', '--subset-inds-into-groups', type = int, default = None, required = False, help = 'Split the population into groups of N individuals, and run the analysis on each group separately.')
subset_parser.add_argument('-subset-not-random', '--subset-not-random', action = 'store_true')

output_option_parser = argparse.ArgumentParser(add_help = False)
output_option_parser.add_argument('-output', '--output-format', default = 'win', choices = ['win', 'snps', 'snpbed'], help = 'Output format')

ancestry_parser = argparse.ArgumentParser(add_help = False)
ancestry_parser_group = ancestry_parser.add_mutually_exclusive_group(required = False)
ancestry_parser_group.add_argument('-anc', '--ancestry_file', required = False, type=argparse.FileType('r'), default = None)
ancestry_parser.add_argument('-anc-col', '--ancestral_column', required = False, type=int, default=8)
ancestry_parser_group.add_argument('-anc-bsg', '--ancestry_file_bsg', required = False, action=BinarySeqFileAction, default = None)


time_to_primate_divergence_parser = argparse.ArgumentParser(add_help = False)
time_to_primate_divergence_parser.add_argument('-pd', '--time-to-primate-divergence', nargs=1, required=False, type=int, default=6000000)

window_file_parser = argparse.ArgumentParser(add_help = False)
window_file_parser_group = window_file_parser.add_mutually_exclusive_group(required=True)
window_file_parser_group.add_argument('-wf', '--window-file', type=argparse.FileType('r'))
window_file_parser_group.add_argument('-ms-win', '--ms-win-params', nargs=2, type=int, metavar=['winlen', 'winstep'])
window_file_parser.add_argument('-win4', '--merge-windows_by-fourth-column', action = 'store_true')

optional_window_file_parser = argparse.ArgumentParser(add_help = False)
optional_window_file_parser_group = optional_window_file_parser.add_mutually_exclusive_group(required=False)
optional_window_file_parser_group.add_argument('-wf', '--window-file', required=False, type=argparse.FileType('r'), default = None)
optional_window_file_parser_group.add_argument('-ms-win', '--ms-win-params', nargs=2, type=int, metavar=['winlen', 'winstep'])
optional_window_file_parser.add_argument('-win4', '--merge-windows_by-fourth-column', action = 'store_true')
#optional_window_file_parser.add_argument('-wbbg', '--window-bbg-file', required=False, action=BinaryBedFileAction, default = None)

window_file_output_parser = argparse.ArgumentParser(add_help = False)
window_file_output_parser.add_argument('-wf-odir', '--window-file-output-dir', required = False, help="Create an output directory structure, where each window gets its own output file - useful if you're going to be reading in snps for each window (i.e., to make a plot for a region).  I think this is only used for S*..", default = None)
window_file_output_parser.add_argument('-wf-append', '--window-file-append', required = False, help="Append window output to this file - the headers will be messed up, especially if snps are also output, but it makes for a less complicated directory structure.", default = None, type=argparse.FileType('w'))
window_file_output_parser.add_argument('-wf-odir-tmp', '--window-file-tmp-output-dir', required = False, default = None)
window_file_output_parser.add_argument('-wf-odir-tmp-buffer', '--window-file-tmp-output-dir-buffer', required = False, default = 1)


summary_only_option_parser = argparse.ArgumentParser(add_help = False)
summary_only_option_parser.add_argument('-sum', '--summary-only', required=False, action='store_true')
summary_only_option_parser.add_argument('-no-sum', '--suppress-summary', required=False, action='store_true')

snp_report_option_parser = argparse.ArgumentParser(add_help = False)
snp_report_option_parser.add_argument('-no-snp-report', '--no-snp-report', required=False, action='store_true')

folded_sfs_option_parser = argparse.ArgumentParser(add_help = False)
folded_sfs_option_parser.add_argument('-folded-sfs', '--folded-sfs', required=False, action='store_true')

calc_pearson_option_parser = argparse.ArgumentParser(add_help = False)
calc_pearson_option_parser.add_argument('-no-pearson', '--no-pearson', required=False, action='store_false', default=True, dest='calc_pearson', help = 'By default we do a pearson calculation - use this argument to turn that off (runs much faster, but doesn\'t give you linked vars).')

merge_windows_option_parser = argparse.ArgumentParser(add_help = False)
merge_windows_option_parser.add_argument('-merge-windows', '--merge-windows', required=False, action='store_true')

s_star_args_parser = argparse.ArgumentParser(add_help = False)
s_star_args_parser.add_argument('-allowable-diffs', default = 5, type=int)
s_star_args_parser.add_argument('-max-perf', '--max-diffs-for-perfect-score', default = 0, type=int)
s_star_args_parser.add_argument('-impch', '--imperfect-charge', default = -10000, type=int)
s_star_args_parser.add_argument('-perf-bonus', '--perfect-score-bonus', default = 5000, type=int)
s_star_args_parser.add_argument('-exclude-singletons', action = 'store_true')
s_star_args_parser.add_argument('-num-snps', '--pop-specific-snp-count', default = None, type = int)

report_snps_option_parser = argparse.ArgumentParser(add_help = False)
report_snps_option_parser.add_argument('-rs', '--report-snps', required=False, action='store_true')

is_male_file_parser = argparse.ArgumentParser(add_help = False)
is_male_file_parser.add_argument('-sex', '--is-male-file', nargs=1, required=False, default='/net/akey/vol1/home/bvernot/tishkoff/cg_coverage/all_inds/var-IND.tsv/var-IND.tsv.nocalls_norefs.is_male', dest='is_male')

resample_option_parser = argparse.ArgumentParser(add_help = False)
resample_option_parser.add_argument('-resample', '--resample', required=False, action='store_true')
resample_option_parser.add_argument('-print-resamples', '--print-resamples', required=False, action = 'store_true')

s4_option_parser = argparse.ArgumentParser(add_help = False)
s4_option_parser.add_argument('-s4', action = 'store_true', help = 'Get window name from fourth column.')

method_number_parser = argparse.ArgumentParser(add_help = False)
method_number_parser.add_argument('-m', '--method-number', default=['0'], nargs='+', help = 'Select a method number (for methods with more than one base option - this is mostly used for experimental methods).')

# _parser = argparse.ArgumentParser(add_help = False)
# _parser.add_argument('-pd', '--time-to-primate-divergence', nargs=1, required=True, type=int, dest=ancestry_file)
# _parser = argparse.ArgumentParser(add_help = False)
# _parser.add_argument('-pd', '--time-to-primate-divergence', nargs=1, required=True, type=int, dest=ancestry_file)

test_anc_req = argparse.ArgumentParser(add_help = False)
test_anc_req.add_argument('-anc', '--ancestry_req', required=False, help=argparse.SUPPRESS, default = 'hey')
test_anc_req.add_argument('-anc2', '--ancestry_req2', required=False, default = 'hey')

individual_regions_parser = argparse.ArgumentParser(add_help = False)
individual_regions_parser.add_argument('-ind-reg-names', '--individual-regions-names', nargs='+', required=False, default = None)
individual_regions_parser.add_argument('-ind-regs', '--individual-regions', required=False, default = None, action = MultBinaryBedFileAction, nargs='+')

# test parsing settings
test_parser = subparsers.add_parser('test', parents = [ancestry_parser, base_parser, is_male_file_parser, summary_only_option_parser, optional_window_file_parser, multiple_population_parser])

# freqs arguments
freqs_parser = subparsers.add_parser('freqs', parents = [vars_parser, base_parser, is_male_file_parser, summary_only_option_parser, optional_window_file_parser, multiple_population_parser, ancestry_parser, compare_calls_parser])

# sites arguments
sites_parser = subparsers.add_parser('sites', parents = [vars_parser, base_parser, is_male_file_parser, multiple_population_parser])

# summary statistics arguments
window_file_output_pstats_parser = subparsers.add_parser('stats', parents = [vars_parser, base_parser, is_male_file_parser, summary_only_option_parser, optional_window_file_parser, multiple_population_parser, ancestry_parser, compare_calls_parser])
windowed_stats_parser = subparsers.add_parser('stats-win', parents = [vars_parser, base_parser, is_male_file_parser, summary_only_option_parser, window_file_parser, multiple_population_parser, ancestry_parser, compare_calls_parser, s_star_args_parser, extra_population_parser])

# output ms file arguments
msfile_parser = subparsers.add_parser('msfile', parents = [vars_parser, ancestry_parser, base_parser, is_male_file_parser, optional_window_file_parser])

# runs of homozygosity arguments
roh_parser = subparsers.add_parser('roh', parents = [vars_parser, base_parser, is_male_file_parser, summary_only_option_parser, optional_window_file_parser, multiple_population_parser])

# snp counts arguments
snpcount_parser = subparsers.add_parser('snpcount', parents = [vars_parser, base_parser, is_male_file_parser, optional_window_file_parser])

# snp counts arguments
site_freq_parser = subparsers.add_parser('site_freq', parents = [vars_parser, base_parser, is_male_file_parser, optional_window_file_parser, summary_only_option_parser, folded_sfs_option_parser])

# pi arguments
pi_parser = subparsers.add_parser('pi', parents = [vars_parser, base_parser, optional_window_file_parser, is_male_file_parser, summary_only_option_parser, resample_option_parser])

# dist matrix arguments
dist_matrix_parser = subparsers.add_parser('dist_matrix', parents = [vars_parser, base_parser, optional_window_file_parser, is_male_file_parser, summary_only_option_parser, resample_option_parser, ancestry_parser])

# mk arguments
mk_parser = subparsers.add_parser('mk', parents = [vars_parser, base_parser, optional_window_file_parser, is_male_file_parser, resample_option_parser, chimp_unmapped_parser, chimp_snps_parser, s4_option_parser, selection_parser])

# fst arguments
fst_parser = subparsers.add_parser('fst', parents = [vars_parser, base_parser, optional_window_file_parser, is_male_file_parser, summary_only_option_parser, resample_option_parser, s4_option_parser, multiple_population_parser])

# mutation_rate arguments
mut_rate_parser = subparsers.add_parser('mutrate', parents = [base_parser, optional_window_file_parser, chimp_unmapped_parser, chimp_snps_parser, summary_only_option_parser, resample_option_parser, s4_option_parser])

# region length arguments
#region_length_parser = subparsers.add_parser('reglen', parents = [base_parser, optional_window_file_parser, summary_only_option_parser, s4_option_parser, vars_parser])
region_length_parser = subparsers.add_parser('reglen', parents = [base_parser, optional_window_file_parser, summary_only_option_parser, s4_option_parser])

# region ratio arguments
two_region_ratio_parser = subparsers.add_parser('region-ratio', parents = [base_parser, optional_window_file_parser, summary_only_option_parser, s4_option_parser, second_region_parser])

# counts arguments
counts_parser = subparsers.add_parser('counts', parents = [mult_vars_parser, base_parser])

# plink arguments
plink_parser = subparsers.add_parser('plink', parents = [vars_parser, base_parser, is_male_file_parser, individual_regions_parser])

# PED/MAP arguments
ped_map_parser = subparsers.add_parser('ped_map', parents = [vars_parser, base_parser, is_male_file_parser, merge_windows_option_parser, optional_window_file_parser])

# visual genome arguments
visual_genotype_parser = subparsers.add_parser('visual_genotype', parents = [vars_parser, base_parser, is_male_file_parser, second_region_parser, optional_window_file_parser])

# visual fasta arguments
visual_fasta_parser = subparsers.add_parser('visual_fasta', parents = [vars_parser, base_parser, is_male_file_parser, second_region_parser, optional_window_file_parser, mult_bsg_parser])

# tmrca arguments
tmrca_parser = subparsers.add_parser('tmrca', parents = [base_parser,
                                                         vars_parser,
                                                         chimp_snps_parser,
                                                         chimp_unmapped_parser,
                                                         ancestry_parser,
                                                         time_to_primate_divergence_parser,
                                                         window_file_parser,
                                                         is_male_file_parser])


# tmrca subset arguments
tmrca_subset_parser = subparsers.add_parser('tmrca_subsets', parents = [base_parser,
                                                                        vars_parser,
                                                                        chimp_snps_parser,
                                                                        chimp_unmapped_parser,
                                                                        ancestry_parser,
                                                                        time_to_primate_divergence_parser,
                                                                        window_file_parser,
                                                                        method_number_parser,
                                                                        extra_population_parser,
                                                                        compare_calls_parser,
                                                                        output_option_parser,
                                                                        is_male_file_parser,
                                                                        report_snps_option_parser])

# region_s_star arguments
region_s_star_parser = subparsers.add_parser('region_s_star', parents = [base_parser,
                                                                         vars_parser,
                                                                         window_file_parser,
                                                                         method_number_parser,
                                                                         extra_population_parser,
                                                                         compare_calls_parser,
                                                                         output_option_parser,
                                                                         is_male_file_parser,
                                                                         report_snps_option_parser,
                                                                         s_star_args_parser,
                                                                         window_file_output_parser,
                                                                         summary_only_option_parser,
                                                                         calc_pearson_option_parser,
                                                                         snp_report_option_parser,
                                                                         subset_parser])


args = parser.parse_args()

# print 'argssss', args

## set debug arguments
if args.debug_3: args.debug_2 = True
if args.debug_2: args.debug = True

if args.debug: print args

setattr(args, 'check_time', time.time())
setattr(args, 'first_line', True)

setattr(args, 'baselookup', BaseLookup('.', ref_version = args.ref_version))

if args.data_base_name != None:
    setattr(args, 'table_name', args.data_base_name[1])
    args.data_base_name = args.data_base_name[0]
    pass

if 'base_sample_probability' not in args:
    setattr(args, 'base_sample_probability', 1)
    pass

if args.tag != None:
    tag_dict = {}
    for (k,v) in args.tag:
        tag_dict[k] = v
        pass
    args.tag = tag_dict
else:
    setattr(args, 'tag', {})
    pass

if 'compare_calls' in args or 'compare_calls_vcf' in args:
    setattr(args, 'comp', True)
else:
    setattr(args, 'comp', False)
    pass

if 'forward_ms_file' in args:
    ms_file_backlog = defaultdict(lambda : cStringIO.StringIO())
    pass

## check arguments
# if we require an ancestry file, only actually require it with a varfile
if 'varfile' in args or 'ms_file' in args:

    if 'ancestry_file' in args and args.ancestry_file == None and args.ancestry_file_bsg == None and 'varfile' in args and args.varfile != None:
        print "Ancestry File required: -anc"
        sys.exit(-1)
        pass

    if 'ms_file' in args and args.ms_file != None and args.ms_reglen == None:
        print "-ms-reglen is required with -ms-file"
        sys.exit(-1)
        pass

    if args.ms_file == None and 'chimp_snps' in args and args.chimp_snps == None:
        print "--chimp-snps is required."
        sys.exit(-1)
        pass
    
    if args.ms_file == None and 'chimp_unmapped' in args and args.chimp_mapped == None:
        print "--chimp-unmapped is required."
        sys.exit(-1)
        pass
        
    if args.ms_file == None and 'ms_win_params' in args and args.ms_win_params != None:
        print "--window-file is required with -v (--ms-win-params is only valid with -ms-file)."
        sys.exit(-1)
        pass
    
    if args.ms_region_length_distribution != None:
        reglens = []
        for line in args.ms_region_length_distribution:
            reglens.append(int(line.strip()))
            pass
        args.ms_region_length_distribution = reglens
        pass
    
    pass

## set the flag for ms or not ms
if 'ms_file' in args and args.ms_file != None:
    setattr(args, 'ms', True)
else:
    setattr(args, 'ms', False)
    pass

#sys.exit()

## set up variables?
if 'second_region' not in args: setattr(args, 'second_region', None)
    
if 'multiple_populations' not in args:
    setattr(args, 'multiple_populations', None)
    pass

# if 'population1' not in args: setattr(args, 'population1', None)
# if 'population2' not in args: setattr(args, 'population2', None)

# if args.multiple_populations != None and len(multiple_populations) == 2 and args.population1 == None and args.population2 == None:
#     args.population1 = args.multiple_populations[0]
#     args.population2 = args.multiple_populations[1]
#     pass

# if args.population1 != None and args.population2 != None:
#     setattr(args, 'two_populations', True)
# else:
#     setattr(args, 'two_populations', False)
#     pass


if 'method_number' in args:
    setattr(args, 'method_args', args.method_number[1:])
    args.method_number = args.method_number[0]
    pass


# set stdout to the -o file, if applicable
sys.stdout = args.output_file


if 'varfile' in args or args.ms:

    if args.varfile != None and args.varfile_is_esp:
        original_variant_file = args.varfile
        pop = 'EA' if 'EA' in args.varfile.name else 'AA'
        var_header = original_variant_file.readline().strip().split('\t')
        args.varfile.seek(0)
        individuals = range(len(var_header) - 4)
        individuals_labels = ['%s_ind_%d' % (pop, ind) for ind in individuals]
        individuals_pops = [pop for ind in individuals]

    elif args.varfile != None and args.varfile_is_1kg:
        original_variant_file = args.varfile
        # get sample information
        path = '' if os.path.dirname(args.varfile.name) == '' else os.path.dirname(args.varfile.name) + os.path.sep
        sample_file = open(path + 'SampleInfo.txt', 'r')
        sys.stderr.write('sampleinfo: %s\n' % sample_file.name)
        sample_file.readline()
        pops = []
        individual_id = []
        for line in sample_file:
            line = line.strip().split()
            # sys.stderr.write('split line: %s\n' % line)
            pops.append(line[6])
            individual_id.append(line[1])
            pass
        var_header = original_variant_file.readline().strip().split('\t')
        args.varfile.seek(0)
        individuals = range(len(var_header) - 4)
        individuals_labels = ['%s_ind_%d' % (pops[ind], ind) for ind in individuals]
        if args.use_1kg_id: individuals_labels = [individual_id[ind] for ind in individuals]
        individuals_pops = [pops[ind] for ind in individuals]

    elif args.varfile != None and args.varfile_is_tair10:
        original_variant_file = args.varfile
        var_header = original_variant_file.readline().strip().split('\t')
        args.varfile.seek(0)
        individuals = range(len(var_header) - 6)
        individuals_labels = ['ind_%d' % ind for ind in individuals]
        individuals_pops = ['tair' for ind in individuals]

    elif args.varfile != None:
        original_variant_file = args.varfile
        var_header = original_variant_file.readline().strip().split('\t')
        individuals = var_header[8:]
        # keep chimp here?
        individuals_labels = [args.baselookup.pop_mapping[ind] for ind in individuals]
        individuals_pops = [args.baselookup.pop_mapping_general[ind] for ind in individuals]

    elif args.ms:
        original_ms_file = args.ms_file
        ms_header = original_ms_file.readline().strip()
        if args.forward_ms_file: print ms_header
        ms_header = ms_header.split()
        if args.debug: print "ms header:", ms_header
        num_chrs = int(ms_header[1])
        setattr(args, 'ms_num_sims', int(ms_header[2]))
        setattr(args, 'ms_sim_labels', ['%ssim%%0%dd' % (args.sim_prefix, len(str(args.ms_num_sims))) % i for i in range(args.ms_num_sims)])
        if args.debug: print "simulation labels:", args.ms_sim_labels

        if '-I' in ms_header:
            i_index = ms_header.index('-I')
            num_pops = int(ms_header[i_index+1])
            pop_sizes_from_header = [int(i) for i in ms_header[i_index+2:i_index+2+num_pops]]
            pass
#         if num_chrs % 2 != 0:
#             print "even number of chrs required!", num_chrs, num_chrs % 2
#             sys.exit(-1)
#             pass

        if args.ms_pops == None:
            args.ms_pops = [num_chrs]
            pass
        elif sum(args.ms_pops) != num_chrs:
            print "sum(ms_pops) must equal num_chrs:", sum(args.ms_pops), '!=', num_chrs
            sys.exit(-1)
            pass
        elif '-I' in ms_header and args.ms_pops != pop_sizes_from_header and args.ms_pops_force_match:
            print "ms_pops must equal ms_header pops:", args.ms_pops, '!=', pop_sizes_from_header
            sys.exit(-1)
            pass

        if args.debug: print "splitting ms chrs into individuals"
        individuals = []
        individuals_labels = []
        individuals_pops = []
        setattr(args, 'ms_chr_pairs', [])
        sum_popsize = 0
        for pop_num, popsize in enumerate(args.ms_pops):
            if args.debug: print "current popsize:", popsize
            if args.debug: print "chr range:", range(num_chrs)[sum_popsize:sum_popsize+popsize]
            if args.ms_no_randomize and args.ms_random_seed != None:
                print "Cannot set both ms-no-randomize and ms-random-seed!"
                sys.exit(-1)
            elif args.ms_no_randomize:
                rsample = range(num_chrs)[sum_popsize:sum_popsize+popsize]
            elif args.ms_random_seed != None:
                random.seed(args.ms_random_seed)
                rsample = random.sample(range(num_chrs)[sum_popsize:sum_popsize+popsize], popsize)
            else:
                rsample = random.sample(range(num_chrs)[sum_popsize:sum_popsize+popsize], popsize)
                pass
            sum_popsize = sum_popsize + popsize
            
            if args.debug: print "randomized chrs:", rsample
            # don't allow non-even pop sizes (because popsize is actually #chrs).  should we only allow this in special circumstances?  like if popsize == 1, and it's the last chr?
            if False and popsize % 2 != 0 and not (popsize == 1 and pop_num+1 == len(args.ms_pops)):
                print "Uneven population sizes only allowed if popsize == 1 and it's the last population."
                print pop_num, len(args.ms_pops), popsize
                sys.exit(-1)
                pass
            if popsize % 2 != 0: continue
            for i in range(popsize//2):
                args.ms_chr_pairs.append((rsample[i*2], rsample[i*2+1]))
                individuals.append('pop%d_%d' % (pop_num+1, i+1))
                individuals_labels.append('pop%d_%d' % (pop_num+1, i+1))
                individuals_pops.append('%d' % (pop_num+1))
                pass
            
            if args.debug: print "paired chrs:", args.ms_chr_pairs
            pass

        if args.ms_outgroup != None:
            if args.debug: print "setting ms archaic info", args.ms_outgroup
            args.ms_outgroup = (int(args.ms_outgroup[0]), float(args.ms_outgroup[1]))
            setattr(args, 'ms_arc', True)
        else:
            setattr(args, 'ms_arc', False)
            pass

        ## set up ms data
        def set_up_next_ms_sim():
            if args.debug: print "getting next ms_simulation"

            ## save current chr
            if 'ms_current_chr_num' in args:
                args.ms_current_chr_num += 1
            else:
                setattr(args, 'ms_current_chr_num', 0)
                pass
            setattr(args, 'ms_chr', '%ssim%%0%dd' % (args.sim_prefix, len(str(args.ms_num_sims))) % args.ms_current_chr_num)

            ## set up region length
            if args.base_sample_probability_normal != None:
                reg_len = int(random.normalvariate(args.base_sample_probability_normal[0], args.base_sample_probability_normal[1]))
                if args.debug: print 'setting random region length from normal distribution', reg_len, args.base_sample_probability_normal
                if reg_len > args.ms_reglen: reg_len = args.ms_reglen
                if reg_len < 0: reg_len = 0
                #if args.ms_current_chr_num == 1: reg_len = 0
                args.base_sample_probability = reg_len / args.ms_reglen
                setattr(args, 'ms_region_length', reg_len)
                pass
            elif args.ms_region_length_distribution != None:
                reg_len = random.sample(args.ms_region_length_distribution, 1)[0]
                if args.debug: print 'setting random region length from distribution file', reg_len
                if reg_len > args.ms_reglen: reg_len = args.ms_reglen
                if reg_len < 0: reg_len = 0
                #if args.ms_current_chr_num == 1: reg_len = 0
                args.base_sample_probability = reg_len / args.ms_reglen
                setattr(args, 'ms_region_length', reg_len)
                pass

            ## find start of simulation (//)
            ms_line = original_ms_file.readline().strip()
            if args.forward_ms_file: print >> ms_file_backlog[args.ms_chr], ms_line
            if args.debug: print "starting ms file backlog for:", args.ms_chr, ms_file_backlog.keys()
            while len(ms_line) == 0 or ms_line[0] != '//':
                ms_line = original_ms_file.readline().strip()
                if args.forward_ms_file: print >> ms_file_backlog[args.ms_chr], ms_line
                ms_line = ms_line.split()
                pass

            ## is there a tree?
            ## currently only interested in the tree if we're looking at archaics
            setattr(args, 'ms_trees', False)
            if '-T' in ms_header and args.ms_arc:
                if 'ms_trees' not in args: setattr(args, 'ms_trees', {})
                if 'ms_outgroup_branch_len' not in args: setattr(args, 'ms_outgroup_branch_len', {})
                if 'ms_outgroup_branch_len_by_win' not in args: setattr(args, 'ms_outgroup_branch_len_by_win', {})
                if 'ms_outgroup_branch_len_by_win_by_chr' not in args: setattr(args, 'ms_outgroup_branch_len_by_win_by_chr', {})
                #if 'ms_archaic_regions' not in args: setattr(args, 'ms_archaic_regions', defaultdict(lambda : mydefaultdict(lambda : False)))
                #if 'ms_archaic_regions_lens' not in args: setattr(args, 'ms_archaic_regions_lens', defaultdict(lambda : 0))
                ## reset this every time we're reading a new ms experiment
                default_bitarray = bitarray(args.ms_reglen)
                default_bitarray.setall(0)
                #setattr(args, 'ms_archaic_regions', defaultdict(lambda : mydefaultdict(lambda : False)))
                setattr(args, 'ms_trees', True)
                setattr(args, 'ms_archaic_regions', defaultdict(lambda : default_bitarray.copy()))
                #setattr(args, 'ms_archaic_regions', {})
                setattr(args, 'ms_archaic_regions_lens', defaultdict(lambda : 0))

                if '-r' in ms_header:
                    setattr(args, 'ms_multiple_trees', True)
                    tree_bases_count = 0
                    tree_bases_introgressed_count = 0
                    tree_bases_introgressed_count_by_win = {}
                    tree_bases_introgressed_count_by_win_by_chr = {}
                    while tree_bases_count < args.ms_reglen:
                        full_ms_tree = original_ms_file.readline().strip()
                        if args.forward_ms_file: print >> ms_file_backlog[args.ms_chr], full_ms_tree

                        #bases2 = re.findall(r'\[\d+\]', full_ms_tree)[0]
                        bases = full_ms_tree[:full_ms_tree.index(']')+1]
                        #print 'findall', bases, bases2
                        ms_tree = full_ms_tree[len(bases):]
                        bases = int(bases[1:-1])
                        prev_bases_count = tree_bases_count
                        tree_bases_count += bases


                        if args.ms_outgroup != None and 'ms_win_params' in args:
                            # if this region has been introgressed in any individuals, select the chrs that have been introgressed
                            # ms_archaic_chrs = pct_arc.extract_introgressed_chrs_from_tree(full_ms_tree, args.ms_outgroup[0], args.ms_outgroup[1])
                            ms_archaic_chrs = pct_arc.extract_introgressed_chrs_from_tree(ms_tree, args.ms_outgroup[0], args.ms_outgroup[1])
                            # now mark all relevant positions as being introgressed for those chrs
                            for chr in ms_archaic_chrs:
                            #if (chr-1) not in args.ms_archaic_regions: args.ms_archaic_regions[chr-1] = default_bitarray.copy()
                                for i in range(prev_bases_count, tree_bases_count):
                                ## the chrs returned by pct_arc are 1-based, but everywhere in my code they're 0-based
                                    args.ms_archaic_regions[chr-1][i] = True
                                    pass
                                pass
                            
                            # print 'read a recomb region,', prev_bases_count, tree_bases_count, 'archaic chrs:', ms_archaic_chrs
                            # print 'read a recomb region, archaic regions:', args.ms_archaic_regions

                            split_tree = ms_tree.split('%d:' % args.ms_outgroup[0])
                            #print ms_tree
                            #print split_tree
                            end_dist = min((x if x >=0 else 10000000000 for x in (split_tree[1].find(','), split_tree[1].find('('), split_tree[1].find(')'))))
                            bls = split_tree[1][:end_dist]
                            bl = float(bls)
                            #print 'findall', float(re.findall(r'[\d\.]+', split_tree[1])[0]), bl
                            if args.debug_2: print "\ninferred outgroup branch length:", bl
                            if args.debug_2: print "span for these bases %d-%d:" % (prev_bases_count, tree_bases_count)
                            if bl < args.ms_outgroup[1]: tree_bases_introgressed_count += bases

                            # get the number of bases with low branch len for each window
                            wincountstart = max(0, (prev_bases_count - args.ms_win_params[0] + args.ms_win_params[1]) // args.ms_win_params[1])
                            wincountend = (tree_bases_count) // args.ms_win_params[1]
                            if args.debug_2: print "windows %d-%d" % (wincountstart, wincountend)
                            for wincount in range(wincountstart, wincountend+1):
                                winstart = wincount * args.ms_win_params[1]
                                winend = wincount * args.ms_win_params[1] + args.ms_win_params[0]
                                hash = '%d-%d' % (winstart, winend)
                                if hash not in tree_bases_introgressed_count_by_win:
                                    tree_bases_introgressed_count_by_win[hash] = 0
                                    tree_bases_introgressed_count_by_win_by_chr[hash] = defaultdict(lambda : 0)
                                    pass
                                if bl < args.ms_outgroup[1]:
                                    tree_bases_introgressed_count_by_win[hash] += min(winend, tree_bases_count) - max(winstart, prev_bases_count)
                                    for chr in ms_archaic_chrs:
                                        # have to subtract one, to get chrs on zero based numbering
                                        tree_bases_introgressed_count_by_win_by_chr[hash][chr-1] += min(winend, tree_bases_count) - max(winstart, prev_bases_count)
                                        pass
                                    pass
                                if args.debug_2: print "window count, hash, basecount:", wincount, hash, tree_bases_introgressed_count_by_win[hash]
                                pass
                            pass
                        
                        #print bases, bl, tree_bases_introgressed_count
                        pass
                    ## save info on the amount of sequence introgressed for this simulation
                    ## (have to do it in a hash by chr num, because data is reset before any calc_* function is called)
                    args.ms_outgroup_branch_len[args.ms_chr] = tree_bases_introgressed_count
                    args.ms_outgroup_branch_len_by_win[args.ms_chr] = tree_bases_introgressed_count_by_win
                    args.ms_outgroup_branch_len_by_win_by_chr[args.ms_chr] = tree_bases_introgressed_count_by_win_by_chr
                    pass
                else:
                    ms_tree = original_ms_file.readline().strip()
                    if args.forward_ms_file: print >> ms_file_backlog[args.ms_chr], ms_tree

                    if args.debug: print "ms_tree %s" % args.ms_chr, ms_tree
                    args.ms_trees[args.ms_chr] = ms_tree
                    if args.ms_outgroup != None:
                        split_tree = ms_tree.split('%d:' % args.ms_outgroup[0])
                        if args.debug: print "tree split on outgroup (%d)" % args.ms_outgroup[0], split_tree[1]
                        args.ms_outgroup_branch_len[args.ms_chr] = re.findall(r'[\d\.]+', split_tree[1])[0]
                        if args.debug: print "inferred outgroup branch length:", args.ms_outgroup_branch_len
                        pass
                    pass
                pass

            ## get the number of variant sites
            while ms_line[0] != 'segsites:':
                ms_line = original_ms_file.readline().strip()
                if args.forward_ms_file: print >> ms_file_backlog[args.ms_chr], ms_line
                ms_line = ms_line.split()
                pass
            setattr(args, 'ms_num_sites', int(ms_line[1]))
            setattr(args, 'ms_current_site', 0)

            ## get the positions of the variant sites
            ms_line = original_ms_file.readline().strip()
            if args.forward_ms_file: print >> ms_file_backlog[args.ms_chr], ms_line
            ms_line = ms_line.split()

            setattr(args, 'ms_snp_poses', [float(p) for p in ms_line[1:]])
            ## 0 and 1 are possible positions in the ms file, so e.g. with a reglen of 200 you could get poses 0 and 200, which is actually 201 possible positions.
            # to protect against this, multiply p by reglen-1 (e.g., 0 to 199)
            args.ms_snp_poses = [int((args.ms_reglen-1) * p) for p in args.ms_snp_poses]
            for i,p in enumerate(args.ms_snp_poses):
                if i+1 < len(args.ms_snp_poses) and args.ms_snp_poses[i+1] == p:
                    args.ms_snp_poses[i+1] = p+1
                    if p+1 >= args.ms_reglen:
                        sys.stderr.write('bad snp position (%d, out of range of the sequence length, 0-%d) due to very dense ms sampling\n' % (p+1, args.ms_reglen-1))
                        #sys.exit(-1)
                        pass
                    pass
                pass
            #if args.debug: print "ms snp positions:", args.ms_snp_poses

            ## save ms snp data
            setattr(args, 'ms_chrs', [])
            for i in range(num_chrs):
                ms_line = original_ms_file.readline().strip()
                if args.forward_ms_file: print >> ms_file_backlog[args.ms_chr], ms_line
                args.ms_chrs.append(ms_line)
                pass

            # print ms_line
            # ## go through the rest of the file until a blank line, which should signify the next simulation
            # if args.forward_ms_file: ms_line = original_ms_file.readline().strip()
            # print ms_line
            # while args.forward_ms_file and len(ms_line) != 0:
            #     if ms_line[0] == '//':
            #          print "error in ms file logic! looking for forwarded ms data, and then a blank line, but found the startf ms simulations(//)"
            #          sys.exit(-1)
            #          pass
            #     if args.forward_ms_file: print >> ms_file_backlog[args.ms_chr], ms_line
            #     ms_line = original_ms_file.readline().strip()
            #     pass

            ## end of set_up_next_ms_sim
            return

        set_up_next_ms_sim()
        
        pass

    original_ind_labels = list(individuals_labels)
    original_ind_pops = list(individuals_pops)
    
    # if args.filter_inds != None:
    #     # remove inds from header, set up list of sites to remove (or keep) in read_snp
    #     pass

    num_genotypes = len(individuals)

    setattr(args, 'filter_inds', None)

    if args.keep_pops != None:
        # select individuals to keep
        # This is a little confusing, but FILTER INDS should be the only thing that's indexed to the original list of individuals.
        # All other lists need to be indexed to the updated list of individuals.
        args.filter_inds = [ind for ind in range(len(individuals)) if individuals_pops[ind] in args.keep_pops]
        # recalculate shit
        ## THIS ONLY WORKS THE FIRST TIME?  I.E., NOT IN THE COMMENTED OUT STUFF BELOW?
        individuals = [individuals[ind] for ind in args.filter_inds]
        num_genotypes = len(individuals)
        #individuals_labels = [args.baselookup.pop_mapping[ind] for ind in individuals]
        #individuals_pops = [args.baselookup.pop_mapping_general[ind] for ind in individuals]
        individuals_labels = [original_ind_labels[ind] for ind in args.filter_inds]
        individuals_pops = [original_ind_pops[ind] for ind in args.filter_inds]
        pass

    if args.keep_inds != None or args.keep_inds_by_id != None:
        if args.keep_inds_by_id != None:
            ## keep inds is one-based (to make it easier to give numbers on the command line)
            setattr(args, 'keep_inds', [i+1 for i,ind in enumerate(individuals_labels) if ind in args.keep_inds_by_id])
        else:
            args.keep_inds = sorted(args.keep_inds)
            pass
        args.filter_inds   = [args.filter_inds[i-1] for i in args.keep_inds]
        individuals[:]        = [individuals[i-1] for i in args.keep_inds]
        individuals_labels[:] = [individuals_labels[i-1] for i in args.keep_inds]
        individuals_pops[:]   = [individuals_pops[i-1] for i in args.keep_inds]
        num_genotypes = len(individuals)
        pass
    
    if args.add_inds != None or args.add_inds_by_id != None:
        if args.add_inds_by_id != None:
            ## add inds is zero-based (but one based on the command line)
            setattr(args, 'add_inds', [i for i,ind in enumerate(original_ind_labels) if ind in args.add_inds_by_id])
        else:
            ## add inds is zero-based (but one based on the command line)
            args.add_inds = sorted([i-1 for i in args.add_inds])
            pass
        if args.filter_inds == None:
            args.filter_inds = sorted(args.add_inds)
        else:
            args.filter_inds = sorted(list(set(args.filter_inds).union(set(args.add_inds))))
            pass

        #print
        #print 'debug!!!!!!!!!!!!!!'
        #print args.add_inds
        #print args.filter_inds

        individuals = list(args.filter_inds)
        individuals_labels = [original_ind_labels[ind] for ind in args.filter_inds]
        individuals_pops = [original_ind_pops[ind] for ind in args.filter_inds]
        num_genotypes = len(individuals)
        
        #print individuals

        pass
    
    sys.stderr.write('Retained only %s\n' % (str(args.keep_pops)))
    sys.stderr.write('Number of inds: %d\n' % len(individuals))
    sys.stderr.write('etc: %s %s\n' % (str(individuals_labels), str(individuals_pops)))
        
    # if args.individual_limit != None:
    #     # select individuals to keep
    #     # This is a little confusing, but FILTER INDS should be the only thing that's indexed to the original list of individuals.
    #     # All other lists need to be indexed to the updated list of individuals.
    #     args.filter_inds = [ind for ind in range(len(individuals)) if ind < args.individual_limit and (args.filter_inds == None or original_ind_labels[ind] not in args.filter_inds)]
    #     # recalculate shit
    #     individuals = [individuals[ind] for ind in args.filter_inds]
    #     num_genotypes = len(individuals)
    #     #individuals_labels = [args.baselookup.pop_mapping[ind] for ind in individuals]
    #     #individuals_pops = [args.baselookup.pop_mapping_general[ind] for ind in individuals]
    #     individuals_labels = [original_ind_labels[ind] for ind in args.filter_inds]
    #     individuals_pops = [original_ind_pops[ind] for ind in args.filter_inds]
    #     sys.stderr.write('NOT SURE IF THIS WORKS WITH -POPS')
    #     sys.stderr.write('Retained only %s\n' % (str(args.keep_pops)))
    #     sys.stderr.write('Number of inds: %d\n' % len(individuals))
    #     sys.stderr.write('etc: %s %s\n' % (str(individuals_labels), str(individuals_pops)))
    #     pass

#     for i in range(num_genotypes):
#         print i, individuals[i]
#         pass

#     sys.exit()

# THIS IS CURRENTLY BROKEN, BECAUSE IT UPDATES FILTER INDS BASED ON THE MODIFIED LIST OF INDIVIDUALS, NOT THE ORIGINAL LIST OF INDIVIDUALS
#     if args.remove_pops != None:
#         # select individuals to keep
#         args.filter_inds = [ind for ind in range(len(individuals)) if individuals_pops[ind] not in args.remove_pops]
#         # recalculate shit
#         individuals = [individuals[ind] for ind in args.filter_inds]
#         num_genotypes = len(individuals)
#         individuals_labels = [args.baselookup.pop_mapping[ind] for ind in individuals]
#         individuals_pops = [args.baselookup.pop_mapping_general[ind] for ind in individuals]
#         sys.stderr.write('Removed pos %s\n' % (str(args.remove_pops)))
#         sys.stderr.write('Number of inds: %d\n' % len(individuals))
#         sys.stderr.write('etc: %s %s\n' % (str(individuals_labels), str(individuals_pops)))
#         pass

    if 'target_pops' in args and args.target_pops != None:
        # select individuals for the target population
        pops = [ind for ind in range(len(individuals)) if individuals_pops[ind] in args.target_pops]
        non_target_pop = [ind for ind in range(len(individuals)) if individuals_pops[ind] not in args.target_pops]
        sys.stderr.write('Target Pop: Number of inds: %d\n' % len(pops))
        sys.stderr.write('Non-Target Pop: Number of inds: %d\n' % len(non_target_pop))
        #sys.stderr.write('Target Pop: %s\n' % str(pops))
        #sys.stderr.write('Non-Target Pop: %s\n' % str(non_target_pop))
        #sys.stderr.write('etc: %s %s\n' % (str([args.baselookup.pop_mapping[individuals[ind]] for ind in pops]), str([args.baselookup.pop_mapping_general[individuals[ind]] for ind in pops])))
        args.target_pops = pops
        setattr(args, 'non_target_pop', non_target_pop)
        pass

    # if args.two_populations:
    #     sys.stderr.write('Splitting individuals into two population groups: %s and %s\n' % (str(args.population1), str(args.population2)))
    #     args.population1 = p1 = [ind for ind in range(len(individuals)) if individuals_pops[ind] in args.population1]
    #     args.population2 = p2 = [ind for ind in range(len(individuals)) if individuals_pops[ind] in args.population2]
    #     # print individuals_pops
    #     # print p1
    #     # print p2
    #     sys.stderr.write('Number of people in each population: %d and %d\n' % (len(p1), len(p2)))
    #     pass
    
    if args.multiple_populations:
        sys.stderr.write('Splitting individuals into several population groups: %s\n' % (str(args.multiple_populations)))
        setattr(args, 'mpops', args.multiple_populations)
        args.multiple_populations = [p.split(',') for p in args.multiple_populations]
        for i, p in enumerate(args.multiple_populations):
            #print i,p
            args.multiple_populations[i] = [ind for ind in range(len(individuals)) if individuals_pops[ind] in p]
            pass
        sys.stderr.write('Number of people in each population: %s\n' % (', '.join((str(len(p)) for p in args.multiple_populations))))
        pass
        
    if 'is_male' in args and args.varfile_is_tair10:
        is_male = [False for ind in individuals]
        num_males = sum(is_male)

    elif 'is_male' in args and args.varfile_is_1kg:
        is_male = []
        path = '' if os.path.dirname(args.varfile.name) == '' else os.path.dirname(args.varfile.name) + os.path.sep
        sample_file = open(path + 'SampleInfo.txt', 'r')
        sample_file.readline()
        pops = []
        for line in sample_file:
            line = line.strip().split()
            is_male.append(line[4] == '1')
            pass
        num_males = sum(is_male)

    elif 'is_male' in args and args.varfile_is_esp:
        pop = 'EA' if 'EA' in args.varfile.name else 'AA'
        is_male = []
        path = '' if os.path.dirname(args.varfile.name) == '' else os.path.dirname(args.varfile.name) + os.path.sep
        is_male_file = open(path + 'SampleInfo.txt', 'r')
        is_male_file.readline()
        for line in is_male_file:
            line = line.strip().split()
            if line[1] != pop: continue
            is_male.append(line[2] == '0')
            pass
        num_males = sum(is_male)
        pass

    elif 'is_male' in args and not args.ms:
        is_male = []
        for ind in individuals:
            is_male_file = open(args.is_male.replace('IND', ind), 'r')
            is_male_num = int(is_male_file.readline().strip())
            is_male.append( is_male_num > 0 )
            pass
        num_males = sum(is_male)
        pass


    pass

if args.regions != None: sys.stderr.write('Initial number of bases considered (regions): %s\n' % locale.format('%d', args.regions.total_length(), grouping=True))

# pad before everything else!
if args.regions != None and args.pad_regions != None:
    #print 'pre padding regions', args.regions.total_length()
    args.regions.pad_beds_by(args.pad_regions)
    #print 'post padding regions', args.regions.total_length()
    pass

if args.regions != None and args.intersect_region != None:
#     for ir in args.intersect_region.bases:
#         args.regions.bases &= ir
#         pass
    args.regions.bases &= args.intersect_region.bases
    args.intersect_region = None
    # print 'post intersect regions', args.regions.total_length()
    pass
# if args.regions != None and args.second_region != None:
#     args.second_region.bases &= args.regions.bases
#     pass
if args.regions != None and args.exclude_region != None:
    sys.stderr.write('Excluding %s bases.\n' % locale.format('%d', args.exclude_region.total_length(), grouping=True))
    args.regions.bases &= ~args.exclude_region.bases
    args.exclude_region = None
elif args.regions == None and args.exclude_region != None:
    sys.stderr.write('Excluding %s bases.\n' % locale.format('%d', args.exclude_region.total_length(), grouping=True))
    args.regions = args.exclude_region
    args.exclude_region = None
    args.regions.bases.invert()
    pass

# only consider chimp bases that overlap with regions
if args.regions != None and 'chimp_mapped' in args:
    args.chimp_mapped.bases &= args.regions.bases
    pass

# exclude nocall data after masking chimp (because this data is specific to human data - for most applications we don't want to apply it to the chimp data)
if args.exclude_all_nocalls != None and len(individuals) > 0:
    sys.stderr.write('Excluding 100% of nocalls (and setting missingness threshold to 0%)\n')
    #sys.stderr.write('Excluding 100% of nocalls (AND NOOOOOOOOOOOOOOOOOOOOOOT SETTING MISSINGNESS THRESHOLD TO 0%)\n')
    bt = myBedTools(myBedTools.binarybedfilegenome, initialize = False)
    for ind in individuals:
        f = args.exclude_all_nocalls.replace('IND', ind)
        bt.read_as_regions(f)
        pass
    if args.regions != None:
        args.regions.bases &= ~bt.bases
    else:
        args.regions = bt
        args.regions.bases.invert()
        pass
    bt = None
    args.filter_missingness = 1
    pass

# make sure to do this after all modifications to regions are complete - intersection, etc.
if args.regions != None and args.second_region != None:
    args.second_region.bases &= args.regions.bases
    pass
if args.regions != None and args.func == 'mk':
    # only take the neutral or selected regions that are also in the regions file
    # print 'pre neutral/regions intersect', args.neutral_regions.total_length()
    args.neutral_regions.bases &= args.regions.bases
    args.selected_regions.bases &= args.regions.bases
    # print 'post neutral/regions intersect', args.neutral_regions.total_length()
    pass
if args.func == 'mk':
    # adjust neutral and selected regions

    # remove selected regions from neutral regions
    args.neutral_regions.bases &= ~args.selected_regions.bases
    # print 'post remove selected from neutral', args.neutral_regions.total_length()

    # remove chimp unmapped regions
    args.neutral_regions.bases &= args.chimp_mapped.bases
    args.selected_regions.bases &= args.chimp_mapped.bases
    # print 'post remove chimp unmapped from neutral', args.neutral_regions.total_length()

    # remove exclude regions
    if args.exclude_region != None:
        args.neutral_regions.bases &= ~args.exclude_region.bases
        args.selected_regions.bases &= ~args.exclude_region.bases
        # print 'post remove exclude from neutral', args.neutral_regions.total_length()
        pass

    # remove neutral-specific exclude regions
    if args.neutral_regions_exclude != None:
        args.neutral_regions.bases &= ~args.neutral_regions_exclude.bases
        # print 'post remove neutral-specific exclude from neutral', args.neutral_regions.total_length()
        pass
        
    #print 'after anti-intersect', args.neutral_regions.total_length(), args.selected_regions.total_length()
    #not_sel = ~args.selected_regions.bases
    #intersected = args.neutral_regions.bases & not_sel
    setattr(args, 'snp_regions', myBedTools(myBedTools.binarybedfilegenome, invert = False, initialize = True))
    pass

if args.regions != None: sys.stderr.write('Final number of bases considered (regions, exclude, intersect, etc): %s\n' % locale.format('%d', args.regions.total_length(), grouping=True))


consider_ancestry = False
consider_ancestry_bsg = False
if 'ancestry_file' in args:
    ancestral_state_file = args.ancestry_file
    ancestral_column = args.ancestral_column
    consider_ancestry = True
    pass
if 'ancestry_file_bsg' in args:
    consider_ancestry = True
    consider_ancestry_bsg = True
    pass

if 'chimp_snps' in args:
    chimp_snps_file = args.chimp_snps
    pass

if 'chimp_mapped' in args:
    chimp_mapped_regions = args.chimp_mapped
    pass

if 'time_to_primate_divergence' in args:
    time_to_primate_divergence = args.time_to_primate_divergence
    pass

if 'nocall_file' in args:
    nocall_file_name = args.nocall_file
    pass

if 'variant_type' not in args:
    # This is not currently a setting - all of our analyses use snps and only snps.
    # But if we wanted to include alternative variant types, we'd need to add this setting.
    only_snps = True
    pass


def add_to_db(keys, vals):

    keys += ['date']
    vals += [time.time()]
    conn = sqlite3.connect(args.data_base_name, timeout=360)
    cur = conn.cursor()
    
    fieldnamelist = ','.join(keys)
    cur.execute( 'create table if not exists %s (%s);' % (args.table_name, fieldnamelist) )

    placeholder = '?'
    placeholderlist = ','.join([placeholder] * len(keys))
    
    insert = 'insert into %s (%s) values (%s);' % (args.table_name, fieldnamelist, placeholderlist)
    cur.execute(insert, vals)

    conn.commit()
    conn.close()
    pass


def sample_wr(population, k):
    "Chooses k random elements (with replacement) from a population"
    n = len(population)
    _random, _int = random.random, int  # speed hack
    return [population[_int(_random() * n)] for i in itertools.repeat(None, k)]

def calc_chr_missing_etc(chr, pos, subset = None):

    if args.ms:
        if subset == None: return [num_genotypes*2, 0, 'A']
        return [len(subset) * 2, 0, 'A']
    
    if subset == None:
        num_g = num_genotypes
        num_m = num_males
    else:
        num_g = len(subset)
        num_m = sum([is_male[i] for i in subset])
        pass
    
    if chr == 'chrX' and pos > args.baselookup.chrX_par1 and pos <= args.baselookup.chrX_par2:
        num_chrs = num_g * 2 - num_m
        sex_chr = 'X'
    elif chr == 'chrY':
        num_chrs = num_m
        sex_chr = 'Y'
    else:
        num_chrs = num_g * 2
        sex_chr = 'A'
        pass
    max_missing_chrs = (1-args.filter_missingness) * num_chrs

    return [num_chrs, max_missing_chrs, sex_chr]

# can be used as a generator to read all snps in a file, as a list
def read_snps_gen(read_fn):
    s = read_fn()
    while s:
        yield s
        s = read_fn()
        pass
    return


def read_snp():

    s = read_snp_helper()

    while s == 'bad snp':
        s = read_snp_helper()
        pass

    if s == 'eof':
        return False

    return s


## calculate vars, refs, etc
def get_refs_vars_haps_hets(s, subset = None):
    # print 'in get_refs_vars_haps_hets', s
    if subset != None:
        #print s, subset
        s = [s[i] for i in subset]
        pass
    refs = sum([x[1] for x in s])
    vars = sum([x[2] for x in s])
    haps = refs + vars
    hets = sum([1 for x in s if x[1] == 1 and x[2] == 1])
    chrs = sum([sum(t) for t in s])
    return (refs, vars, haps, hets, chrs)


def calc_x(s, subset = None):
    if subset == None:
        v = s['vars']
        r = s['refs']
    else:
        (r, v, _, _, _) = get_refs_vars_haps_hets(s['snp_details'], subset)
        pass

    # x, the putative number of derived aleles at this site
    # unless r is 0 or v is 0, in which case this site is fixed for the individuals in subset
    if r == 0 or v == 0:
        x = 0
    elif s['mrca'] == 0:
        x = v
    elif s['mrca'] == 1:
        x = r
    else:
        print 'mrca is not zero or one!', s, v, r
        pass

    return x

setattr(args, 'rejected', 0)
def read_snp_helper_ms():

    if args.debug_3: print
    
    if args.debug_3 or args.debug_read_snps: print 'reading a snp', args.ms_current_site, args.ms_num_sites, args.ms_current_chr_num, args.ms_num_sims-1

    if args.ms_current_site >= args.ms_num_sites and args.ms_current_chr_num == args.ms_num_sims-1:
        if args.debug_3 or args.debug_read_snps: print 'ms sim eof'
        return 'eof'
    elif args.ms_current_site >= args.ms_num_sites:
        if args.debug_3 or args.debug_read_snps: print 'setting up next ms simulation', args.ms_current_site, args.ms_num_sites, args.ms_current_chr_num, args.ms_num_sims-1
        set_up_next_ms_sim()
        return 'bad snp'
    elif args.ms_snp_poses[args.ms_current_site] >= args.ms_reglen and args.ms_current_chr_num == args.ms_num_sims-1:
        sys.stderr.write('cutting off a snp because of too-dense site sampling at the end of a simulation (ALSO AT END OF FILE, SO EXITING): %d\n' % args.ms_snp_poses[args.ms_current_site])
        if args.debug_3 or args.debug_read_snps: print 'setting up next ms simulation'
        return 'eof'
    elif args.ms_snp_poses[args.ms_current_site] >= args.ms_reglen:
        sys.stderr.write('cutting off a snp because of too-dense site sampling at the end of a simulation: %d\n' % args.ms_snp_poses[args.ms_current_site])
        if args.debug_3 or args.debug_read_snps: print 'setting up next ms simulation'
        set_up_next_ms_sim()
        return 'bad snp'
    elif args.base_sample_probability != 1:
        r = random.random()
        if args.debug_3 or args.debug_read_snps: print '%f, %f' % (r, args.base_sample_probability)
        #args.snp_count += 1
        if r >= args.base_sample_probability:
            args.ms_current_site += 1
            if args.debug_2: print 'SKIP rejected due to base sample probability'
            args.rejected += 1
            #print args.rejected / args.snp_count
            #print 'rejecting snp'
            return 'bad snp'
        pass
    elif args.limit != None and args.ms_current_site > args.limit:
        return 'eof'
    # this essentially tests is_archaic
    if args.filt_retain_arc and args.ms_arc and not args.ms_chrs[args.ms_outgroup[0]-1][args.ms_current_site] == '1':
        if args.debug_3 or args.debug_read_snps: print 'SKIPPING snp because it is not found in the archaic chr:', args.ms_chr, args.ms_snp_poses[args.ms_current_site]
        args.ms_current_site += 1
        return 'bad snp'
    
    ret_data = {}

    just_snp = [0] * len(args.ms_chr_pairs)
    ## also save which ms chr each variant is on (or both, if it's homozygous derived in a given individual)
    ms_snp_haplotypes = [0] * len(args.ms_chr_pairs)
    archaic_inds = mydefaultdict(lambda : False)
    # one based position
    pos = args.ms_snp_poses[args.ms_current_site]
    # infer snp
    # print "getting snp state for %d pairs" % (len(args.ms_chr_pairs))
    for i in range(len(args.ms_chr_pairs)):
        c1 = args.ms_chr_pairs[i][0]
        c2 = args.ms_chr_pairs[i][1]
        num_anc = int(args.ms_chrs[c1][args.ms_current_site] == '0') + int(args.ms_chrs[c2][args.ms_current_site] == '0')
        num_der = int(args.ms_chrs[c1][args.ms_current_site] == '1') + int(args.ms_chrs[c2][args.ms_current_site] == '1')

        # if num_der > 0:
        #     print 'ms_chrs', i, c1, c2, len(args.ms_chrs)
        #     if c1 > 66: print 'ASIAN'
        #     print '   vars', pos, c1, int(args.ms_chrs[c1][args.ms_current_site] == '1')
        #     print '   vars', pos, c2, int(args.ms_chrs[c2][args.ms_current_site] == '1')
        #     pass

        just_snp[i] = (0, num_anc, num_der)
        ms_snp_haplotypes[i] = [c for c in args.ms_chr_pairs[i] if args.ms_chrs[c][args.ms_current_site] == '1']
        ## if this variant is introgressed in one of the two chrs, mark it as such for the individual
        ## then snp['archaic_inds'][ind] will be true if snp is introgressed in ind
        ## remember that c1 and c2 are zero based, whereas the ms chr numbers are 1-based (for debuging purposes)
        # print pos
        #if (chr-1) not in args.ms_archaic_regions: args.ms_archaic_regions[chr-1] = default_bitarray.copy()
        #if c1 in args.ms_archaic_regions and args.ms_archaic_regions[c1][pos] and args.ms_chrs[c1][args.ms_current_site] == '1': archaic_inds[i] = True
        #if c2 in args.ms_archaic_regions and args.ms_archaic_regions[c2][pos] and args.ms_chrs[c2][args.ms_current_site] == '1': archaic_inds[i] = True
        if args.ms_arc and args.ms_trees and args.ms_archaic_regions[c1][pos] and args.ms_chrs[c1][args.ms_current_site] == '1': archaic_inds[i] = True
        if args.ms_arc and args.ms_trees and args.ms_archaic_regions[c2][pos] and args.ms_chrs[c2][args.ms_current_site] == '1': archaic_inds[i] = True
        pass

    # print 'archaic chrs for snp', pos, archaic_inds

    # set up basic snp info
    chr = args.ms_chr
    ref = '0'
    var = '1'

    if args.debug_3 or args.debug_read_snps: print 'in read snp helper', chr, pos

    ret_data['chr'] = chr
    ret_data['pos'] = pos
    ret_data['ref'] = ref
    ret_data['var'] = var
    ret_data['snp_details'] = just_snp
    ret_data['snp_details_ms_haplotypes'] = ms_snp_haplotypes

    if args.ms_arc:
        ret_data['is_archaic'] = [args.ms_chrs[args.ms_outgroup[0]-1][args.ms_current_site] == '1']
        ret_data['archaic_snps'] = [args.ms_chrs[args.ms_outgroup[0]-1][args.ms_current_site]]
        ret_data['archaic_inds'] = archaic_inds
        ## what?
        args.ms_outgroup[0]
        pass

    if args.debug_3 or args.debug_read_snps: print 'snp details', ret_data

    ## if looking at a subset of individuals, remove unwanted inds now
    ## important!  if using pop1 pop2 AND filter_inds, then args.population1, etc, must be in terms of the filtered list of individuals.
    if args.debug_3: print "filter inds, number of inds", args.filter_inds, range(len(just_snp))
    if args.filter_inds:
        just_snp[:] = [just_snp[i] for i in range(len(just_snp)) if i in args.filter_inds]
        pass

    (refs, vars, haps, hets, chrs) = get_refs_vars_haps_hets(just_snp)

    ret_data['refs'] = refs
    ret_data['vars'] = vars
    ret_data['haps'] = haps
    ret_data['anc'] = '0'
    ret_data['mrca'] = 0
    ret_data['mrca_known'] = True
    ret_data['x'] = calc_x(ret_data)


    ## DON'T RETURN BEFORE THIS POINT UNLESS IT'S EOF OR YOU INCREMENT MS_CURRENT_SITE!!
    args.ms_current_site += 1
    

    


    ## should we skip this snp?
    if refs == 0 or vars == 0:
        # population is fixed for a snp
        if args.debug_2: print "SKIP population is fixed for a snp", refs, vars, just_snp, ret_data
        return 'bad snp'
    if hets == haps / 2 and num_genotypes > 1:
        # totally heterozygous - potentially paralogous
        # only skip if we're looking at more than one ind, though
        if args.debug_2: print "SKIP totally heterozygous - potentially paralogous", ret_data
        return 'bad snp'
    if args.min_freq > vars / haps:
        if args.debug_2: print "SKIP rejecting snp because of too low freq", ret_data
        return 'bad snp'
    if args.max_freq < vars / haps:
        if args.debug_2: print "SKIP rejecting snp because of too high freq", ret_data
        return 'bad snp'

    ## get two population values (add here if new values are needed - currently just for fst)
    # if args.two_populations:
    #     pop1_snp = [just_snp[i] for i in range(len(just_snp)) if i in args.population1]
    #     pop2_snp = [just_snp[i] for i in range(len(just_snp)) if i in args.population2]
    #     (p1_refs, p1_vars, p1_haps, p1_hets, p1_chrs) = get_refs_vars_haps_hets(pop1_snp)
    #     (p2_refs, p2_vars, p2_haps, p2_hets, p2_chrs) = get_refs_vars_haps_hets(pop2_snp)
    #     ret_data['p1_refs'] = p1_refs
    #     ret_data['p1_vars'] = p1_vars
    #     ret_data['p1_haps'] = p1_haps
    #     ret_data['p2_refs'] = p2_refs
    #     ret_data['p2_vars'] = p2_vars
    #     ret_data['p2_haps'] = p2_haps
    #     print "seems to be messing up in the Y chromosome (see fst results?)"
    #     sys.exit(-1)
    #     if p1_chrs - p1_haps > (1-args.filter_missingness) * p1_chrs: return 'bad snp'
    #     if p2_chrs - p2_haps > (1-args.filter_missingness) * p2_chrs: return 'bad snp'
    #     pass

    if args.multiple_populations != None:
        ret_data['mpops'] = []
        for pop in args.multiple_populations:
            pop_data = {}
            ret_data['mpops'].append(pop_data)
            pop1_snp = [just_snp[i] for i in range(len(just_snp)) if i in pop]
            (p1_refs, p1_vars, p1_haps, p1_hets, p1_chrs) = get_refs_vars_haps_hets(pop1_snp)
            pop_data['refs'] = p1_refs
            pop_data['vars'] = p1_vars
            pop_data['haps'] = p1_haps
            pass
        pass

    if args.report_progress != None and args.snp_count % args.report_progress == 0:
        sys.stderr.write('..at snp\t%d:\t%s\t%d;\t%.2f\n' % (args.snp_count, chr, pos, time.time() - args.check_time))
        args.check_time = time.time()
        pass
    args.snp_count += 1
        
    return ret_data

setattr(args, 'snp_count', 0)
def read_snp_helper_varfile():
    ret_data = {}

    snp = original_variant_file.readline().strip().split('\t')
    if args.debug: print snp[0:5]
    
    if len(snp) == 1:
        if args.debug or args.debug_read_snps: print "eof in read snps"
        return 'eof'
    elif args.base_sample_probability != 1:
        r = random.random()
        if r >= args.base_sample_probability:
            args.rejected += 1
            if args.debug or args.debug_read_snps: print "randomly rejecting snp", snp[0], snp[1]
            return 'bad snp'
        pass
    
    if args.varfile_is_esp:

        chr = 'chr%s' % snp[0]
        # one based position
        pos = int(snp[1])
        ref = snp[2]
        var = snp[3]
        snp_start = 4

        def get_genotype(genotype):
            if genotype == '0':
                return (0,2,0)
            if genotype == '1':
                return (0,1,1)
            if genotype == '2':
                return (0,0,2)
            if genotype == '9':
                return (2,0,0)
            print "error in esp genotype: should be 0,1,2,9:", genotype
            sys.exit(-1)
            return

    elif args.varfile_is_tair10:

        chr = snp[0]
        # one based position
        pos = int(snp[1])
        ref = snp[2]
        var = snp[3]
        snp_start = 6

        def get_genotype(genotype):
            if genotype == ref:
                return (0,2,0)
            elif genotype == var:
                return (0,0,2)
            else:
                return (2,0,0)
            print "error in tair genotype: should be 0,1,2,9:", genotype
            sys.exit(-1)
            return
        
    elif args.varfile_is_1kg:

        chr = 'chr' + snp[0]
        # one based position
        pos = int(snp[1])
        ref = snp[2]
        # for some reason, variants are encoded as A/T, etc, but the ref and var are not always in the same order
        r_v = snp[3].split('/')
        var = r_v[0] if r_v[1] == ref else r_v[1]
        
        num_alleles = 2
        num_calls = int(snp[5])
        snp_start = 4

        def get_genotype(genotype):
            if genotype == '0':
                # hom ref
                return (0,2,0)
            elif genotype == '3':
                # hom alt
                return (0,0,2)
            elif genotype == '1':
                # het 0|1
                return (0,1,1)
            elif genotype == '2':
                # het 1|0
                return (0,1,1)
            print "error in 1kg genotype: should be 0,1,2,3:", genotype
            sys.exit(-1)
            return

    else:
        
        if only_snps and snp[4] != 'snp':
            if args.debug_2 or args.debug_read_snps: print "non-snp in read snps", snp[4]
            return 'bad snp'
        
        chr = snp[1]
        # one based position
        pos = int(snp[3])
        ref = snp[5]
        var = snp[6]
        snp_start = 8
        
        def get_genotype(genotype):
            g = (genotype.count('N'), genotype.count('0'), genotype.count('1'))
            # filter out partial calls
            if args.fpc:
                if g[0] > 0: g = (sum(g),0,0)
                pass
            return g
        
        pass

    if args.debug_3 or args.debug_read_snps: print 'in read snp helper', chr, pos


    # rough check to see if we want this snp (speeds up calculation for a small region, like mhc)
    # also doesn't work for things like the global calculations.. I say skip it.
    # This function is provided by the regions function in the next statement
    # if args.baselookup.chr_nums[chr] < args.baselookup.chr_nums[current_chr]:
    #     return ret_data

    # see if this is in the regions file
    if args.regions != None and not args.regions.in_region_one_based(chr, pos):
        if args.debug_3 or args.debug_read_snps: print "SKIP not in regions", snp
        return 'bad snp'

    # see if this is in the exclude file
    if args.exclude_region != None and args.exclude_region.in_region_one_based(chr, pos):
        if args.debug_3 or args.debug_read_snps: print "SKIP in exclude", snp
        return 'bad snp'

    # see if this should be excluded by random chance
    # if random.random() >= args.base_sample_probability:
    #     return 'bad snp'
    

    ret_data['chr'] = chr
    ret_data['pos'] = pos
    ret_data['ref'] = ref
    ret_data['var'] = var
    #snp[8:] = [(genotype.count('N'), genotype.count('0'), genotype.count('1')) for genotype in snp[8:]]
    just_snp = [get_genotype(genotype) for genotype in snp[snp_start:]]
    ret_data['snp_details'] = just_snp

    if args.comp and args.compare_calls != None:
        ## this only checks to see if the VAR is archaic, not the REF
        ret_data['is_archaic'] = [var.upper() == bsg_file.get_base_one_based(chr, pos) for bsg_file in args.compare_calls]
        ret_data['archaic_snps'] = [bsg_file.get_base_one_based(chr, pos) for bsg_file in args.compare_calls]
    elif args.comp and args.compare_calls_vcf != None:
        ## this only checks to see if the VAR is archaic, not the REF
        if chr not in args.compare_calls_vcf:
            print "vcf file chrs do not match data chrs"
            print chr, args.compare_calls_vcf
            sys.exit(-1)
            pass
        ret_data['is_archaic'] = [var.upper() in args.compare_calls_vcf[chr][pos]]
        ret_data['archaic_snps'] = [args.compare_calls_vcf[chr][pos][0]]
        #print 'archaic', chr, pos, ref, var, ret_data['is_archaic'], ret_data['archaic_snps']

    else:
        ret_data['is_archaic'] = [False]
        ret_data['archaic_snps'] = ['N']
        pass

    # this essentially tests is_archaic
    # but we look at both REF and VAR
    # we do *another* test after mrca calculation, if necessary, to restrict down to matches with just ref or just var
    if args.filt_retain_arc and args.comp and ref not in ret_data['archaic_snps'] and var not in ret_data['archaic_snps']:
        if args.debug_3 or args.debug_read_snps: print 'SKIP snp because it is not found in the archaic chr:', chr, pos, ref, var, ret_data['archaic_snps']
        return 'bad snp'

    ## if looking at a subset of individuals, remove unwanted inds now
    ## important!  if using pop1 pop2 AND filter_inds, then args.population1, etc, must be in terms of the filtered list of individuals.
    if args.filter_inds:
        just_snp[:] = [just_snp[i] for i in range(len(just_snp)) if i in args.filter_inds]
        pass

    (refs, vars, haps, hets, chrs) = get_refs_vars_haps_hets(just_snp)
    #print ','.join([str(9+i) for i in args.filter_inds])
    #print chr, pos, (refs, vars, haps, hets)
    ret_data['refs'] = refs
    ret_data['vars'] = vars
    ret_data['haps'] = haps

    ## number of effective chrs, etc
    [num_chrs, max_missing_chrs, sex_chr] = calc_chr_missing_etc(chr, pos)
    #print sum([sum(t) for t in just_snp]), num_chrs


    ## should we skip this snp?
    if (refs == 0 and args.remove_fixed_vars) or vars == 0:
        # population is fixed for a snp
        if args.debug_2 or args.debug_read_snps: print "SKIP population is fixed for a snp", refs, vars, just_snp, snp
        return 'bad snp'
    if (num_chrs - haps) > int(max_missing_chrs + .000000001):
        # too much missing!
        if args.debug_2 or args.debug_read_snps: print "SKIP too much missing!", num_chrs - haps, '>', '%.50f' % max_missing_chrs, int(max_missing_chrs + .000000001), snp
        return 'bad snp'
    if hets == haps / 2:
        # totally heterozygous - potentially paralogous
        if args.debug_2 or args.debug_read_snps: print "SKIP totally heterozygous - potentially paralogous", snp
        return 'bad snp'
    if args.min_freq > vars / haps:
        if args.debug_2 or args.debug_read_snps: print "SKIP rejecting snp because of too low freq", snp
        return 'bad snp'
    if args.max_freq < vars / haps:
        if args.debug_2 or args.debug_read_snps: print "SKIP rejecting snp because of too high freq", snp
        return 'bad snp'

    ## get two population values (add here if new values are needed - currently just for fst)
    # if args.two_populations:
    #     pop1_snp = [just_snp[i] for i in range(len(just_snp)) if i in args.population1]
    #     pop2_snp = [just_snp[i] for i in range(len(just_snp)) if i in args.population2]
    #     (p1_refs, p1_vars, p1_haps, p1_hets, p1_chrs) = get_refs_vars_haps_hets(pop1_snp)
    #     (p2_refs, p2_vars, p2_haps, p2_hets, p2_chrs) = get_refs_vars_haps_hets(pop2_snp)
    #     ret_data['p1_refs'] = p1_refs
    #     ret_data['p1_vars'] = p1_vars
    #     ret_data['p1_haps'] = p1_haps
    #     ret_data['p2_refs'] = p2_refs
    #     ret_data['p2_vars'] = p2_vars
    #     ret_data['p2_haps'] = p2_haps
    #     print "seems to be messing up in the Y chromosome (see fst results?)"
    #     sys.exit(-1)
    #     if p1_chrs - p1_haps > (1-args.filter_missingness) * p1_chrs: return 'bad snp'
    #     if p2_chrs - p2_haps > (1-args.filter_missingness) * p2_chrs: return 'bad snp'
    #     pass

    if args.multiple_populations != None:
        ret_data['mpops'] = []
        for pop in args.multiple_populations:
            pop_data = {}
            ret_data['mpops'].append(pop_data)
            pop1_snp = [just_snp[i] for i in range(len(just_snp)) if i in pop]
            (p1_refs, p1_vars, p1_haps, p1_hets, p1_chrs) = get_refs_vars_haps_hets(pop1_snp)
            pop_data['refs'] = p1_refs
            pop_data['vars'] = p1_vars
            pop_data['haps'] = p1_haps
            pass
        pass


    # calculations that use ancestry
    if consider_ancestry:
        if consider_ancestry_bsg != None:
            ancestral_allele = args.ancestry_file_bsg.get_base_one_based(chr, pos)
        else:
            ancestry = ancestral_state_file.readline().strip().split('\t')
            ancestry[2] = int(ancestry[2])
            ancestral_allele = ancestry[ancestral_column]
            
            # check that chr and position match
            if chr != ancestry[1] or pos != ancestry[2] or ref != ancestry[3] or var != ancestry[4]:
                print "ERROR!!  chr and pos do not match:"
                print snp
                print ancestry
                sys.exit(-1)
                pass
            pass
        
        # don't take the upper of this, because if it's inferred through dnaml, we want to know if it's low confidence
        ret_data['anc'] = ancestral_allele
        
        ## check against the ancestral allele, to get total number of changes
        if ancestral_allele == ref:
            # if ancestral_allele matches ref
            mrca = 0
            ret_data['mrca_known'] = True
        elif ancestral_allele == var:
            # if ancestral_allele matches either var
            mrca = 1
            ret_data['mrca_known'] = True
        elif args.filt_retain_known_anc:
            if args.debug_2 or args.debug_read_snps: print "SKIP ancestry not known for this snp, or it doesn't match ref or var", chr, pos, ref, var, ancestral_allele
            return 'bad snp'
        else:
            # otherwise, take the major allele
            ret_data['mrca_known'] = False
            if vars > refs:
                mrca = 1
            else:
                mrca = 0
                pass
            pass
        ret_data['mrca'] = mrca

        ret_data['x'] = calc_x(ret_data)

        ## if we're filtering based on archaic AND ancestry, then only keep snps that are derived and present in the archaic (comp) bsg
        if args.filt_retain_arc and args.filt_retain_known_anc and args.comp and mrca == 0 and var not in ret_data['archaic_snps']:
            if args.debug_3 or args.debug_read_snps: print 'SKIP snp because var is derived, but is not found in the archaic chr:', chr, pos, ref, var, ret_data['archaic_snps']
            return 'bad snp'
        elif args.filt_retain_arc and args.filt_retain_known_anc and args.comp and mrca == 1 and ref not in ret_data['archaic_snps']:
            if args.debug_3 or args.debug_read_snps: print 'SKIP snp because ref is derived, but is not found in the archaic chr:', chr, pos, ref, var, ret_data['archaic_snps']
            return 'bad snp'
        
        pass

    ## why is this at the end?  I guess we have to be sure that we don't not-return for some other reason before we increment count?
    if args.limit != None and args.snp_count > args.limit:
        return 'eof'
    if args.report_progress != None and args.snp_count % args.report_progress == 0:
        sys.stderr.write('..at snp\t%d:\t%s\t%d;\t%.2f\n' % (args.snp_count, chr, pos, time.time() - args.check_time))
        args.check_time = time.time()
        pass
    args.snp_count += 1
    
    return ret_data

### decide which snp_helper should be called
if args.ms:
    read_snp_helper = read_snp_helper_ms
else:
    read_snp_helper = read_snp_helper_varfile
    pass


def read_chimp():
    s = read_chimp_helper()
    while s == 'bad snp':
        s = read_chimp_helper()
        pass
    if s == 'eof':
        return False
    return s

def read_chimp_helper():
    snp = chimp_snps_file.readline().strip().split('\t')
    if len(snp) == 1:
        return 'eof'
    if snp[3] == 'N':
        return 'bad snp'
    # see if this is in the regions file
    if args.chimp_mapped != None and not args.chimp_mapped.in_region_one_based(snp[0], int(snp[1])):
        return 'bad snp'

    ret_data = {}
    ret_data['chr'] = snp[0]
    ret_data['pos'] = int(snp[1])
    ret_data['hum'] = snp[2]
    ret_data['chimp'] = snp[3]
    return ret_data

def calc_ind_dist(snps, i1, i2, prt = False):
    if prt: print i1, i2
    for s in snps:
        if prt: print s['snp_details']
        pass
    d = [max(0, abs(s['snp_details'][i1][1] - s['snp_details'][i2][1]) - max(s['snp_details'][i1][0], s['snp_details'][i2][0])) for s in snps]
    if prt: print d
    return sum(d)

def calc_ind_dist_from_pop(snps, ind, pop, prt = False):
    if ind in pop:
        pop = list(pop)
        pop.remove(ind)
        pass
    if prt: print ind, pop
    for s in snps:
        if prt: print s['snp_details']
        pass
    d = [calc_ind_dist(snps, ind, p) for p in pop]
    if prt: print d
    return sum(d) / len(pop)


#############
## how does this handle missingness????

from snpdist import calc_dist_cy2_sum, calc_dist_cy

def calc_snp_dist(s1, s2, subset = None, prt = False, invert = False):
    if subset != None:
        #print s, subset
        s1 = [s1[i] for i in subset]
        s2 = [s2[i] for i in subset]
        pass

    # missing = snp[i][0]
    # refs = snp[i][1]
    # alts = snp[i][2]
    ## difference in ref alleles, minus max difference in missingness
    ##
    ## i.e., snp1, T = ref, C = alt
    ## i.e., snp2, T = ref, G = alt
    ##
    ## TT vs TG = max(0, abs(2-1) - max(0,0))
    ## TT vs TG = max(0, 1 - 0)
    ## TT vs TG = 1
    ##
    ## TC vs TG = max(0, abs(1-1) - max(0,0))
    ## TC vs TG = max(0, 0 - 0)
    ## TC vs TG = 0
    ##
    ## TC vs GG = max(0, abs(1-0) - max(0,0))
    ## TC vs GG = max(0, 1 - 0)
    ## TC vs GG = 1
    ##
    ## TT vs TT = max(0, abs(2-2) - max(0,0))
    ## TT vs TT = max(0, 0 - 0)
    ## TT vs TT = 0
    ##
    ## TT vs GG = max(0, abs(2-0) - max(0,0))
    ## TT vs GG = max(0, 2 - 0)
    ## TT vs GG = 2
    ##
    ## TT vs Tx = max(0, abs(2-1) - max(0,1))
    ## TT vs Tx = max(0, 1 - 1)
    ## TT vs Tx = 0
    ##
    ## TT vs xx = max(0, abs(2-0) - max(0,2))
    ## TT vs xx = max(0, 2 - 2)
    ## TT vs xx = 0
    ##
    ## this is an issue, but only for CGI data (that hasn't been filtered to set missingness as homozygous)
    ## Tx vs Gx = max(0, abs(1-0) - max(1,1))
    ## TT vs xx = max(0, 1 - 1)
    ## TT vs xx = 0
    if not args.fpc:
        print "snpdist is broken for snps with partial calls!"
        sys.exit(-1)
        pass

    #d = sum([max(0, abs(s1[i][1] - s2[i][1]) - max(s1[i][0], s2[i][0])) for i in range(len(s1))])
    #d = sum(calc_dist_cy(s1, s2))
    d = calc_dist_cy2_sum(s1, s2)
    #print s1
    #print s2
    #print d
    return d

def print_snps(snps, mark = None, mark_inds = None):
    print
    if mark != None:
        print ' '.join(['  '] + ['*' if s in mark else ' ' for s in range(len(snps))])
        pass
    for i in range(len(individuals)):
        print ' '.join(['i%d' % i] + [str(snps[s]['snp_details'][i][2]) for s in range(len(snps))] + ['*' if mark_inds != None and i in mark_inds else ' '])
        pass
    return

def get_nonvarcounts1(s1, s2, individuals, varpop):
    return [s1['snp_details'][i][2] + s2['snp_details'][i][2] for i in list(set(range(len(individuals))) - set(varpop))]
def get_nonvarcounts1a(s1, s2, individuals, varpop):
    return [s1['snp_details'][i][2] for i in list(set(range(len(individuals))) - set(varpop))]
def get_nonvarcounts1b(s1, s2, individuals, varpop):
    return [s2['snp_details'][i][2] for i in list(set(range(len(individuals))) - set(varpop))]
def get_nonvarcounts2(s1, s2, individuals, varpop):
    return [s1[i] + s2[i] for i in list(set(range(len(individuals))) - set(varpop))]
def get_nonvarcounts3(s1, s2, individuals, varpop):
    return [i for i in range(len(individuals)) if i not in varpop]
def get_nonvarcounts4(s1, s2, individuals, varpop):
    return list(set(range(len(individuals))) - set(varpop))


# doesn't really look for snps in phase - looks for snps that are non-ref in the same individuals
# take two snps (s1, s2)
def calc_S(s1, s2, subset = None, prt = False, varpop = None, allowable_diffs = 5, cache = None, max_diffs_for_perfect_score = 0):

    if cache != None and cache[s1['pos']][s2['pos']] != None:
        return cache[s1['pos']][s2['pos']]

    def store_in_cache(val):
        if cache == None: return
        cache[s1['pos']][s2['pos']] = val
        pass

    # guarantee that snps are at least 10bp apart
    if abs(s1['pos'] - s2['pos']) < 10:
        if prt: print 'too close!', s1, s2
        store_in_cache(-sys.maxint)
        return -sys.maxint

    #s1v = [s[2] for s in s1['snp_details']]
    #s2v = [s[2] for s in s2['snp_details']]


    # if we're requiring snps to be nonref only in varpop..
    ## (this has been moved out to the code that calls calc_s_star)
    if varpop != None:
        print "cannot currently do varpop in calc_S!  This has been moved to the code that calls calc_s_star."
        sys.exit(-1)
        pass

    # actually calculate the distance between these snps (potentially on a subset of individuals)
    # this could use varpop (instead of subset), if varpop != None; we've already checked to make sure that these snps don't exist in !varpop.
    # It would be a bit faster, but not too much faster if !varpop is small 
    d = calc_snp_dist(s1['snp_details'], s2['snp_details'], subset)
    if prt: print 'dist:', d

    # don't allow this if d > 5 or if there are no variants in either of these snps
    if subset == None: subset = range(len(individuals))
    s1_subset_vars = sum([s1['snp_details'][i][2] for i in subset])
    s2_subset_vars = sum([s2['snp_details'][i][2] for i in subset])
    #s1_subset_vars = 1
    #s2_subset_vars = 1
    if args.debug_2: print 'subset var counts:', s1_subset_vars, s2_subset_vars

    if args.debug_sstar: print 'debug sstar: ', s1['pos'], s2['pos'], d, allowable_diffs, max_diffs_for_perfect_score, s1_subset_vars, s2_subset_vars

    if d > max(allowable_diffs, max_diffs_for_perfect_score) or s1_subset_vars == 0 or s2_subset_vars == 0:
        store_in_cache(-sys.maxint)
        if args.debug_sstar: print 'debug sstar: too many diffs, returning maxint. %d > max(%d, %d)\n' % (d, allowable_diffs, max_diffs_for_perfect_score)
        return -sys.maxint

    # perfect!
    if d <= max_diffs_for_perfect_score:
        val = abs(s1['pos'] - s2['pos']) + args.perfect_score_bonus
        store_in_cache(val)
        if args.debug_sstar: print 'debug sstar: perfect match. %d <= %d\n' % (d, max_diffs_for_perfect_score)
        return val

    if args.debug_sstar: print 'debug sstar: slightly imperect with %d diffs, charging %d\n' % (d, args.imperfect_charge)

    # slightly imperfect
    imperfect_charge = args.imperfect_charge
    store_in_cache(imperfect_charge)
    return imperfect_charge

def get_common_in_sets(sets):
    common = set()
    [common.update(s) for s in sets]
    [common.intersection_update(s) for s in sets]
    return sorted(list(common))

def get_majority_in_sets(sets):
    counts = {}
    [counts.setdefault(i, 0) for i in range(num_genotypes)]
    for s in sets:
        for i in s:
            counts[i] += 1
            pass
        pass
    return sorted([i for i in counts if counts[i] >= len(sets) / 2])

def remove_common_in_sets(large_set, sets):
    common = set()
    [common.update(s) for s in sets]
    [common.intersection_update(s) for s in sets]
    return sorted(list(set(large_set) - common))

def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))


# chr1 829000 839000 - test window
def calc_S_star(snps, subset = None, prt = False, ret_inds = False, varpop = None, allowable_diffs = 5, max_diffs_for_perfect_score = 0, force_ind = None, snp_indices = None, cache = None):

    ## these are the same length as snps, not snp_indices, so that we can just use the snp index as the index in these as well
    ## - wasted space, but shouldn't be too bad in reasonably sized regions
    S_star = [0] * len(snps)
    S_snps = [[]] * len(snps)

    if prt: print 'S_star', S_star
    if prt: print 'S_snps', S_snps

    # initialize snp_indices (can be used to restrict the set of snps that are considered)
    if snp_indices == None: snp_indices = range(len(snps))

    # if we're requiring snps to be found in a particular individual..
    if force_ind != None:
        # if subset != None:
        #     print "can't do force_ind and subset!"
        #     sys.exit(-1)
        #     return

        # get the list of snps in the force_ind individual
        force_ind_indices = [i for i,s in enumerate(snps) if s['snp_details'][force_ind][2] > 0]
        # and then get the intersection of those snps and the snps we're considering (could make this faster by just considering the snp_indices to begin with?)
        snp_indices = sorted(set(snp_indices).intersection(set(force_ind_indices)))
        if args.debug: print 'force-ind is %d, vars: %s' % (force_ind, snp_indices)
        pass
    
    if prt: print 'force inds snps', snp_indices
    if prt: print 'force inds snps', [snps[i]['pos'] for i in snp_indices]

    if args.print_sstar_data_for_revisions_on_neand_paper:
        print
        print "SNP dist table for individual", force_ind
        print '.', ' '.join(['snp_%d' % s for s in snp_indices])
        for i in snp_indices:
            dists = [str(calc_snp_dist(snps[i]['snp_details'], snps[j]['snp_details'], subset)) if i != j else '-' for j in snp_indices]
            print "snp_%d" % i, ' '.join(dists)
            pass
        print
        print "S* table for individual", force_ind
        print '.', ' '.join(['snp_%d' % s for s in snp_indices])
        pass

    # look at every snp k and j (no subsetting here)
    for j in snp_indices:
        if j == snp_indices[0]:
            S_snps[snp_indices[0]] = [snp_indices[0]]
            if prt: print 'setting first index', snp_indices[0], S_snps[snp_indices[0]], S_snps
            continue
        # current S_star for snp j, if snp k is the previously chosen snp:
        #  i.e., x_x___k___j ('x' snps do not matter - their cost is included in S_star[k])
        S_tmp = [S_star[k] + calc_S(snps[k], snps[j], subset = subset, varpop = varpop, allowable_diffs = allowable_diffs, max_diffs_for_perfect_score = max_diffs_for_perfect_score, cache = cache) \
                     for k in [x for x in snp_indices if x < j]]
        S_tmp2 = [(S_star[k], calc_S(snps[k], snps[j], subset = subset, varpop = varpop, allowable_diffs = allowable_diffs, max_diffs_for_perfect_score = max_diffs_for_perfect_score, cache = cache)) \
                     for k in [x for x in snp_indices if x < j]]
        #print S_tmp
        #print S_tmp2

        S_max = max(S_tmp)


        # if S_max is negative, then store 0 and just [j] as the snp - allow the algorithm to start over from here?
        if S_max > 0:
            # just gets the first (and thus farthest away) k that maximizes S_star
            # there could be other k, but I don't know how often that occurs
            #S_index = S_tmp.index(S_max)
            S_index = snp_indices[S_tmp.index(S_max)]
            S_star[j] = S_max
            S_snps[j] = S_snps[S_index] + [j]
            if prt: print j, 'S_max, S_index', S_max, S_index#, snp_indices[S_index], S_snps[snp_indices[S_index]]
        else:
            S_star[j] = 0
            S_snps[j] = [j]
            pass

        if args.print_sstar_data_for_revisions_on_neand_paper:
            print 'snp_%d' % j, ' '.join([str(i) if i >= 0 else '-inf' if i == -sys.maxint else '-inf+%d' % (i + sys.maxint) for i in S_tmp] + ['-'] * (len(snp_indices) - snp_indices.index(j)))
            print 'snp_%d' % j, ' '.join([(str(i) if i >= 0 else '-inf' if i == -sys.maxint else '-inf+%d' % (i + sys.maxint)) + '+' + str('-inf' if x == -sys.maxint else x) for i,x in S_tmp2] + ['-'] * (len(snp_indices) - snp_indices.index(j)))
            pass

        # print S values for all snps (k) for this snp (j)
        if prt: print j, 'S_tmp', S_tmp
        if prt: print j, 'S_star', S_star
        if prt: print j, 'S_snps', S_snps
        if prt: print
        pass

    S_snps_optimal = S_snps[S_star.index(max(S_star))] if len(snps) > 1 else []

    if len(snps) > 1 and len(S_snps_optimal) > 1:
        if prt: print 'result:', max(S_star), S_snps_optimal
        if ret_inds:
            # need to account for subset
            if subset == None: subset = range(len(individuals))
            # for every snp in the optimal solution, look at each individual i in subset.
            # if i has a non-ref allele at that snp, then include i in the list of inds for that snp
            S_inds = [[i for i in subset if snps[s]['snp_details'][i][2] > 0] for s in S_snps_optimal]
            return (max(S_star), S_inds, S_snps_optimal)
        return max(S_star)
    if ret_inds:
        return(0, [], [])
    return 0


def get_nonvarcount_sum(s, individuals, varpop):
    ## 2 is vars: (missing, ref, var)
    return sum([s['snp_details'][i][2] for i in list(set(range(len(individuals))) - set(varpop))])
def get_varcount_sum(s, individuals, varpop):
    ## 2 is vars: (missing, ref, var)
    return sum([s['snp_details'][i][2] for i in varpop])

primes = (13,17,19,23,29,31)

def calc_region_s_star(data):

    ## if there's no data in this window, just forget about it
    if len(data['snps']) == 0:
        return

    if args.print_sstar_data_for_revisions_on_neand_paper:
        print "snp pos", ' '.join(['ind_%d' % i for i in args.target_pops]), ' '.join(['ind_%d' % i for i in args.non_target_pop])
        for i,snp in enumerate(data['snps']):
            print "snp_%d" % i, snp['pos'], ' '.join([str(snp['snp_details'][i][2]) for i in args.target_pops]), ' '.join([str(snp['snp_details'][i][2]) for i in args.non_target_pop])
            #print "snp_%d" % i, snp['pos'], ' '.join([str(s[1]) for s in snp['snp_details']])
            pass
        pass
    
    summary_stat_hash = {}
    summary_stat_hash['region_length'] = data['region_length']
    #summary_stat_hash['region_tag'] = data['region_tag'] + ' ' + data['region_tag2']
    summary_stat_hash['allowed_diffs'] = args.allowable_diffs
    summary_stat_hash['max_diffs'] = args.max_diffs_for_perfect_score
    summary_stat_hash['perfect_bonus'] = args.perfect_score_bonus
    summary_stat_hash['imperfect_penalty'] = args.imperfect_charge
    summary_stat_hash['missingness_filt'] = args.filter_missingness
    [summary_stat_hash.setdefault(s, 'NA') for s in ('S_star_varpop', 'varpop_common', 'varpop_maj', 'count_varpop_common', 'count_varpop_maj', 'num_snps', 'num_inds', 'haplen_varpop', 'reglen', 'ind_num', 'all_inds', 'num_snps_all', 'num_snps_all_subpop', 'num_snps_varpop', 'hap_start', 'hap_end')]
    if args.debug: supposed_summary_stat_hash_keys = set(summary_stat_hash)

    ## here we get a list of varpop-only indices, and send that to calc_S_star
    ## and then draw j and k from just those
    ## the problem is that it doesn't allow us to print out info about non-varpop snps
    varpop_indices = [i for i,s in enumerate(data['snps']) if get_nonvarcount_sum(s, individuals, args.target_pops) == 0]
    non_varpop_indices = [i for i,s in enumerate(data['snps']) if get_nonvarcount_sum(s, individuals, args.target_pops) > 0]

    print [get_nonvarcount_sum(s, individuals, args.target_pops) for i,s in enumerate(data['snps'])]
    
    #if True or debug_ms: print 'varpop indices', args.ms_chr, varpop_indices
    #if True or debug_ms: print 'non_varpop indices', args.ms_chr, non_varpop_indices

    if args.subset_inds_into_groups == None:
        run_calc_region_s_star(data, summary_stat_hash, None, varpop_indices, non_varpop_indices)
        return
    
    def mygrouper(iterable, n):
        args = [iter(iterable)] * n
        return ([e for e in t if e != None] for t in itertools.izip_longest(*args))
    random_groups_of_inds = list(mygrouper(random.sample(args.target_pops, len(args.target_pops)), args.subset_inds_into_groups))
    if args.subset_not_random: random_groups_of_inds = list(mygrouper(args.target_pops, args.subset_inds_into_groups))

    

    for subset_count, subset in enumerate(random_groups_of_inds):
        summary_stat_hash['Subset'] = subset_count + 1
        run_calc_region_s_star(data, summary_stat_hash, subset, varpop_indices, non_varpop_indices, subset_count + 1)
        pass
    
    return

def run_calc_region_s_star(data, summary_stat_hash, subset, varpop_indices, non_varpop_indices, subset_count = None):

    if args.forward_ms_file:
        if args.debug: print "printing ms file for", data['window_chr']
        if args.debug: print "currently", len(ms_file_backlog), "ms files in the backlog:", ms_file_backlog.keys()
        sys.stdout.write( ms_file_backlog[data['window_chr']].getvalue() )
        del ms_file_backlog[data['window_chr']]
        pass

    supposed_summary_stat_hash_len = len(summary_stat_hash)

    snps = data['snps']
    target_pop = args.target_pops

    ## also want to report the number of non-filtered snps (when we're
    ## doing something like -x yri.bbg to exclude yoruban snps, these
    ## should be the same - but they are not the same when we're using subset)
    summary_stat_hash['num_snps_all'] = len(snps)
    summary_stat_hash['num_snps_varpop'] = len(varpop_indices)
    summary_stat_hash['num_inds'] = len(target_pop)

    window_chr = data['window_chr']
    window_start = data['window_start']
    window_end = data['window_end']
    window_label = '%s.%d.%d' % (window_chr, window_start, window_end)


    if args.debug: print 'varpop positions pre subset:', [snps[i]['pos'] for i in varpop_indices]


    ## if we're just looking at a subset, update target_pop and varpop_indices
    ## we have to do this after initializing varpop indices, because we don't want to treat target-population individuals that just
    ## are not in the subset as non-target-pop when initializing varpop indices
    if subset != None:
        subpop_indices = set([i for i,s in enumerate(snps) if sum([s['snp_details'][j][2] for j in subset]) > 0])
        varpop_indices = sorted(list(subpop_indices.intersection(set(varpop_indices))))
        target_pop = subset

        ## we have to make sure that there aren't any fixed variants in this new subset - mostly for consistancy, but also because in the 
        ## ms simulations we are just generating sequence for 20 inds, and are definitely removing fixed variants there
        ## basically, we just require that those variants have *some* reference site in the subpop inds (we already know they have some var site)
        varpop_indices = [i for i in varpop_indices if sum([snps[i]['snp_details'][j][1] for j in subset]) > 0]
        pass

    if args.debug: print 'target_pop:', target_pop
    if args.debug: print 'varpop positions:', [snps[i]['pos'] for i in varpop_indices]
    if debug_ms: print 'varpop target_pop:', target_pop
    if debug_ms: print 'varpop positions:', [snps[i]['pos'] for i in varpop_indices]
    if debug_ms: print 'varpop freqs:', [snps[i]['haps'] for i in varpop_indices]


    ## if there aren't any snps for this subset, return
    if len(varpop_indices) == 0: return
    

    ## restrict the number of pop specific snps (this is mostly used for simulated data, to make it look like a particular region and get a null value)
    if args.pop_specific_snp_count != None:
        if len(varpop_indices) < args.pop_specific_snp_count:
            print "too few snps: %d vs %d" % (len(varpop_indices), args.pop_specific_snp_count)
            #print varpop_indices
            return
        varpop_indices = random.sample(varpop_indices, args.pop_specific_snp_count)
        pass

    ## this is the number of snps after all filtering - this could be less than varpop indices, because of subsetting or pop_specific_snp_count
    summary_stat_hash['num_snps'] = len(varpop_indices)
    summary_stat_hash['num_snps_all_subpop'] = len(set(varpop_indices).union(set(non_varpop_indices)))


    ## get branch length for outgroup (only works with ms)
    if args.ms and args.ms_outgroup != None:
        ms_outgroup_bl_tot = args.ms_outgroup_branch_len[window_chr]
        hash = '%d-%d' % (window_start, window_end)
        if hash in args.ms_outgroup_branch_len_by_win[window_chr]:
            ms_outgroup_bl = args.ms_outgroup_branch_len_by_win[window_chr][hash]
            ms_outgroup_bl_by_chr = args.ms_outgroup_branch_len_by_win_by_chr[window_chr][hash]
        else:
            ms_outgroup_bl = 0
            ms_outgroup_bl_by_chr = mydefaultdict(lambda : 0)
            pass
    else:
        ms_outgroup_bl_tot = 'NA'
        ms_outgroup_bl = 'NA'
        ms_outgroup_bl_by_chr = mydefaultdict(lambda : 0)
        pass


    if not args.summary_only:
        win_out = cStringIO.StringIO()
        pass
    
    window_stat_hash = {}
    [window_stat_hash.setdefault(s, 'NA') for s in ('S_star_varpop', 'ind_num', 'num_snps_ind', 'num_snps_tag', 'hap_start', 'hap_end')]
    if args.ms:
        window_stat_hash['ms_c1'] = 'NA'
        window_stat_hash['ms_c2'] = 'NA'
        window_stat_hash['ms_c1_tag_snp_count'] = 'NA'
        window_stat_hash['ms_c2_tag_snp_count'] = 'NA'
        window_stat_hash['ms_region'] = window_chr
        pass
    if args.debug: supposed_window_stat_hash_keys = set(window_stat_hash)
    supposed_window_stat_hash_len = len(window_stat_hash)

    ## this is just a way to store the region information in the header - kind of janky
    window_stat_hash['ind_num'] = len(target_pop)
    window_stat_hash['num_snps_ind'] = len(varpop_indices)
    window_stat_hash['S_star_varpop'] = window_chr
    window_stat_hash['hap_start'] = window_start
    window_stat_hash['hap_end'] = window_end


    
    if not args.summary_only:
        ## print the header and first informational line to the window file
        if args.no_snp_report:
            win_out_snp_header = []
            win_out_snp_info = []
        else:
            win_out_snp_header = ['snp%d' % (i+1) for i in range(len(varpop_indices))]
            win_out_snp_info = [str(snps[s]['pos']) for s in varpop_indices]
            pass
        win_out.write('\t'.join(sorted(window_stat_hash.keys()) + ['ms_outgroup_bl', 'ms_outgroup_bl_tot'] + win_out_snp_header) + '\n')
        win_out.write('\t'.join([str(k) for k in [window_stat_hash[i] for i in sorted(window_stat_hash.keys())] + 
                                 [ms_outgroup_bl, ms_outgroup_bl_tot]] + win_out_snp_info) + '\n')
        pass




    ## set up the cache
    ## this is a redundant matrix (the dist between i and j is the same as j and i)

    if 'cache' not in data['persist']:
        data['persist']['cache'] = {}
        pass
    ## hack to keep the cache from growing out of control
    if len(data['persist']['cache'].keys()) > 5000: data['persist']['cache'] = {}
    ## also reset the cache if we're doing subsetting - you can't reuse values of S if they've been computed on different sets of individuals!
    if args.subset_inds_into_groups != None: data['persist']['cache'] = {}
    cache = data['persist']['cache']
    for i in varpop_indices:
        i_pos = snps[i]['pos']
        if i_pos not in cache: cache[i_pos] = {}
        for j in varpop_indices:
            j_pos = snps[j]['pos']
            if j_pos not in cache[i_pos]: cache[i_pos][j_pos] = None
            pass
        pass

    #### old cache code
    # cache = {}
    # for i in varpop_indices:
    #     i_pos = snps[i]['pos']
    #     cache[i_pos] = {}
    #     for j in varpop_indices:
    #         j_pos = snps[j]['pos']
    #         cache[i_pos][j_pos] = None
    #         pass
    #     pass

    # sys.stderr.write('cache size: %d\n' % (len(varpop_indices) * len(varpop_indices)))
    max_sstar = 0
    pearson_cache = {}

    ## loop through every individual in target_pop
    for current_ind in target_pop:
        
        # ('S_star_varpop', 'varpop_common', 'varpop_maj', 'count_varpop_common', 'count_varpop_maj', 'num_snps', 'num_inds', 'reglen', 'ind_num', 'all_inds')

        if args.keep_inds != None:
            ## HAVE TO SUBTRACT ONE BECAUSE TARGET_POP IS ZERO BASED, BUT KEEP_INDS IS ONE BASED!
            window_stat_hash['ind_num'] = args.keep_inds[current_ind] - 1
        else:
            window_stat_hash['ind_num'] = current_ind
            pass
        #print "current ind", current_ind, window_stat_hash['ind_num']
        # window_stat_hash['ms_chr1'] = args.ms_chr_pairs[current_ind][0]
        # window_stat_hash['ms_chr2'] = args.ms_chr_pairs[current_ind][1]
        if args.ms:
            c1 = args.ms_chr_pairs[current_ind][0]
            c2 = args.ms_chr_pairs[current_ind][1]
            window_stat_hash['ms_c1'] = c1
            window_stat_hash['ms_c2'] = c2
            intr_for_ind = max(ms_outgroup_bl_by_chr[c1], ms_outgroup_bl_by_chr[c2])
            # intr_for_ind = max(args.ms_archaic_regions[c1].count(), args.ms_archaic_regions[c2].count())
            # print 'intr_for_ind1', intr_for_ind
            # print 'intr_for_ind2', intr_for_ind
        else:
            intr_for_ind = 'NA'
            pass


        # calculate various versions of S_star
        if args.debug: print 'about to calc sstar'
        ad = args.allowable_diffs
        md = args.max_diffs_for_perfect_score

        #####################
        ## TO DO SUBSET, WOULDN'T IT MAKE MORE SENSE TO FITLER THE SET OF SNPS ONCE (OUT HERE), RATHER THAN FILTER THEM WITHIN CALC_SNP_DIST MANY TIMES?  JUST LIKE WE DO FOR VARPOP, CURRENTLY.
        #####################

        [S_star_varpop, S_inds_varpop, S_snps_varpop] = calc_S_star(snps, ret_inds = True, prt=args.debug, force_ind = current_ind, snp_indices = varpop_indices, cache = cache, allowable_diffs = ad, max_diffs_for_perfect_score = md, subset = subset)
        # [S_star_varpop, S_inds_varpop, S_snps_varpop] = calc_S_star(snps, ret_inds = True, prt=args.debug, force_ind = current_ind, snp_indices = varpop_indices, cache = cache, allowable_diffs = args.allowable_diffs, max_diffs_for_perfect_score = args.max_diffs_for_perfect_score)
        # [S_star_varpop, S_inds_varpop, S_snps_varpop] = calc_S_star(snps, ret_inds = True, prt=args.debug, force_ind = current_ind, snp_indices = varpop_indices, cache = cache, allowable_diffs = 1000, max_diffs_for_perfect_score = 5)
        # print args.allowable_diffs, args.max_diffs_for_perfect_score
        if args.debug: print 'calced!'


        window_stat_hash['S_star_varpop'] = S_star_varpop
        ############# i have no idea why we used to use snps instead of varpop_indices - I'm actually not sure that this is ever used downstream
        window_stat_hash['num_snps_ind'] = sum([snps[s]['snp_details'][current_ind][2] > 0 for s in varpop_indices])
        # window_stat_hash['num_snps_ind'] = sum([snp['snp_details'][current_ind][2] > 0 for snp in snps])
        window_stat_hash['num_snps_tag'] = len(S_snps_varpop)
        if args.debug: print 'S varpop:', S_star_varpop, S_snps_varpop, S_inds_varpop, get_common_in_sets(S_inds_varpop)
        if args.ms:
            ms_chrs_for_tag_snps = [snps[s]['snp_details_ms_haplotypes'][current_ind] for s in S_snps_varpop]
            ms_chrs_for_tag_snps_flat = [item for sublist in ms_chrs_for_tag_snps for item in sublist]
            c1 = args.ms_chr_pairs[current_ind][0]
            c2 = args.ms_chr_pairs[current_ind][1]
            ms_chrs_for_tag_snps_counter = Counter(ms_chrs_for_tag_snps_flat)
            c1_c = ms_chrs_for_tag_snps_counter[c1] if c1 in ms_chrs_for_tag_snps_counter else 0
            c2_c = ms_chrs_for_tag_snps_counter[c2] if c2 in ms_chrs_for_tag_snps_counter else 0
            window_stat_hash['ms_c1_tag_snp_count'] = c1_c
            window_stat_hash['ms_c2_tag_snp_count'] = c2_c
            if args.debug: print 'S varpop ms chrs:', len(S_snps_varpop), ms_chrs_for_tag_snps, ms_chrs_for_tag_snps_flat, ms_chrs_for_tag_snps_counter, c1_c, c2_c
            pass

        # get involved individuals
        #vp_com_inds = get_common_in_sets(S_inds_varpop)
        #vp_maj_inds = get_majority_in_sets(S_inds_varpop)
        #window_stat_hash['varpop_common'] = ' '.join([str(s) for s in vp_com_inds])
        #window_stat_hash['varpop_maj'] = ' '.join([str(s) for s in vp_maj_inds])
        #window_stat_hash['count_varpop_common'] = len(vp_com_inds)
        #window_stat_hash['count_varpop_maj'] = len(vp_maj_inds)
        
        # get snps for involved individuals
        # vp_com_snps = [s for i,s in enumerate(S_snps_varpop) if len(set(S_inds_varpop[i]).intersection(set(vp_com_inds))) != 0]
        # window_stat_hash['reglen'] = snps[vp_com_snps[-1]]['pos'] - snps[vp_com_snps[0]]['pos'] if len(vp_com_snps) > 1 else 0

        ## get variants in high LD with the tag variants
        # print current_ind, S_star_varpop, 'varpop indices!', varpop_indices
        # print current_ind, S_star_varpop, 'sstar snps!!', S_snps_varpop
        # print current_ind, S_star_varpop, 'reduced set: ', len(S_snps_varpop), len(varpop_indices), len(non_tag_varpop_indices)
        pear_snps = []
        if len(S_snps_varpop) > 0 and not args.summary_only and args.calc_pearson:
            non_tag_varpop_indices = [v for v in varpop_indices if v not in S_snps_varpop]
            for s in non_tag_varpop_indices:
                pear = [calc_pearson_freal(snps, t, s, subset = target_pop if len(target_pop) != num_genotypes else None, \
                                               cache = pearson_cache, \
                                               check_missingness = not args.varfile_is_1kg) for t in S_snps_varpop]
                # print 'pearson', s, sum([p > .9 for p in pear]), pear
                if sum(pear) / len(pear) > .9: pear_snps.append(s)
                pass
            
            ## get limits of haplotype
            window_stat_hash['hap_start'] = snps[min(pear_snps + S_snps_varpop)]['pos']
            window_stat_hash['hap_end'] = snps[max(pear_snps + S_snps_varpop)]['pos']
            pass
        elif len(S_snps_varpop) > 0:
            ## get limits of haplotype
            window_stat_hash['hap_start'] = snps[min(S_snps_varpop)]['pos']
            window_stat_hash['hap_end'] = snps[max(S_snps_varpop)]['pos']
        else:
            ## set fake limits of haplotype
            window_stat_hash['hap_start'] = 'NA'
            window_stat_hash['hap_end'] = 'NA'
            pass
        
        # print 'final pearson', pear_snps


        if args.debug: print 'varpop_common inds:', S_inds_varpop, S_star_varpop
        if args.debug: print 'varpop snps:', len(snps), len(S_snps_varpop), S_snps_varpop

        #### Save max values for several stats
        if S_star_varpop > max_sstar: max_sstar = S_star_varpop

        #### Report
        if not args.summary_only:
            
            if len(window_stat_hash) != supposed_window_stat_hash_len:
                print "we have incorrect stat hash!"
                print sorted(window_stat_hash)
                print [window_stat_hash[k] for k in sorted(window_stat_hash)]
                print len(window_stat_hash), supposed_window_stat_hash_len
                if args.debug: print set(window_stat_hash).difference(supposed_window_stat_hash_keys)
                sys.exit(-1)
                return
            
            report = '\t'.join([str(k) for k in [window_stat_hash[i] for i in sorted(window_stat_hash.keys())] + [intr_for_ind, ms_outgroup_bl_tot]])
            #present_in_nontarget = [sum([s['snp_details'][j][2] > 0 for j in args.non_target_pop]) > 0 for s in snps]
            
            def calc_arc_codes(snp):
                code = 1
                var_present = snp['snp_details'][current_ind][2] > 0
                for i,arc in enumerate(snp['is_archaic']):
                    if arc and var_present: code *= primes[i]
                    pass
                return code

            if args.no_snp_report:
                snp_report = []
            else:
                snp_report = ['%d' % ((2 if s in S_snps_varpop else 1) * \
                                          # two if tag snp \
                                          # three if het \
                                          (3 if snps[s]['snp_details'][current_ind][2] == 1 else 1) * \
                                          # five if it's hom \
                                          (5 if snps[s]['snp_details'][current_ind][2] == 2 else 1) * \
                                          # seven if it's an introgressed snp \
                                          (7 if args.ms and args.ms_arc and snps[s]['archaic_inds'][current_ind] else 1) * \
                                          # eleven if this snp is in high LD with the tag snps \
                                          (11 if s in pear_snps else 1) * \
                                          # thirteen if this snp matches a variant found in the legit archaic chr \
                                          # (remember that the chr number in ms_outgroup is 1-based, but everything else is 0-based) \
                                          #(13 if snps[s]['is_archaic'] and snps[s]['snp_details'][current_ind][2] > 0 else 1) * \
                                          (calc_arc_codes(snps[s]) if (args.ms and args.ms_arc) or not args.ms else 1) * \
                                          1 ) for s in varpop_indices]
                # 1 ) for s in range(len(snps))]
                pass
            report = '\t'.join([report] + snp_report) + '\n'
            win_out.write(report)
            
            ## this option is only used for the very strange circumstance of printing the S* value for every ms chromosome
            ## so that they can be later interpreted by get_pct_arc_per_ind_from_ms_file_new.py to annotate results with the S* threshold required to produce those results
            if args.forward_ms_chr_sstar:
                if len(args.keep_pops) != 2 or args.keep_pops[0] != '1' or args.keep_pops[1] not in ('2', '3'):
                    print "We must be keeping pops 1 and 2 or 1 and 3 for forward_ms_chr_sstar to work.", args.keep_pops
                    sys.exit(-1)
                    pass
                if window_stat_hash['ms_c1'] < args.ms_pops[0]:
                    print "Must have pop 1 be the reference pop for forward_ms_chr_sstar to work.", args.keep_pops
                    sys.exit(-1)
                    pass
                ms_chr_incrementer = 1 if args.keep_pops[1] == '2' else 1 + args.ms_pops[1] if args.keep_pops[1] == '3' else None
                print 'ms_chr_sstar', window_stat_hash['ms_region'], window_stat_hash['S_star_varpop'], \
                    window_stat_hash['ms_c1'] + ms_chr_incrementer, window_stat_hash['hap_start'], window_stat_hash['hap_end'], \
                    summary_stat_hash['num_snps_all'], summary_stat_hash['num_snps_all_subpop']
                print 'ms_chr_sstar', window_stat_hash['ms_region'], window_stat_hash['S_star_varpop'], \
                    window_stat_hash['ms_c2'] + ms_chr_incrementer, window_stat_hash['hap_start'], window_stat_hash['hap_end'], \
                    summary_stat_hash['num_snps_all'], summary_stat_hash['num_snps_all_subpop']
                pass
            pass

        pass

    summary_stat_hash['S_star_varpop'] = max_sstar

    ## print a line to the summary file
    if args.first_line and not args.suppress_summary:
        print '\t'.join(['window_chr', 'window_start', 'window_end'] \
                            + sorted(summary_stat_hash.keys()) + \
                            ['ms_outgroup_bl', 'ms_outgroup_bl_tot'])
        args.first_line = False
        pass
    report = '\t'.join([str(k) for k in [window_chr, window_start, window_end] + [summary_stat_hash[i] for i in sorted(summary_stat_hash.keys())] + [ms_outgroup_bl, ms_outgroup_bl_tot]])
    if not args.suppress_summary: print report
    sys.stdout.flush()

    ## don't save the detailed window file if sstar is zero.
    ## could try to check sstar before writing the file at all, but then we'd have to go back through all of the individuals
    if max_sstar == 0 or args.summary_only:
        # do nothing
        pass
    # elif args.window_file_tmp_output_dir != None:
    #     ## set up the directory structure for output files, if necessary
    #     if args.window_file_tmp_output_dir != None and not os.path.isdir(args.window_file_tmp_output_dir + '/' + data['window_chr']):
    #         os.makedirs(args.window_file_tmp_output_dir + '/' + data['window_chr'])
    #         pass
    #     win_ofile = open('%s/%s/%s.%d-%d' % (args.window_file_tmp_output_dir, window_chr, window_chr, window_start, window_end), 'w')
    #     win_ofile.write(win_out.getvalue())
    #     win_ofile.close()
    #     win_out.close()
    #     final_file_name = '%s/%s/%s.%d-%d' % (args.window_file_output_dir, window_chr, window_chr, window_start, window_end)
    #     ## only write to tmp and then move so that we can implement buffering, if necessary - probably should just get rid of the tmp file functionality altogether
    #     shutil.move(win_ofile.name, final_file_name)
    elif args.window_file_output_dir != None:
        ## set up the directory structure for output files, if necessary
        window_file_odir = args.window_file_output_dir
        window_file_odir += '' if subset == None else ('.%d' % subset_count)
        window_file_odir +=  '/' + data['window_chr']
        if not os.path.isdir(window_file_odir):
            os.makedirs(window_file_odir)
            pass
        win_ofile = open('%s/%s.%d-%d' % (window_file_odir, window_chr, window_start, window_end), 'w')
        win_ofile.write(win_out.getvalue())
        win_ofile.close()
        win_out.close()
    elif args.window_file_append != None:
        args.window_file_append.write(win_out.getvalue())
        win_out.close()
    else:
        #sys.stderr.write("Should provide either -wf-odir or -wf-append!\n")
        pass
    
    
    return



def calc_tmrca_subsets(data):

    stat_hash = {}
    stat_hash['region_length'] = data['region_length']
    [stat_hash.setdefault(s, 'NA') for s in ('S_star_varpop', 'S_star_varpop_ad0', 'tmrca_all', 'tmrca_target_pop', 'tmrca_varpop', 'tmrca_varpop_maj', \
                                                 'varpop_common', 'varpop_maj', 'count_varpop_common', 'count_varpop_maj', 'tmrca_non_target_pop', 'tmrca_varpop_ad0', \
                                                 'varpop_ad0_common', 'count_varpop_ad0_common', 'num_snps', 'haplen_varpop', 'target_inds', 'ref_inds', 'num_tag_snps_varpop')]
    for i in range(len(args.target_pops) * 2):
        stat_hash['sfs_%d' % (i+1)] = 0
        pass
    for i in range(len(args.target_pops)):
        stat_hash['tmrca_varpop_com_ind_%d' % (i+1)] = 0
        pass
    if args.debug: supposed_stat_hash_keys = set(stat_hash)
    supposed_stat_hash_len = len(stat_hash)

    if args.first_line:
        print '\t'.join(['window_chr', 'window_start', 'window_end'] + sorted(stat_hash.keys()) + ['ms_outgroup_bl', 'ms_outgroup_bl_tot'])
        args.first_line = False
        pass

    snps = data['snps']
    stat_hash['num_snps'] = len(snps)
    window_chr = data['window_chr']
    window_start = data['window_start']
    window_end = data['window_end']
    window_label = '%s.%d.%d' % (window_chr, window_start, window_end)
    target_pop = args.target_pops

    stat_hash['target_inds'] = len(target_pop)
    stat_hash['ref_inds'] = num_genotypes - len(target_pop)
     
    ## get branch length for outgroup (only works with ms)
    if args.ms and args.ms_outgroup != None:
        ms_outgroup_bl_tot = args.ms_outgroup_branch_len[window_chr]
        hash = '%d-%d' % (window_start, window_end)
        if hash in args.ms_outgroup_branch_len_by_win[window_chr]:
            ms_outgroup_bl = args.ms_outgroup_branch_len_by_win[window_chr][hash]
        else:
            ms_outgroup_bl = 0
            pass
    else:
        ms_outgroup_bl_tot = 'NA'
        ms_outgroup_bl = 'NA'
        pass

    ## if there's no data in this window, just forget about it
    if len(snps) == 0:
        report = '\t'.join([str(k) for k in [window_chr, window_start, window_end] + [stat_hash[i] for i in sorted(stat_hash.keys())] + [ms_outgroup_bl, ms_outgroup_bl_tot]])
        print report
        return



    ## calculate various versions of S_star
    if args.debug: print_snps(snps)
    #[S_star, S_inds, S_snps] = calc_S_star(snps, ret_inds = True, prt=args.debug)
    [S_star_varpop, S_inds_varpop, S_snps_varpop] = calc_S_star(snps, varpop = target_pop, ret_inds = True, prt=args.debug)
    # [S_star_varpop_ti, S_inds_varpop_ti, S_snps_varpop_ti] = calc_S_star_i(snps, varpop = target_pop, ret_inds = True, prt=args.debug)
    # [S_star_ad0, S_inds_ad0, S_snps_ad0] = calc_S_star(snps, ret_inds = True, prt=args.debug, allowable_diffs = 0)
    [S_star_varpop_ad0, S_inds_varpop_ad0, S_snps_varpop_ad0] = calc_S_star(snps, varpop = target_pop, ret_inds = True, prt=args.debug, allowable_diffs = 0)
#     [S_star_varpop_ad1, S_inds_varpop_ad1, S_snps_varpop_ad1] = calc_S_star(snps, varpop = target_pop, ret_inds = True, prt=args.debug, allowable_diffs = 1)
#     [S_star_varpop_ad2, S_inds_varpop_ad2, S_snps_varpop_ad2] = calc_S_star(snps, varpop = target_pop, ret_inds = True, prt=args.debug, allowable_diffs = 2)
#     [S_star_varpop_ad3, S_inds_varpop_ad3, S_snps_varpop_ad3] = calc_S_star(snps, varpop = target_pop, ret_inds = True, prt=args.debug, allowable_diffs = 3)
#     [S_star_varpop_ad4, S_inds_varpop_ad4, S_snps_varpop_ad4] = calc_S_star(snps, varpop = target_pop, ret_inds = True, prt=args.debug, allowable_diffs = 4)
#     [S_star_varpop_ad5, S_inds_varpop_ad5, S_snps_varpop_ad5] = calc_S_star(snps, varpop = target_pop, ret_inds = True, prt=args.debug, allowable_diffs = 5)


    #stat_hash['S_star'] = S_star
    stat_hash['S_star_varpop'] = S_star_varpop
    stat_hash['S_star_varpop_ad0'] = S_star_varpop_ad0
    # stat_hash['S_star_varpop_ad1'] = S_star_varpop_ad1
    # stat_hash['S_star_varpop_ad2'] = S_star_varpop_ad2
    # stat_hash['S_star_varpop_ad3'] = S_star_varpop_ad3
    # stat_hash['S_star_varpop_ad4'] = S_star_varpop_ad4
    # stat_hash['S_star_varpop_ad5'] = S_star_varpop_ad5

    #if args.debug: print 'S :', S_star, S_snps, S_inds, get_common_in_sets(S_inds), max([len(s) for s in S_inds]) == len(get_common_in_sets(S_inds))
    if args.debug: print 'S varpop:', S_star_varpop, S_snps_varpop, S_inds_varpop, get_common_in_sets(S_inds_varpop)
    # if args.debug: print 'S varpop ti:', S_star_varpop_ti, S_snps_varpop_ti, S_inds_varpop_ti, get_common_in_sets(S_inds_varpop_ti), len(S_inds_varpop_ti[0]) == len(get_common_in_sets(S_inds_varpop_ti))
    if args.debug: print 'S varpop ad0:', S_star_varpop_ad0, S_snps_varpop_ad0, S_inds_varpop_ad0, get_common_in_sets(S_inds_varpop_ad0)

    ## get involved individuals
    vp_com_inds = get_common_in_sets(S_inds_varpop)
    vp_maj_inds = get_majority_in_sets(S_inds_varpop)
    ad0_com_inds = get_common_in_sets(S_inds_varpop_ad0)
    #ad0_maj_inds = get_majority_in_sets(S_inds_varpop_ad0)

    ## get snps for involved individuals
    #print vp_com_inds
    #print S_inds_varpop, len(S_inds_varpop)
    #print S_snps_varpop, len(S_snps_varpop)
    vp_com_snps = [s for i,s in enumerate(S_snps_varpop) if len(set(S_inds_varpop[i]).intersection(set(vp_com_inds))) != 0]
    #print vp_com_snps
    #print [len(set(S_inds_varpop[i]).intersection(set(vp_com_inds))) for i,s in enumerate(S_snps_varpop)]
    #print snps[vp_com_snps[-1]]['pos'] - snps[vp_com_snps[0]]['pos']
    stat_hash['haplen_varpop'] = snps[vp_com_snps[-1]]['pos'] - snps[vp_com_snps[0]]['pos'] if len(vp_com_snps) > 1 else 0


    ## get basic tmrcas
    stat_hash['tmrca_all'] = calc_tmrca(data, prt=args.debug)[1]
    stat_hash['tmrca_target_pop'] = calc_tmrca(data, subset = target_pop, prt=args.debug)[1]
    stat_hash['tmrca_varpop'] = calc_tmrca(data, subset = vp_com_inds, prt=args.debug)[1]
    tmp_tmrcas = [calc_tmrca(data, subset = [ind], prt=args.debug)[1] for ind in vp_com_inds]
    for i in range(len(tmp_tmrcas)):
        stat_hash['tmrca_varpop_com_ind_%d' % (i+1)] = tmp_tmrcas[i]
        pass
    stat_hash['tmrca_varpop_maj'] = calc_tmrca(data, subset = vp_maj_inds, prt=args.debug)[1]
    #stat_hash['tmrca_varpop_bar'] = calc_tmrca(data, subset = remove_common_in_sets(range(num_genotypes), S_inds_varpop), prt=args.debug)[1]
    stat_hash['varpop_common'] = ' '.join([str(s) for s in vp_com_inds])
    stat_hash['varpop_maj'] = ' '.join([str(s) for s in vp_maj_inds])
    stat_hash['count_varpop_common'] = len(vp_com_inds)
    stat_hash['count_varpop_maj'] = len(vp_maj_inds)

    stat_hash['tmrca_varpop_ad0'] = calc_tmrca(data, subset = ad0_com_inds, prt=args.debug)[1]
    #stat_hash['tmrca_varpop_ad0_maj'] = calc_tmrca(data, subset = ad0_maj_inds, prt=args.debug)[1]
    #stat_hash['tmrca_varpop_ad0_bar'] = calc_tmrca(data, subset = remove_common_in_sets(range(num_genotypes), S_inds_varpop_ad0), prt=args.debug)[1]
    stat_hash['varpop_ad0_common'] = ' '.join([str(s) for s in ad0_com_inds])
    #stat_hash['varpop_ad0_maj'] = ' '.join([str(s) for s in ad0_maj_inds])
    stat_hash['count_varpop_ad0_common'] = len(ad0_com_inds)
    #stat_hash['count_varpop_ad0_maj'] = len(ad0_maj_inds)
    stat_hash['tmrca_non_target_pop'] = calc_tmrca(data, subset = args.non_target_pop, prt=args.debug)[1]

    sfs_varpop = [sum([snps[s]['snp_details'][j][2] for j in target_pop]) for s in S_snps_varpop]
    if args.debug: print 'varpop_common inds:', S_inds_varpop, S_star_varpop
    if args.debug: print 'varpop snps:', S_snps_varpop, len(snps)
    if args.debug: print 'sfs:', sfs_varpop
    if args.debug: print 'tmrca, tmrca avg:', stat_hash['tmrca_varpop'], [stat_hash['tmrca_varpop_com_ind_%d' % (i+1)] for i in range(len(target_pop))]
    sfs_varpop = Counter(sfs_varpop)
    if args.debug: print 'sfs dict:', sfs_varpop
    if args.debug: print 'max sfs dict key:', max(sfs_varpop.keys()) if len(sfs_varpop.keys()) > 0 else 0
    for i in range(len(args.target_pops * 2)):
        stat_hash['sfs_%d' % (i+1)] = sfs_varpop[i+1]
        pass
    stat_hash['num_tag_snps_varpop'] = len(S_snps_varpop)

    #### Report

    if len(stat_hash) != supposed_stat_hash_len:
        print "we have incorrect stat hash!"
        print sorted(stat_hash)
        print [stat_hash[k] for k in sorted(stat_hash)]
        print len(stat_hash), supposed_stat_hash_len
        if args.debug: print set(stat_hash).difference(supposed_stat_hash_keys)
        sys.exit(-1)
        return
    

    report = '\t'.join([str(k) for k in [window_chr, window_start, window_end] + [stat_hash[i] for i in sorted(stat_hash.keys())] + [ms_outgroup_bl, ms_outgroup_bl_tot]])
    if args.report_snps:
        #print snps
        present_in_nontarget = [sum([s['snp_details'][j][2] > 0 for j in args.non_target_pop]) > 0 for s in snps]
        report = '\t'.join([report, 'snps', str(len(snps))] + ['%d %s %d %d' % (snps[s]['pos'], snps[s]['var'], s in S_snps_varpop, present_in_nontarget[s]) for s in range(len(snps))])
        pass
    print report
    sys.stdout.flush()
    pass
    
#     else:
#         print "oh shit invalid method!!"
#         print args
#         print args.method_number
#         print args.method_args
#         sys.exit(-1)
#         pass

    return



    

def calc_tmrca(data, subset = None, prt = True):
    debug_tmrca = False

    snps = data['snps']
    window_chr = data['window_chr']
    window_start = data['window_start']
    window_end = data['window_end']

    ## x is the number of derived variants present in the sample, for this region (summed over all individuals)
    if subset == None:
        x = sum([snp['x'] for snp in snps])
    else:
        x = sum([calc_x(snp, subset) for snp in snps])
        pass
    # print subset, x, [calc_x(snp, subset) for snp in snps]
    # x_dropone = [sum([snp['subx'][i] for snp in snps]) for i in range(num_genotypes)]
    
    [num_chrs, max_missing_chrs, sex_chr] = calc_chr_missing_etc(window_chr, window_start + 1, subset)
    [num_chrs_e, max_missing_chrs_e, sex_chr_e] = calc_chr_missing_etc(window_chr, window_end, subset)

    if num_chrs != num_chrs_e or max_missing_chrs != max_missing_chrs_e:
        print "ERROR IN NUMBER OF HAPLOTYPE LOGIC"
        print num_chrs, num_chrs_e, max_missing_chrs, max_missing_chrs_e
        print window_chr, window_start, window_end
        sys.exit(-1)
        pass

    ##### calculate T_hat

    putative_window_size = window_end - window_start
    # nocalls has to be 100% right now, so I just exclude them
    # total_nocalls = calc_nocalls( window_chr, window_start, window_end, num_chrs, max_missing_chrs )

    ## calculate D and mutation rate (if we're using chimp data)
    if args.chimp_snps != None:
        chimp_snps = data['chimp_snps']
        chimp_mapped = chimp_mapped_regions.amount_in_region(window_chr, window_start, window_end)
        if chimp_mapped == 0:
            D = 0
        else:
            D = len(chimp_snps) / chimp_mapped
            pass
        mut_rate = D / time_to_primate_divergence / 2
    else:
        mut_rate = 2.5e-8 / 25 # mutation rate per year
        D = mut_rate * time_to_primate_divergence * 2
        chimp_mapped = window_end - window_start
        chimp_snps = [0] * int(D*chimp_mapped)
        if args.debug: print 'ms chimp stats:', D, chimp_mapped, mut_rate, D / time_to_primate_divergence / 2, len(chimp_snps)
        pass

    ## calculate T_hat
    N_adj = data['region_length']
    if N_adj == 0 or num_chrs == 0:
        T_hat = 0
    else:
        T_hat = x / num_chrs / N_adj
        pass
    
    if D != 0:
        tmrca = 2 * T_hat * time_to_primate_divergence / D
    else:
        tmrca = 0
        pass

    if N_adj < putative_window_size * args.required_window_fraction or chimp_mapped < putative_window_size * args.required_window_fraction:
        tmrca_filtered = 0
        mut_rate_filtered = 0
        T_hat_filtered = 0
    else:
        tmrca_filtered = tmrca
        T_hat_filtered = T_hat
        mut_rate_filtered = mut_rate
        pass

    if debug_tmrca: print
    if debug_tmrca: print "x = %f,    num_chrs = %d,    N_put = %f,    avg_nocalls = %f,    N_after_nocalls = %f,    T_hat = %f,     d = %f,     N_d = %f,     D = %f" % (x, num_chrs, putative_window_size, total_nocalls[0] / num_chrs, N_adj[0], T_hat[0], len(chimp_snps), chimp_mapped, D)

    stat = (tmrca, tmrca_filtered, num_chrs, len(snps), float(N_adj), chimp_mapped, mut_rate, mut_rate_filtered, T_hat, T_hat_filtered,
            "2 * %d * %d / ((%d / %f) * %d * %f)" % (x, time_to_primate_divergence, len(chimp_snps), chimp_mapped, num_chrs, N_adj)
            #2 * x * time_to_primate_divergence / ((len(chimp_snps) / chimp_mapped) * num_chrs * N_adj[i]),
            )
    
    if prt: print window_chr, window_start, window_end, ' '.join([str(i) for i in stat])
    return stat

def snp_str(snps):
    ret = "\n"
    for s in snps:
        #ret += str(s[:5])
        ret += str(s)
        ret += '\n'
        pass
    return ret


def print_dist_matrix(matrix):
    print num_genotypes+1
    for i, ind in enumerate(individuals_labels + ['chimp']):
        s = "%-10s " % ind
        s += ' '.join([str(n) for n in matrix[i]])
        print s
        sys.stdout.flush()
        pass
    return

def print_ind_dists(matrix):
    for i, ind in enumerate(individuals_labels):
        s = "%-10s " % ind
        s += str(sum(matrix[i][:-1]))
        print "filter:", s
        sys.stdout.flush()
        pass
    return



#next_snp = None
#if if 'chimp_snps' in args: next_snp_chimp = read_chimp()


####
## metric functions

# distance matrix
one_dist_matrix = None
def calc_dist_matrix(data):

    #print "THIS DOES NOT ACCOUNT FOR CHIMP SNPS (LOOK TO MK CODE TO SEE ONE WAY TO DO IT).  IT ALSO MIGHT NOT PROPERLY ACCOUNT FOR CHIMP MAPPED REGIONS, ETC."
    #sys.exit(-1)

    null_matrix = array([[0] * (num_genotypes+1)] * (num_genotypes+1))
    # if data['is_window'] and one_dist_matrix == None: one_dist_matrix = array(null_matrix)
    chimp_dist_matrix = array(null_matrix)
    counted_snps = 0

    if args.resample:
        fileh = tables.openFile("tmp.dat.%d.h5" % random.randint(1,1000000), mode = "w")
        all_dist_matrices = fileh.createTable("/", 'arrays', {'a':tables.Int32Col(shape=(num_genotypes+1,num_genotypes+1))})
        pass

    chr = ''
    
    for next_snp in data['snps']:

        ancestral_allele = next_snp['anc']
        if ancestral_allele.islower():
            chimp_snp = (2, 0, 0)
        elif ancestral_allele == next_snp['ref']:
            chimp_snp = (0, 2, 0)
        elif ancestral_allele == next_snp['var']:
            chimp_snp = (0, 0, 2)
        else:
            chimp_snp = (2, 0, 0)
            pass
        snp_dist_matrix = array([[max(0, abs(s1[1] - s2[1]) - max(s1[0], s2[0])) for s1 in next_snp['snp_details'] + [chimp_snp]] for s2 in next_snp['snp_details'] + [chimp_snp]])
        chimp_dist_matrix += snp_dist_matrix
        if args.resample:
            all_dist_matrices.row['a'] = snp_dist_matrix
            all_dist_matrices.row.append()
            pass
        counted_snps += 1
        #print "snp:", counted_snps % 1000
        #if counted_snps % 1000 == 0: print 'debug', "snp:", counted_snps
        pass
    if args.resample: all_dist_matrices.flush()
    
#     if data['is_window']:
#         one_dist_matrix += dist_matrix
#         pass
    
    print_dist_matrix(chimp_dist_matrix)

    if args.resample:
        print 'debug', "number of matrices:", counted_snps
        for i in xrange(1000):
            #samples = sample_wr(xrange(counted_snps), 100000)
            samples = sample_wr(xrange(counted_snps), counted_snps)
            #print 'debug', i, samples, all_dist_matrices
            #resample_matrix = sum([t['a'] for t in all_dist_matrices[samples]])
            resample_matrix = array(null_matrix)
            for s in samples:
                resample_matrix += all_dist_matrices[s]['a']
                pass
            print_dist_matrix(resample_matrix)
            pass
        pass
    

#     print_dist_matrix(one_dist_matrix)
#             sys.stderr.write( "starting " + next_snp['chr'] + 'from' + current_chr + '\n' )
#             sys.stderr.flush()
#             current_chr = next_snp['chr']
#             pass

#         pass
#     print "filter:", "whole genome!"
#     print_dist_matrix(one_dist_matrix + chimp_dist_matrix)
#     if calc_windowed_dist_matrix:
#         snp_pos = [snp['pos'] for snp in list_of_snps]
#         valid_chimp_snp_pos = [snp['pos'] for snp in list_of_snps_chimp if snp['pos'] not in snp_pos]
#         valid_chimp_snp_len = len(valid_chimp_snp_pos)
#         #print snp_pos
#         #print valid_chimp_snp_pos
#         #print valid_chimp_snp_len
#         window_matrix = sum([null_matrix] + [snp['dist-matrix'] for snp in list_of_snps])
#         chimp_matrix = array(null_matrix)
#         chimp_matrix[-1] = valid_chimp_snp_len
#         for i in range(num_genotypes+1):
#             chimp_matrix[i][num_genotypes] = chimp_matrix[num_genotypes][i]
#             pass
#         window_matrix += chimp_matrix
#         # keep track of total chimp dists - remember that each window overlaps with, on average, 10 other windows
#         ##  THIS IS NOT EXACT, BUT PROBABLY GOOD ENOUGH!!
#         chimp_dist_matrix += (chimp_matrix / 10)
#         #print "chimp matrices"
#         #print_dist_matrix(chimp_matrix)
#         #print_dist_matrix(chimp_dist_matrix)
#         #print window_matrix
#         #print [[null_matrix] + snp['dist-matrix'] for snp in list_of_snps]
#         [num_chrs, max_missing_chrs, sex_chr] = calc_chr_missing_etc(window_chr, window_start + 1)
#         putative_window_size = window_end - window_start
#         total_nocalls = calc_nocalls( window_chr, window_start, window_end, num_chrs, max_missing_chrs )
#         #print total_nocalls
#         chimp_mapped = chimp_mapped_regions.amount_in_region(window_chr, window_start, window_end)
#         if len(list_of_snps) >= 10 and total_nocalls[0] / num_chrs < putative_window_size / 2 and (putative_window_size - chimp_mapped) < putative_window_size / 2:
#             print "filter-window1:", window_chr, window_start, window_end
#             print "filter-window2:", len(list_of_snps), total_nocalls[0] / num_chrs, putative_window_size - chimp_mapped
#             print_dist_matrix(window_matrix)
#             print_ind_dists(window_matrix)
#             pass
#         pass
    return

def add_1(i, lst):
    if i not in lst: lst.append(i)
    return

def add_2(i, lst):
    if i not in lst: lst.append(i)
    return


def calc_summary_stats_windowed(data):
    
    snp_list = data['snps']
    if len(snp_list) < 30 or data['region_length'] < 20000: return

    if 'start_time' not in data['persist']:
        data['persist']['start_time'] = time.time()
        pass

    if 'dc_freq_results' not in data['persist']:
        data['persist']['dc_freq_results'] = init_freqs([len(pop)*2 for pop in args.multiple_populations], prefix = 'dcf')
        pass
    dc_freq_results = data['persist']['dc_freq_results']

    if 'all_freq_results' not in data['persist']:
        data['persist']['all_freq_results'] = init_freqs([len(pop)*2 for pop in args.multiple_populations], prefix = 'f')
    all_freq_results = data['persist']['all_freq_results']

    if 'sstar_freq_results' not in data['persist']:
        data['persist']['sstar_freq_results'] = init_freqs([len(pop)*2 for pop in args.multiple_populations], prefix = 'sstar')
    sstar_freq_results = data['persist']['sstar_freq_results']

    if 'sstar_freq_results2' not in data['persist']:
        data['persist']['sstar_freq_results2'] = defaultdict(lambda : init_freqs([len(pop)*2 for pop in args.multiple_populations], prefix = 'sstar'))
    sstar_freq_results2 = data['persist']['sstar_freq_results2']

    if 'fst_results' not in data['persist']:
        store_fst = True
        data['persist']['fst_results'] = None
    else:
        store_fst = False
        pass
    fst_results = data['persist']['fst_results']

    if 'pi1_results' not in data['persist']:
        store_pi = True
        data['persist']['pi1_results'] = None
        data['persist']['pi2_results'] = None
    else:
        store_pi = False
        pass
    pi1_results = data['persist']['pi1_results']
    pi2_results = data['persist']['pi2_results']

    if 's_star_snps' not in data['persist']:
        data['persist']['s_star_snps'] = []
        data['persist']['s_star_snps2'] = []
        data['persist']['s_star_snps_c'] = Counter()
        pass

    if 's_stars' not in data['persist']:
        data['persist']['s_stars'] = []
        pass

    if 'reglens' not in data['persist']:
        data['persist']['reglens'] = []
        pass
    data['persist']['reglens'].append(data['region_length'])

    for next_snp in snp_list:

        fst_results = calc_fst_freal(next_snp, fst_results, drop_singletons = True)
        pi1_results = calc_pi_freal(next_snp['mpops'][0], pi1_results)
        pi2_results = calc_pi_freal(next_snp['mpops'][1], pi2_results)

        if not next_snp['mrca_known']: continue

        ref_as_minor = consider_ancestry and next_snp['mrca'] == 1
        if (ref_as_minor and next_snp['ref'] in next_snp['archaic_snps']) or (not ref_as_minor and next_snp['var'] in next_snp['archaic_snps']):
            dc_freq_results = calc_freqs_freal(next_snp['mpops'], ref_as_minor = ref_as_minor, results = dc_freq_results, prefix = 'dcf')
            pass
        all_freq_results = calc_freqs_freal(next_snp['mpops'], ref_as_minor = ref_as_minor, results = all_freq_results, prefix = 'f')

        pass

    ## fix this
    varpop_indices = [i for i,s in enumerate(snp_list) if get_nonvarcount_sum(s, individuals, args.target_pops) == 0]
    # varpop_indices = [i for i,s in enumerate(snp_list) if get_nonvarcount_sum(s, individuals, args.target_pops) == 0 and s['mrca'] == 0]

    #cache = None
    if 'cache' not in data['persist']:
        data['persist']['cache'] = {}
        pass
    ## hack to keep the cache from growing out of control
    if len(data['persist']['cache'].keys()) > 5000: data['persist']['cache'] = {}
    cache = data['persist']['cache']
    #print 'cache', len(cache.keys())
#    for i in range(len(snp_list)):
    for i in varpop_indices:
        i_pos = snp_list[i]['pos']
        if i_pos not in cache: cache[i_pos] = {}
        for j in range(len(snp_list)):
            j_pos = snp_list[j]['pos']
            if j_pos not in cache[i_pos]: cache[i_pos][j_pos] = None
            pass
        pass

    [S_star_varpop, S_inds_varpop, S_snps_varpop] = calc_S_star(snp_list, snp_indices = varpop_indices, ret_inds = True, prt=args.debug, cache = cache)
    data['persist']['s_stars'].append(S_star_varpop)

    for i in S_snps_varpop:
        if 's_star' not in snp_list[i]: snp_list[i]['s_star'] = 0
        snp_list[i]['s_star'] = max(S_star_varpop, snp_list[i]['s_star'])
        # add_1( snp_list[i], data['persist']['s_star_snps'] )
        add_2( snp_list[i], data['persist']['s_star_snps2'] )
        # print 's_star snps:', len(data['persist']['s_star_snps'])
        pass
    # data['persist']['s_star_snps'].sort(key = lambda s : s['s_star'])
    # print len(data['persist']['s_star_snps']), [s['s_star'] for s in data['persist']['s_star_snps']]
    # print S_star_varpop
    # print S_inds_varpop
    # print S_snps_varpop
    #print 'adding for window', data['window_chr'], data['window_start'], S_star_varpop
    #print 'with s_star_snps', data['persist']['s_star_snps2']
    for next_snp in [s for s in data['persist']['s_star_snps2'] if s['chr'] != data['window_chr'] or s['pos'] < data['window_start']]:
        if not next_snp['mrca_known'] or next_snp['mrca'] == 1: continue
        sstar_freq_results = calc_freqs_freal(next_snp['mpops'], ref_as_minor = False, results = sstar_freq_results, prefix = 'sstar')
        sstar_freq_results2[S_star_varpop] = calc_freqs_freal(next_snp['mpops'], ref_as_minor = False, results = sstar_freq_results2[S_star_varpop], prefix = 'sstar')
        # print 'debug sstar snps', S_star_varpop, sstar_freq_results2
        pass
    data['persist']['s_star_snps_c'].update([s['s_star'] for s in data['persist']['s_star_snps2'] if s['chr'] != data['window_chr'] or s['pos'] < data['window_start']])
    data['persist']['s_star_snps2'][:] = [s for s in data['persist']['s_star_snps2'] if s['chr'] == data['window_chr'] and s['pos'] >= data['window_start']]

    if store_fst:
        data['persist']['fst_results'] = fst_results
        pass

    if store_pi:
        data['persist']['pi1_results'] = pi1_results
        data['persist']['pi2_results'] = pi2_results
        pi1_results['reglen'] = 0
        pi2_results['reglen'] = 0
        pass
    pi1_results['reglen'] += data['region_length']
    pi2_results['reglen'] += data['region_length']

    if 'num_snps' not in data['persist']:
        data['persist']['num_snps'] = []
        #print '\t'.join(['stats', 'snps', 'reglen', 'sstar', 'time'])
        pass
    data['persist']['num_snps'].append(len(snp_list))

    elapsed_time = time.time() - data['persist']['start_time']
    data['persist']['start_time'] = time.time()

    #print '\t'.join([str(s) for s in ['stats', len(snp_list), data['region_length'], S_star_varpop, elapsed_time]])
    
    return


def report_summary_stats(data):
    # total_snps = sum(dc_freq_results['sfs'].values())
    # print '\t'.join(sorted(dc_freq_results['sfs'].keys()) + ['fst'])
    # print '\t'.join([str(dc_freq_results['sfs'][f] / total_snps) for f in sorted(dc_freq_results['sfs'].keys())] + ['%f' % fst_results['fst']])
    # total_snps = sum(all_freq_results['sfs'].values())
    # print '\t'.join(sorted(all_freq_results['sfs'].keys()) + ['fst'])
    # print '\t'.join([str(all_freq_results['sfs'][f] / total_snps) for f in sorted(all_freq_results['sfs'].keys())] + ['%f' % fst_results['fst']])

    header = []
    results = []

    dc_freq_results = data['persist']['dc_freq_results']
    all_freq_results = data['persist']['all_freq_results']
    sstar_freq_results = data['persist']['sstar_freq_results']
    sstar_freq_results2 = data['persist']['sstar_freq_results2']
    fst_results = data['persist']['fst_results']
    pi1_results = data['persist']['pi1_results']
    pi2_results = data['persist']['pi2_results']
    S_star_varpop = data['persist']['s_stars'][-1]

    #data['persist']['s_star_snps'].sort(key = lambda s : s['s_star'])

    ## have to capture the s_star values from the last window
    ## because we count the maximum sstar for a snp, so we can't capture it until the window is past it
    #print 'finally adding for window', data['window_chr'], data['window_start']
    #print 'with s_star_snps', data['persist']['s_star_snps2']
    data['persist']['s_star_snps_c'].update([s['s_star'] for s in data['persist']['s_star_snps2']])
    for next_snp in [s for s in data['persist']['s_star_snps2']]:
        if not next_snp['mrca_known'] or next_snp['mrca'] == 1: continue
        sstar_freq_results = calc_freqs_freal(next_snp['mpops'], ref_as_minor = False, results = sstar_freq_results, prefix = 'sstar')
        sstar_freq_results2[S_star_varpop] = calc_freqs_freal(next_snp['mpops'], ref_as_minor = False, results = sstar_freq_results2[S_star_varpop], prefix = 'sstar')
        pass

    header.append('fst')
    results.append(fst_results['fst'])
    header.append('fst_h')
    results.append(fst_results['fst_h'])
    #print 'fst', fst_results['fst']
    header.append('pi1')
    results.append(pi1_results['pi_sum'] / pi1_results['reglen'])
    #print 'pi1', pi1_results['pi_sum'] / pi1_results['reglen']
    header.append('pi2')
    results.append(pi2_results['pi_sum'] / pi2_results['reglen'])
    #print 'pi2', pi2_results['pi_sum'] / pi2_results['reglen']


    # s_star_results = init_freqs([len(pop)*2 for pop in args.multiple_populations], prefix = 'sstar')
    # for next_snp in data['persist']['s_star_snps']:
    #     ref_as_minor = consider_ancestry and next_snp['mrca'] == 1
    #     s_star_results = calc_freqs_freal(next_snp['mpops'], ref_as_minor = ref_as_minor, results = s_star_results, prefix = 'sstar')
    #     pass

    #s_star_results_subset = init_freqs([len(pop)*2 for pop in args.multiple_populations], prefix = 'sstar_s')
    #s_star_subset = data['persist']['s_star_snps'][int(len(data['persist']['s_star_snps']) * .9):]
    #header.append('s_star_snps_count')
    #results.append(len(data['persist']['s_star_snps']))
    #print len(s_star_subset)

    # for next_snp in s_star_subset:
    #     ref_as_minor = consider_ancestry and next_snp['mrca'] == 1
    #     s_star_results_subset = calc_freqs_freal(next_snp['mpops'], ref_as_minor = ref_as_minor, results = s_star_results_subset, prefix = 'sstar_s')
    #     pass

    total_snps = sum(dc_freq_results['sfs'].values())
    header.extend(sorted(dc_freq_results['sfs'].keys()))
    results.extend([dc_freq_results['sfs'][f] / total_snps for f in sorted(dc_freq_results['sfs'].keys())])

    total_snps = sum(all_freq_results['sfs'].values())
    header.extend(sorted(all_freq_results['sfs'].keys()))
    results.extend([all_freq_results['sfs'][f] / total_snps for f in sorted(all_freq_results['sfs'].keys())])

    total_snps = sum(sstar_freq_results['sfs'].values())
    header.extend(sorted(sstar_freq_results['sfs'].keys()))
    results.extend([sstar_freq_results['sfs'][f] / total_snps for f in sorted(sstar_freq_results['sfs'].keys())])

    # print 'sstar debug', total_snps
    # print sum([sstar_freq_results['sfs'][f] / total_snps for f in sorted(sstar_freq_results['sfs'].keys())])
    # print [sstar_freq_results['sfs'][f] / total_snps for f in sorted(sstar_freq_results['sfs'].keys())]

    # print sum(data['persist']['sstar_freq_results']['sfs'].values()), [data['persist']['sstar_freq_results']['sfs'][f] / total_snps for f in sorted(data['persist']['sstar_freq_results']['sfs'].keys())]
    # print ['%s %d' % (f, s_star_results['sfs'][f]) for f in sorted(s_star_results['sfs'].keys())]
    # print ['%s %d' % (f, data['persist']['sstar_freq_results']['sfs'][f]) for f in sorted(data['persist']['sstar_freq_results']['sfs'].keys())]
    # print s_star_results
    # print data['persist']['sstar_freq_results']

    #total_snps = sum(s_star_results_subset['sfs'].values())
    #header.extend(sorted(s_star_results_subset['sfs'].keys()))
    #results.extend([s_star_results_subset['sfs'][f] / total_snps for f in sorted(s_star_results_subset['sfs'].keys())])

    quantiles = [0, .1, .2, .25, .3, .4, .5, .6, .7, .75, .8, .9, 1]
    
    header.extend(['reglen_quant_%.2f' % q for q in quantiles])
    data['persist']['reglens'].sort()
    results.extend([data['persist']['reglens'][int((len(data['persist']['reglens'])-1) * q)] for q in quantiles])

    sstar_snps = [item for k in sorted(data['persist']['s_star_snps_c'].keys()) for item in [k for i in xrange(data['persist']['s_star_snps_c'][k])]]
    header.append('s_star_snps_count')
    results.append(len(sstar_snps))
    header.extend(['s_star_quant_%.2f' % q for q in quantiles])
    results.extend([sstar_snps[int((len(sstar_snps)-1) * q)] for q in quantiles])

    # print data['persist']['s_star_snps_c']
    # print len(data['persist']['s_star_snps']), [s['s_star'] for s in data['persist']['s_star_snps']]
    # print len(sstar_snps), sstar_snps
    # # print [data['persist']['s_star_snps'][int((len(data['persist']['s_star_snps'])-1) * q)]['s_star'] for q in quantiles]
    # print [sstar_snps[int((len(sstar_snps)-1) * q)] for q in quantiles]

    header.extend(['snps_quant_%.2f' % q for q in quantiles])
    data['persist']['num_snps'].sort()
    results.extend([data['persist']['num_snps'][int((len(data['persist']['num_snps'])-1) * q)] for q in quantiles])

    print '\t'.join(header)
    print '\t'.join([str(s) for s in results])

    return

def calc_summary_stats(data):
    
    snp_list = data['snps']
    dc_freq_results = init_freqs([len(pop)*2 for pop in args.multiple_populations], prefix = 'dcf')
    all_freq_results = init_freqs([len(pop)*2 for pop in args.multiple_populations], prefix = 'f')
    fst_results = None

    for next_snp in snp_list:

        ref_as_minor = consider_ancestry and next_snp['mrca'] == 1
        if (ref_as_minor and next_snp['ref'] in next_snp['archaic_snps']) or (not ref_as_minor and next_snp['var'] in next_snp['archaic_snps']):
            dc_freq_results = calc_freqs_freal(next_snp['mpops'], ref_as_minor = ref_as_minor, results = dc_freq_results, prefix = 'dcf')
            pass
        
        ## for multiple sfs, keep more snps, and check before adding a snp to some sfs calcs (i.e., does it match the comparison genome?)
        ## if doing multiple sfs, then add the ability to change the prefix from 'f_' to something else
        all_freq_results = calc_freqs_freal(next_snp['mpops'], ref_as_minor = ref_as_minor, results = all_freq_results, prefix = 'f')
        fst_results = calc_fst_freal(next_snp, fst_results)

        pass

    # total_snps = sum(dc_freq_results['sfs'].values())
    # print '\t'.join(sorted(dc_freq_results['sfs'].keys()) + ['fst'])
    # print '\t'.join([str(dc_freq_results['sfs'][f] / total_snps) for f in sorted(dc_freq_results['sfs'].keys())] + ['%f' % fst_results['fst']])
    # total_snps = sum(all_freq_results['sfs'].values())
    # print '\t'.join(sorted(all_freq_results['sfs'].keys()) + ['fst'])
    # print '\t'.join([str(all_freq_results['sfs'][f] / total_snps) for f in sorted(all_freq_results['sfs'].keys())] + ['%f' % fst_results['fst']])

    total_snps = sum(dc_freq_results['sfs'].values())
    print '\t'.join(sorted(dc_freq_results['sfs'].keys()))
    print '\t'.join([str(dc_freq_results['sfs'][f] / total_snps) for f in sorted(dc_freq_results['sfs'].keys())])
    total_snps = sum(all_freq_results['sfs'].values())
    print '\t'.join([str(all_freq_results['sfs'][f] / total_snps) for f in sorted(all_freq_results['sfs'].keys())])

    return




def get_freq_cat(fqs, prefix):
    return '%s_%.4f' % (prefix, fqs[0]) if len(fqs) == 1 else '%s_%.4f_%.4f' % (prefix, fqs[0], fqs[1])

def calc_freqs_freal(snp_list, ref_as_minor = False, results = None, get_freq_cat = get_freq_cat, prefix = 'f'):
    freqs = [snp['refs'] / snp['haps'] if ref_as_minor else snp['vars'] / snp['haps'] for snp in snp_list]

    if min([snp['haps'] for snp in snp_list]) != 18: print 'haps', snp_list

    if results == None:
        results = {'sfs' : Counter()}
        pass

    results['freqs'] = freqs
    results['sfs'].update([get_freq_cat(freqs, prefix)])
    return results

    
def init_freqs(haps, results = None, get_freq_cat = get_freq_cat, prefix = 'f'):

    if results == None:
        results = {'sfs' : Counter()}
        pass

    if len(haps) == 1:
        for f in range(haps[0]+1):
            results['sfs'][get_freq_cat([f/haps[0]], prefix)] = 0
            pass
        pass
    else:
        for f1 in range(haps[0]+1):
            for f2 in range(haps[1]+1):
                if (f1 == 0 and f2 == 0) or (f1 == haps[0] and f2 == haps[1]): continue
                results['sfs'][get_freq_cat([f1/haps[0],f2/haps[1]], prefix)] = 0
                pass
            pass
        pass
    return results
    

def calc_freqs(data):

    snp_list = data['snps']
    freq_results = init_freqs([len(pop)*2 for pop in args.multiple_populations])
    #print [len(pop) for pop in args.multiple_populations]
    
    for next_snp in snp_list:

        ref_as_minor = consider_ancestry and next_snp['mrca'] == 1
        if args.debug_read_snps: print 'ref as minor', consider_ancestry, ref_as_minor, next_snp['mrca'], next_snp

        if args.multiple_populations == None:
            freq_results = calc_freqs_freal([next_snp], ref_as_minor = ref_as_minor, results = freq_results)
            if not args.summary_only:
                print 'implement this'
                #print '\t'.join([str(s) for s in [next_snp['chr'], next_snp['pos']-1, next_snp['pos'], next_snp['ref'], next_snp['var'], next_snp['mrca'], f, next_snp['archaic_snps'], next_snp['anc']]])
                #print "%s\t%d\t%d\t%s\t%f" % (next_snp['chr'], next_snp['pos']-1, next_snp['pos'], next_snp['var'], f)
                pass
            pass
        else:
            freq_results = calc_freqs_freal(next_snp['mpops'], ref_as_minor = ref_as_minor, results = freq_results)
            if not args.summary_only:
                print 'implement this'
                # print '\t'.join([str(s) for s in [next_snp['chr'], next_snp['pos']-1, next_snp['pos'], next_snp['ref'], next_snp['var'], next_snp['mrca'], freqs, next_snp['archaic_snps'], next_snp['anc']]])
                #print "\t".join((str(s) for s in [next_snp['chr'], next_snp['pos']-1, next_snp['pos'], next_snp['var']] + freqs + haps))
                pass
            pass
        pass

    if args.summary_only:
        total_snps = sum(freq_results['sfs'].values())
        print '\t'.join(sorted(freq_results['sfs'].keys()))
        print '\t'.join([str(freq_results['sfs'][f] / total_snps) for f in sorted(freq_results['sfs'].keys())])
        pass
    
    return


def calc_sites(data):

    snp_list = data['snps']
    
    for next_snp in snp_list:
        print '%s\t%d\t%d\t%.4f\t%s' % (next_snp['chr'], next_snp['pos']-1, next_snp['pos'], next_snp['vars'] / next_snp['haps'], '\t'.join([str(s) for s in next_snp['snp_details']]))
        pass

    return


def calc_msfile(data):

    if data['is_window'] and args.output_file.name != '<stdout>':
        new_file = args.output_file.name + '.' + data['region_tag'].replace('\t', '.')
        sys.stdout = open(new_file, 'w')
        print >> sys.stderr, 'changing output file to', new_file
    else:
        #print >> sys.stderr, 'not a window'
        pass

    # print data['region_tag']

    snp_list = data['snps']

    segsites = 0
    positions = []
    lines = [['',''] for i in range(num_genotypes)]

    for next_snp in snp_list:
        segsites += 1
        positions += [next_snp['pos']]

        for ind in range(num_genotypes):

            ## var is derived
            if next_snp['mrca'] == 0:
                ## one var?
                lines[ind][0] += str(int(next_snp['snp_details'][ind][2] > 0))
                ## two vars?
                lines[ind][1] += str(int(next_snp['snp_details'][ind][2] > 1))
            else:
                ## ref is derived
                ## one ref?
                lines[ind][0] += str(int(next_snp['snp_details'][ind][1] > 0))
                ## two refs?
                lines[ind][1] += str(int(next_snp['snp_details'][ind][1] > 1))
                pass
            
            if next_snp['snp_details'][i][0] != 0:
                print 'oh shit you need to set -fm 1!'
                print next_snp, next_snp['snp_details'][i]
                sys.exit(-1)
                pass
            pass
        pass

    if args.first_line:
        print 'ms %d 1000000000 -t 10' % (num_genotypes * 2)
        print '12789 46813 42616'
        args.first_line = False
        pass
    
    print
    print '//'
    #print data['region_length']
    print 'segsites: %d' % segsites
    if segsites > 0:
        print 'positions: ' + ' '.join(['%.15f' %  (p / args.baselookup.chr_lens[data['window_chr']]) for p in positions])
        for pop in args.keep_pops:
            # print pop
            for ind in range(num_genotypes):
                if individuals_pops[ind] == pop:
                    print lines[ind][0]
                    print lines[ind][1]
                    pass
                pass
            pass
        pass
    else:
        print
        pass
    return

#lookup = args.baselookup('/net/akey/vol1/home/bvernot/tishkoff/primate_sequences')

setattr(args, 'catlines', ['%-10s' % s for s in ('chimp', 'denisovan', 'archaic', 'human')])
def calc_visual_fasta(data):

    #     while True:
    #         print 'pos?'
    #         i = int(sys.stdin.readline().strip())
    #         den_base = args.bsg_files[1].get_base_one_based(data['window_chr'], i)
    #         print den_base
    #         pass
    
    # print 1, args.bsg_files[1].get_base_one_based(data['window_chr'], 51511999)

    
    lines = ['%-10s' % s for s in ('chimp', 'denisovan', 'archaic', 'human')]
    snps = list(data['snps'])
    tag_snps = [s for s in snps if args.second_region.in_region_one_based(s['chr'], s['pos'])]

    if args.debug: print data['region_tag'].replace('\t', '.')

    ## is there an individual in this window that's excessively homozygous for tag snps?
    mlen = 0
    roh_ind = -1
    for i in xrange(num_genotypes):
        #homozygous = ''.join(['t' if s['snp_details'][i][2] == 2 else 'f' for s in tag_snps])
        #roh = [m.span(0) for m in re.finditer("(tttt+t)", homozygous)]
        homozygous = ''.join(['T' if s in tag_snps and s['snp_details'][i][2] == 2 else 't' if s['snp_details'][i][1] == 2 or s['snp_details'][i][2] == 2 else '_' for s in snps])
        roh = [m.span(0) for m in re.finditer("(T[Tt]*T[Tt]*T[Tt]*T[Tt]*T)", homozygous)]
        roh_len =  sum([m[1] - m[0] for m in roh]) if len(roh) > 0 else 0

        if args.debug: print homozygous, roh, roh_len
        if args.debug and roh_len > 0: print ''.join([' ' if len([m for m in roh if m[0] <= i <= (m[1]-1)]) == 0 else '*' for i in range(len(homozygous))])
        
        if mlen < roh_len:
            mlen = roh_len
            roh_ind = i
            roh_pos = roh
            pass
        pass

    if mlen == 0: return

    # print 2, args.bsg_files[1].get_base_one_based(data['window_chr'], 51511999), args.bsg_files[0].get_base_one_based(data['window_chr'], 177346887)

    ## go through each roh window
    lookup = BaseLookup('/net/akey/vol1/home/bvernot/tishkoff/primate_sequences/')
    base_output_count = 0
    for s1, s2 in roh_pos:
        start = snps[s1]['pos']
        end = snps[s2-1]['pos']
        #if start == 177333159: continue
        if args.debug: print 'roh:', start, end
        ## all snps in that window (tag or not)
        roh_snps = [s for s in snps if s['pos'] >= start and s['pos'] <= end]
        
        for i in range(start, end+1):
            if args.debug and base_output_count > 200: break
            if not args.regions.in_region_one_based(data['window_chr'], i):
                #print 'skipping'
                continue

            #print i
            human_base = lookup.getBase('hg19', data['window_chr'], i)
            chimp_base = args.bsg_files[0].get_base_one_based(data['window_chr'], i)
            den_base = args.bsg_files[1].get_base_one_based(data['window_chr'], i)

            s = [s for s in snps if s['pos'] == i]
            if len(s) == 1:
                arc_base = s[0]['var']
            elif len(s) > 1:
                print "error = should be only one matching snp!"
                print s
                sys.exit(-1)
                pass
            else:
                arc_base = 'N'
                pass


            not_diff = (human_base, 'N')
            not_diff_base = human_base
            if args.debug: not_diff_base = '_'

            if args.debug: print 'site', i, human_base, arc_base, den_base, chimp_base, ' '.join(['_' if b in not_diff else p for (p, b) in [('r', arc_base), ('d', den_base), ('p', chimp_base)]])
            if args.debug and chimp_base in not_diff and den_base in not_diff and arc_base in not_diff: continue
            
            lines[0] += chimp_base if chimp_base not in not_diff else not_diff_base
            #lines[0] += chimp_base if chimp_base != 'N' else human_base
            
            lines[1] += den_base if den_base not in not_diff else not_diff_base
            #lines[1] += den_base if den_base != 'N' else human_base

            lines[2] += arc_base if arc_base not in not_diff else not_diff_base
            #lines[1] += den_base if den_base != 'N' else human_base

            lines[3] += human_base

            base_output_count += 1
            pass
        pass

    print 4, len(lines[0]) - 10
    for i in range(len(lines)):
        print lines[i]
        pass

    sys.exit()
    return


def calc_visual_genotype(data):

    lines = []
    header = '%s\tsnp\tpos' % (data['region_tag'].replace('\t', '.'))
    header2 = '%s\ttag\tsnp' % (data['region_tag'].replace('\t', '.'))

    ## build up the first few columns
    for i, ind in enumerate(individuals_labels):
        s = "%s\t%s\t%s" % (data['region_tag'].replace('\t', '.'), ind, individuals_pops[i])
        lines.append(s)
        pass

    for s in data['snps']:

        header += '\t%d' % (s['pos'])
        if args.second_region.in_region_one_based(s['chr'], s['pos']):
            header2 += ' 1'
        else:
            #continue
            header2 += ' 0'
            pass

        for i in xrange(num_genotypes):
            if args.debug: print 'looking at ind %d %s;' % (i, individuals_labels[i])

            ## get the genotype settings
            if s['snp_details'][i][0] > 0:
                lines[i] += ' -9'
            else:
                lines[i] += ' %d' % s['snp_details'][i][2]
                pass
            pass
        pass
    
    print header
    print header2
    for i in range(len(lines)):
        print lines[i]
        pass

    return


def calc_roh(data):
    if args.print_header and args.multiple_populations != None:
        print '\t'.join(['chr', 'start', 'end', 'var'] + ['%s_freq' % p for p in args.mpops] + ['%s_haps' % p for p in args.mpops])
        pass
    snp_list = data['snps']

    for next_snp in snp_list:
        f = next_snp['vars'] / next_snp['haps']
        if args.multiple_populations != None:
            freqs = [p['vars'] / p['haps'] if p['haps'] != 0 else 0 for p in next_snp['mpops']]
            haps = [p['haps'] for p in next_snp['mpops']]
            print "\t".join((str(s) for s in [next_snp['chr'], next_snp['pos']-1, next_snp['pos'], next_snp['var']] + freqs + haps))
        else:
            print "%s\t%d\t%d\t%s\t%f" % (next_snp['chr'], next_snp['pos']-1, next_snp['pos'], next_snp['var'], f)
            #print "%s\t%d\t%d\t%f" % (next_snp['chr'], next_snp['pos']-1, next_snp['pos'], min(f, 1-f))
            pass
        pass
    return


def calc_plink(data):

    plink_lines = []
    pops = {}
    header = ''
    
    for i,p in enumerate(list(set(individuals_pops))):
        pops[p] = i
        pass

    for i, ind in enumerate(individuals_labels):
        s = "%s\t%s" % (ind, pops[individuals_pops[i]])
        plink_lines.append([s,s])
        pass

    ## set up virtual genome stuff
    if args.individual_regions_names != None:
        if len(args.individual_regions_names) != len(args.individual_regions):
            print "lengths of ind regions names and ind regions has to be the same!  One to one!"
            sys.exit(-1)
            pass
        for i in range(len(args.individual_regions_names)):
            if args.individual_regions_names[i] not in individuals_labels:
                print "bad individual label!"
                print args.individual_regions_names[i], 'should be in:', individuals_labels
                sys.exit(-1)
                pass
            if args.individual_regions_names[i] not in args.individual_regions[i].filename:
                print "individual label does not match filename:"
                print args.individual_regions_names[i], 'should be in:', args.individual_regions[i].filename
                sys.exit(-1)
                pass
            pass
        

        ## maps an individual's position in individuals_labels to their position in args.individual_regions_names
        reg_mapping = {}
        virt_mapping = {}
        virt_pops = []
        for i in xrange(num_genotypes):
            if individuals_labels[i] in args.individual_regions_names:

                # map each individual to its regions index
                reg_mapping[i] = args.individual_regions_names.index(individuals_labels[i])

                # establish the virtual individual
                if individuals_pops[i] not in virt_pops:
                    virt_pops.append(individuals_pops[i])
                    s = "virt_%s\t%s" % (individuals_pops[i], len(virt_pops) + len(pops))
                    plink_lines.append([s,s])
                    if args.debug: print 'added virtual genome for pop %s:' % individuals_pops[i], plink_lines
                    pass
                
                # map each individual to its virtual individual
                virt_mapping[i] = virt_pops.index(individuals_pops[i])
            else:
                # some individuals don't map to virtual individuals.  :(
                virt_mapping[i] = None
                pass

            if args.debug: print 'processed ind %d:' % i, reg_mapping, virt_mapping, virt_pops
            pass
        pass

    #print pops

    for s in data['snps']:
        if args.debug: print 'outputing plink, looking at snp %s %d;' % (s['chr'], s['pos'])
        if args.individual_regions_names != None: virt_marked = dict((p,False) for p in virt_mapping.values() if p != None)

        header += '\t%s.%d' % (s['chr'], s['pos'])

        for i in xrange(num_genotypes):
            if args.debug: print 'looking at ind %d %s;' % (i, individuals_labels[i])

            ## get the genotype settings
            if s['snp_details'][i][1] == 0 and s['snp_details'][i][2] == 2:
                l1 = '\t2'
                l2 = '\t2'
            elif s['snp_details'][i][1] == 1 and s['snp_details'][i][2] == 1:
                l1 = '\t1'
                l2 = '\t2'
            elif s['snp_details'][i][1] == 2 and s['snp_details'][i][2] == 0:
                l1 = '\t1'
                l2 = '\t1'
            else:
                print 'oh shit you need to set -fm 1!'
                print s['snp_details'][i]
                pass

            ## considering virtual genomes
            if args.individual_regions_names != None:
                if args.debug: print 'virt_marked:', virt_marked
                # only use the virtual genome if 1) this individual maps to a virtual genome, 2) that virtual genome has not been marked for this snp, and 3) this snp falls in that virtual genome's regions
                if virt_mapping[i] != None and not virt_marked[virt_mapping[i]] and args.individual_regions[reg_mapping[i]].in_region_one_based(s['chr'], s['pos']):
                    virt_marked[virt_mapping[i]] = True
                    if args.debug: l1 += individuals_labels[i]
                    if args.debug: l2 += individuals_labels[i]
                    plink_lines[num_genotypes + virt_mapping[i]][0] += l1
                    plink_lines[num_genotypes + virt_mapping[i]][1] += l2
                    pass
                pass
            
            ## mark non-virtual genomes
            plink_lines[i][0] += l1
            plink_lines[i][1] += l2
            pass
        
        ## add missing data for virtual genomes that do not contain this snp
        if args.individual_regions_names != None:
            for i in virt_marked:
                if args.debug: print 'virt genome marked? %d %s' % (i, virt_marked[i])
                if virt_marked[i]: continue
                plink_lines[num_genotypes + i][0] += '\t-9'
                plink_lines[num_genotypes + i][1] += '\t-9'
                pass
            pass
        pass
    
    if args.print_header: print header
    for i in range(len(plink_lines)):
        print plink_lines[i][0]
        print plink_lines[i][1]
        pass

    return

        
def calc_ped_map(data):

    #### this does not work!  it was stolen from calc_plink
    #### http://www.shapeit.fr/pages/m02_formats/pedmap.html

    # PED file
    # The PED file describes the individuals and the genetic data. The PED file corresponding to the example dataset is:
    # FAM1 IND1  0     0     1 0 A A T T 0 0
    # FAM2 IND2  0     0     1 0 A G T C T A
    # FAM3 TRIOF 0     0     1 0 A G T C A T
    # FAM4 TRIOM 0     0     2 0 A G T C A T
    # FAM5 TRIOC TRIOF TRIOM 1 0 A A C T A T
    # FAM6 DUOP  0     0     2 0 G A T C A A
    # FAM7 DUOC  DUOP  0     2 0 A A T C A A
    # This file can be SPACE or TAB delimited. Each line corresponds to a single individual. The first 6 columns are:
    # Family ID [string]
    # Individual ID [string]
    # Father ID [string]
    # Mother ID [string]
    # Sex [integer]
    # Phenotype [float]
    # Columns 7 & 8 code for the observed alleles at SNP1, columns 9 & 10 code for the observed alleles at SNP2, and so on. 
    # Missing data are coded as "0 0" as for example for SNP3 of IND1. This file should have N lines and 2L+6 columns, 
    # where N and L are the numbers of individuals and SNPs contained in the dataset respectively.
    # Each individual must have an unique ID containing only alphanumeric characters.

    # MAP file
    # The MAP file describes the SNPs. The MAP file corresponding to the example dataset is:
    # 7 SNP1 0 123
    # 7 SNP2 0 456
    # 7 SNP3 0 789
    # This file can be SPACE or TAB delimited. Each line corresponds to a SNP. The 4 columns are:
    # Chromosome number [integer]
    # SNP ID [string]
    # SNP genetic position (cM) [float]
    # SNP physical position (bp) [integer]
    # This file should have L lines and 4 columns, where L is the number of SNPs contained in the dataset.
    # Each SNP must have a unique physical position. All the SNPs must be ordered by physical position.

    ped_lines = []
    map_lines = cStringIO.StringIO()
    
    for i, ind in enumerate(individuals_labels):
        s = cStringIO.StringIO()
        s.write("%s %s %d %d %d %d " % (ind, ind, 0, 0, is_male[i]+1, 0))
        ped_lines.append(s)
        pass

    for s_i,s in enumerate(data['snps']):
        if args.debug: print 'outputing ped, looking at snp %s %d;' % (s['chr'], s['pos'])

        map_lines.write( '%s snp%d %d %d\n' % (s['chr'][3:], s_i, 0, s['pos']) )

        for ind in xrange(num_genotypes):
            if args.debug: print 'looking at ind %d %s;' % (ind, individuals_labels[ind])

            ## get the genotype settings
            gt = ''
            for x in range(s['snp_details'][ind][0]):
                gt += '0 '
                pass
            for x in range(s['snp_details'][ind][1]):
                gt += s['ref'] + ' '
                pass
            for x in range(s['snp_details'][ind][2]):
                gt += s['var'] + ' '
                pass
            if len(gt) != 4:
                print 'incorrect gt length?'
                print gt
                print s['snp_details'][ind]
                sys.exit(-1)
                pass

            ped_lines[ind].write(gt)
            pass
        
        pass

    f = open(args.output_file.name + '.ped', 'w')
    for i in xrange(num_genotypes):
        f.write( ped_lines[i].getvalue() + '\n' )
        pass
    
    f = open(args.output_file.name + '.map', 'w')
    f.write( map_lines.getvalue() )

    return

        
def calc_snp_counts(data):
    snp_count = sum([1 for i in data['snps']])

    snp_per_base = 0 if data['region_length'] == 0 else snp_count / data['region_length']
    print '\t'.join([str(i) for i in (data['region_tag'], snp_count, data['region_length'], snp_per_base)])
    return

    
def calc_site_freq(data):
    snp_list = data['snps']

    summary_hist = {}

    if args.summary_only and args.print_header:
        ## print header
        print '\t'.join([str(i) for i in [data['region_tag'], 'counts', 'num_sites'] + list(sorted(args.tag.keys()))])
        pass

    for next_snp in snp_list:

        if next_snp['mrca_known'] or args.folded_sfs == False:
            # by default get the full SFS - but we need ancestry
            freq = next_snp['vars'] if next_snp['mrca'] == 0 else next_snp['refs']
        else:
            freq = min(next_snp['refs'], next_snp['vars'])
            pass

        if args.summary_only:
            if freq not in summary_hist: summary_hist[freq] = 0
            summary_hist[freq] += 1
        else:
            print '%d\t%d' % (freq, next_snp['haps'])
            pass
        pass

    if args.summary_only:
        for f in sorted(summary_hist.keys()):
            print '\t'.join([str(i) for i in [data['region_tag'], f, summary_hist[f]] + [args.tag[k] for k in sorted(args.tag.keys())]])
            pass
        pass
    
    return

def choose(n, k):
    if 0 <= k <= n:
        p = 1
        for t in xrange(min(k, n - k)):
            p = (p * (n - t)) // (t + 1)
            pass
        return p
    return 0


## THIS IS WEIRD - I DON'T LIKE THE VALUES GIVEN BY THIS COMPUTATION
## ALSO, RIGHT NOW THE CODE HAS SEVERAL VERSIONS INTERLACED AND NOT CLEANED UP
## from page 126 of Bruce Weir's Genetic Data Analysis II
## calculate the composite LD between snps i and j, in the subset of individuals *subset*
def calc_composite_ld_freal(snps, a, b, results = None, subset = None):
    print 'debug composite ld', a, b, subset
    print 'debug composite ld', snps[a]
    print 'debug composite ld', snps[b]

    if subset == None:
        A = snps[a]['snp_details']
        B = snps[b]['snp_details']
    else:
        A = [snps[a]['snp_details'][s] for s in subset]
        B = [snps[b]['snp_details'][s] for s in subset]
        pass

    print 'debug composite ld', A
    print 'debug composite ld', B
    
    ## filter out inds where either snp has missing data
    A_filt = [s for i,s in enumerate(A) if A[i][0] == 0 and B[i][0] == 0]
    B_filt = [s for i,s in enumerate(B) if A[i][0] == 0 and B[i][0] == 0]
    n = len(A_filt)

    print 'debug composite ld', n
    print 'debug composite ld', A_filt
    print 'debug composite ld', B_filt

    # hom_hom = sum([1 for i in xrange(n) if sorted((max(A_filt[i][1:3]), max(B_filt[i][1:3]))) == [2,2]])
    # het_hom = sum([1 for i in xrange(n) if sorted((max(A_filt[i][1:3]), max(B_filt[i][1:3]))) == [1,2]])
    # het_het = sum([1 for i in xrange(n) if sorted((max(A_filt[i][1:3]), max(B_filt[i][1:3]))) == [1,1]])

    hom_hom = sum([1 for i in xrange(n) if (A_filt[i][1], B_filt[i][1]) == (2,2)])
    het_hom = sum([1 for i in xrange(n) if (A_filt[i][1], B_filt[i][1]) == (1,2)])
    hom_het = sum([1 for i in xrange(n) if (A_filt[i][1], B_filt[i][1]) == (2,1)])
    het_het = sum([1 for i in xrange(n) if (A_filt[i][1], B_filt[i][1]) == (1,1)])

    p_A = sum([s[1] for s in A_filt]) / n / 2
    p_B = sum([s[1] for s in B_filt]) / n / 2
    #if p_A > .5: p_A = 1-p_A
    #if p_B > .5: p_B = 1-p_B

    print '2 * %d + %d + %d + %d / 2' % (hom_hom, het_hom, hom_het, het_het)
    n_AB = 2 * hom_hom + het_hom + hom_het + het_het / 2

    print '%f / %d - 2 * %f * %f' % (n_AB, n, p_A, p_B)

    print snps[a]['vars'], snps[a]['haps'],  snps[a]['vars'] / snps[a]['haps']

    delta_AB = n_AB / n - 2 * p_A * p_B

    print 'comp_ld_result', delta_AB * delta_AB

    hom_hom = sum([1 for i in xrange(n) if (A_filt[i][2], B_filt[i][2]) == (2,2)])
    het_hom = sum([1 for i in xrange(n) if (A_filt[i][2], B_filt[i][2]) == (1,2)])
    hom_het = sum([1 for i in xrange(n) if (A_filt[i][2], B_filt[i][2]) == (2,1)])
    het_het = sum([1 for i in xrange(n) if (A_filt[i][2], B_filt[i][2]) == (1,1)])

    p_A = sum([s[2] for s in A_filt]) / n / 2
    p_B = sum([s[2] for s in B_filt]) / n / 2

    n_AB = 2 * hom_hom + het_hom + hom_het + het_het / 2

    delta_AB = n_AB / n - 2 * p_A * p_B
    
    from itertools import imap
    
    def pearsonr(x, y):
        # Assume len(x) == len(y)
        n = len(x)
        sum_x = float(sum(x))
        sum_y = float(sum(y))
        sum_x_sq = sum(map(lambda x: pow(x, 2), x))
        sum_y_sq = sum(map(lambda x: pow(x, 2), y))
        psum = sum(imap(lambda x, y: x * y, x, y))
        num = psum - (sum_x * sum_y/n)
        den = pow((sum_x_sq - pow(sum_x, 2) / n) * (sum_y_sq - pow(sum_y, 2) / n), 0.5)
        if den == 0: return 0
        return num / den
    
    # from scipy.stats.stats import pearsonr
    print 'pearson', pearsonr([A_filt[i][1] for i in xrange(n)], [B_filt[i][1] for i in xrange(n)])

    return delta_AB



# from itertools import imap
# def pearsonr(x, y):
#     # Assume len(x) == len(y)
#     n = len(x)
#     sum_x = float(sum(x))
#     sum_y = float(sum(y))
#     sum_x_sq = sum(map(lambda x: pow(x, 2), x))
#     sum_y_sq = sum(map(lambda x: pow(x, 2), y))
#     psum = sum(imap(lambda x, y: x * y, x, y))
#     num = psum - (sum_x * sum_y/n)
#     den = pow((sum_x_sq - pow(sum_x, 2) / n) * (sum_y_sq - pow(sum_y, 2) / n), 0.5)
#     if den == 0: return 0
#     return num / den

from snpdist import pearson_def, pearsonr

## pearson correlation coefficient ^2 (R^2)
def calc_pearson_freal(snps, a, b, results = None, subset = None, cache = None, check_missingness = True):

    if cache != None and (a,b) in cache:
        return cache[(a,b)]

    if subset == None:
        A = snps[a]['snp_details']
        B = snps[b]['snp_details']
    else:
        def get_subset():
            A = [snps[a]['snp_details'][s] for s in subset]
            B = [snps[b]['snp_details'][s] for s in subset]
            return (A,B)
        (A, B) = get_subset()
        pass

    ## filter out inds where either snp has missing data

    def get_filt():
        n = len(A)
        A_filt = [A[i] for i in xrange(n) if A[i][0] == 0 and B[i][0] == 0]
        B_filt = [B[i] for i in xrange(n) if A[i][0] == 0 and B[i][0] == 0]
        return (A_filt, B_filt)
    if check_missingness:
        (A, B) = get_filt()
        pass

    n = len(A)

    def get_refs():
        return ([A[i][1] for i in xrange(n)], [B[i][1] for i in xrange(n)])
    (A_refs, B_refs) = get_refs()

    # from scipy.stats.stats import pearsonr
    p = pearson_def(A_refs, B_refs)

    if cache != None:
        cache[(a,b)] = p * p
        cache[(b,a)] = p * p
        pass

    return p * p



def calc_tajd_freal(snp, results = None):
    
    # taken from:
    # http://ocw.mit.edu/courses/health-sciences-and-technology/hst-508-quantitative-genomics-fall-2005/study-materials/tajimad1.pdf
    # http://en.wikipedia.org/wiki/Tajima's_D
    # http://www.genetics.org/content/123/3/585.full.pdf


    n = snp['haps']
    pi = 2 * snp['vars'] * snp['refs'] / n / n
    pi = snp['vars'] * snp['refs'] / choose(n, 2)

    a1 = sum([1/i for i in range(1,n)])
    a2 = sum([1/i**2 for i in range(1,n)])
    b1 = (n+1) / 3 / (n-1)
    b2 = 2 * (n**2 + n + 3) / 9 / n / (n-1)
    c1 = b1 - 1/a1
    c2 = b2 - (n+2)/a1/n + a2/a1**2
    e1 = c1 / a1
    e2 = c2 / (a1**2 + a2)

    if results == None:
        results = {'pi' : 0, 'S' : 0}
        pass
    pi += results['pi']
    S = results['S'] + 1

    V = math.sqrt(e1 * S + e2 * S * (S-1))
    D = (pi - S / a1) / V
    results['D'] = D
    results['S'] = S
    results['V'] = V
    results['e1'] = e1
    results['e2'] = e2
    results['a1'] = a1
    results['a2'] = a2
    results['c1'] = c1
    results['c2'] = c2
    results['b1'] = b1
    results['b2'] = b2
    results['pi'] = pi

    return results

def calc_tajd(data):

    snp_list = data['snps']
    results = None
    
    for next_snp in snp_list:
        results = calc_tajd_freal(next_snp, results)
        if not args.summary_only: print '\t'.join([str(s) for s in (next_snp['chr'], next_snp['pos']-1, next_snp['pos'], results['D'])])
        pass

    #if args.summary_only: print '\t'.join([str(s) for s in (fst_sum, fst_count, fst_sum / fst_count, 1-Hs_sum/Ht_sum)])
    print '\t'.join(sorted(results))
    print '\t'.join([str(results[k]) for k in sorted(results)])
        
    return


def calc_pi_freal(snp, results = None):
    pi = 2 * snp['vars'] / snp['haps'] * snp['refs'] / snp['haps']
    if results == None:
        return {'site_pi' : pi, 'pi_sum' : pi}
    results['site_pi'] = pi
    results['pi_sum'] += pi
    return results

def calc_pi(data):

    snp_list = data['snps']
    
    if args.resample: pi_list = []
    
    # do that calculation!
    pi_sum = 0
    num_pi = 0

    # print 'snp list', snp_list

    for next_snp in snp_list:
        pi = 2 * next_snp['vars'] / next_snp['haps'] * next_snp['refs'] / next_snp['haps']
        if args.resample: pi_list.append(pi)
        pi_sum += pi
        num_pi += 1
        if not args.summary_only: print next_snp['chr'], next_snp['pos'], next_snp['vars'], next_snp['refs'], pi
        if args.debug: print 'debug_pi:', next_snp['chr'], next_snp['pos'], pi
        # print pi
        pass

    if args.summary_only:
        region_length = data['region_length']

        if region_length > 0: regions_pi = pi_sum / region_length
        else: regions_pi = 'NA'
        
        if args.resample:
            # print 'resample'
            num_resamples = 10000
            sys.stderr.write('about to resample %d times\n' % num_resamples)
            # get the number of variable sites to sample (* num_resamples)
            num_variable_sites_samples = binomial(region_length, len(pi_list) / region_length, num_resamples)
            resampled_pis = [0] * num_resamples
            for i,num_variable_sites in enumerate(num_variable_sites_samples):
                # sample num_variable_sites sites, sum the pi's
                samp = sample_wr(pi_list, num_variable_sites)
                # print '\t'.join([str(s) for s in samp])
                sample_pi_sum = sum(samp)
                sample_regions_pi = sample_pi_sum / region_length
                resampled_pis[i] = sample_regions_pi
                # print "regions_pi_sample:", args.regions.filename, "%f/%d=" % (pi_sum, region_length), regions_pi, "debug:", pi_sum, num_variable_sites
                if args.print_resamples:
                    print '\t'.join([str(s) for s in (data['region_tag'], sample_pi_sum, region_length, sample_regions_pi, sample_regions_pi, sample_regions_pi, 'resample')])
                    pass
                pass
            resampled_pis.sort()
            print '\t'.join([str(s) for s in (data['region_tag'], pi_sum, region_length, regions_pi, resampled_pis[int(num_resamples * .025)], resampled_pis[int(num_resamples * .975)], 'result')])
            #add_to_db(['region_tag', 'pi_sum', 'region_length', 'regions_pi', 'ci_low', 'ci_high'] + args.tag.keys(), [data['region_tag'], pi_sum, region_length, regions_pi, resampled_pis[int(num_resamples * .025)], resampled_pis[int(num_resamples * .975)]] + args.tag.values())
        else:
            print '\t'.join([str(s) for s in (data['region_tag'], pi_sum, region_length, regions_pi)])
            pass
        pass
    return


def calc_fst_freal(snp, results = None, drop_singletons = False):
    if drop_singletons and snp['mpops'][0]['vars'] + snp['mpops'][1]['vars'] == 1:
        return results
    p1 = snp['mpops'][0]['vars'] / snp['mpops'][0]['haps']
    p2 = snp['mpops'][1]['vars'] / snp['mpops'][1]['haps']
    haps1 = snp['mpops'][0]['haps']
    haps2 = snp['mpops'][1]['haps']
    tot_haps = haps1 + haps2
    p_bar = p1 * haps1 / tot_haps + p2 * haps2 / tot_haps
    Hs = haps1 / tot_haps * 2 * p1 * (1-p1) + haps2 / tot_haps * 2 * p2 * (1-p2)
    Ht = 2 * p_bar * (1-p_bar)
    Fst = 1 - (1 if Ht == 0 else Hs/Ht)
    if results == None:
        return {'site_fst' : Fst, 'fst' : Fst, 'fst_sum' : Fst, 'fst_count' : 1, 'hs' : Hs, 'ht' : Ht, 'fst_h' : Fst}
    results['site_fst'] = Fst
    results['fst_sum'] += Fst
    results['fst_count'] += 1
    results['fst'] = results['fst_sum'] / results['fst_count']
    results['hs'] += Hs
    results['ht'] += Ht
    results['fst_h'] = 1 - (1 if results['ht'] == 0 else results['hs'] / results['ht'])
    return results

def calc_fst(data):

    snp_list = data['snps']
    results = None
    
    for next_snp in snp_list:
        results = calc_fst_freal(next_snp, results, drop_singletons = True)
        if not args.summary_only: print '\t'.join([str(s) for s in (next_snp['chr'], next_snp['pos']-1, next_snp['pos'], results['site_fst'], results['fst'])])
        pass

    #if args.summary_only: print '\t'.join([str(s) for s in (fst_sum, fst_count, fst_sum / fst_count, 1-Hs_sum/Ht_sum)])
    print '\t'.join([str(s) for s in (results['fst_sum'], results['fst_count'], \
                                          results['fst_sum'] / results['fst_count'], \
                                          results['fst'], 1-results['hs']/results['ht'])])
        
    return


def calc_mk(data):

    print "NOT SURE THAT THIS PROPERLY CALCULATES THE REGION LENGTH, OR PROPERLY FILTERS CHIMP/REGULAR SNPS (TAKING INTO ACCOUNT EXCLUDES, INTERSECTIONS, AND CHIMP_MAPPED)"
    sys.exit(-1)

    snp_list = data['snps']
    chimp_snp_list = data['chimp_snps']
    #mapped = data['region_length']
    snp_regions = args.snp_regions

    fixed_neu = 0
    fixed_sel = 0
    poly_neu = 0
    poly_sel = 0
    
    for next_snp in snp_list:
        chr = next_snp['chr']
        pos = next_snp['pos']
        if args.neutral_regions.in_region_one_based(chr, pos): poly_neu += 1
        if args.selected_regions.in_region_one_based(chr, pos): poly_sel += 1
        snp_regions.bases[snp_regions.chr_offset[chr] + pos - 1] = True
        pass
    
    for next_snp in chimp_snp_list:
        chr = next_snp['chr']
        pos = next_snp['pos']
        # make sure this isn't the site of a polymorphic snp
        if snp_regions.bases[snp_regions.chr_offset[chr] + pos - 1]: continue
        if args.neutral_regions.in_region_one_based(chr, pos): fixed_neu += 1
        if args.selected_regions.in_region_one_based(chr, pos): fixed_sel += 1
    	pass
    
    dnpn = 0 if poly_sel == 0 else (fixed_sel / poly_sel)
    dsps = 0 if poly_neu == 0 else (fixed_neu / poly_neu)
    NI = (0 if fixed_sel == 0 or poly_neu == 0 else (fixed_neu * poly_sel) / (fixed_sel * poly_neu))
    alpha = 1 - NI
    if data['is_window']:
        sel_len = args.selected_regions.amount_in_region(data['window_chr'], data['window_start'], data['window_end'])
        neu_len = args.neutral_regions.amount_in_region(data['window_chr'], data['window_start'], data['window_end'])
        if not args.resample: print '\t'.join([str(i) for i in (data['region_tag'], NI, poly_neu, poly_sel, fixed_neu, fixed_sel, dnpn, dsps, dnpn - dsps, alpha, sel_len, neu_len, sel_len + neu_len, min(sel_len, neu_len))])
    else:
        sel_len = args.selected_regions.total_length()
        neu_len = args.neutral_regions.total_length()
        print "mk neutrality index:", NI
        print '(%f * %f) / (%f * %f) = %f' % (fixed_neu,  poly_sel, fixed_sel, poly_neu, NI)
        pass
    
    #print args.neutral_regions.total_length()
    if args.resample:
        num_resamples = 1000
        if not data['is_window']: sys.stderr.write('about to resample %d times\n' % num_resamples)
        resampled_mk = [0] * num_resamples
        if neu_len > 0 and sel_len > 0:
            fixed_neu_samples = binomial(neu_len, fixed_neu / neu_len, num_resamples)
            fixed_sel_samples = binomial(sel_len, fixed_sel / sel_len, num_resamples)
            poly_neu_samples = binomial(neu_len, poly_neu / neu_len, num_resamples)
            poly_sel_samples = binomial(sel_len, poly_sel / sel_len, num_resamples)
            for i in range(num_resamples):
                resampled_mk[i] = (sys.maxint if fixed_sel_samples[i] == 0 or poly_neu_samples[i] == 0 else (fixed_neu_samples[i] * poly_sel_samples[i]) / (fixed_sel_samples[i] * poly_neu_samples[i]))
                pass
            resampled_mk.sort()
            pass
        if data['is_window']:
            print '\t'.join([str(i) for i in (data['region_tag'], NI, poly_neu, poly_sel, fixed_neu, fixed_sel, min(sel_len, neu_len), resampled_mk[int(num_resamples * .025)], resampled_mk[int(num_resamples * .975)])])
        else:
            print "CI:", resampled_mk[int(num_resamples * .025)], resampled_mk[int(num_resamples * .975)]
            
    
    return

def calc_mut_rate(data):

    #print "NOT SURE THAT THIS PROPERLY CALCULATES THE REGION LENGTH (TAKING INTO ACCOUNT EXCLUDES, INTERSECTIONS, AND CHIMP_MAPPED)"
    #sys.exit(-1)

    #snp_list = data['snps']
    chimp_snp_list = data['chimp_snps']

    mapped = chimp_mapped_regions.total_length()
    human_mapped = data['region_length']
    
    if args.resample: print 'resampling not yet implemented!' #pi_list = []

    num_chimp_snps = sum(1 for e in chimp_snp_list)

    if mapped > 0:
        ratio = num_chimp_snps / mapped
    else:
        ratio = 0
        pass

    if args.s4:
        print '\t'.join([str(i) for i in [data['region_tag'], num_chimp_snps, mapped, ratio, ratio / (12000000 / 25), human_mapped, data['window_name']]])
    else:
        print '\t'.join([str(i) for i in [data['region_tag'], num_chimp_snps, mapped, ratio, ratio / (12000000 / 25), human_mapped]])
        pass

    #add_to_db(['region_tag', 'num_chimp_snps', 'mapped', 'human_mapped', 'ratio', 'mut_rate'] + args.tag.keys(), [data['region_tag'], num_chimp_snps, mapped, human_mapped, ratio, ratio / (12000000 / 25)] + args.tag.values())
    return


def calc_reg_len(data):

    length = data['region_length']

    if args.s4:
        print '\t'.join([str(i) for i in [length, data['window_name']]])
    else:
        print '\t'.join([str(i) for i in [data['window_chr'], data['window_start'], data['window_end'], length]])
        pass
    return

def calc_reg_ratio(data):

    length = data['region_length']
    length2 = data['second_region_length']

    if length2 == 0:
        ratio = 0
    else:
        ratio = length2 / length
        pass

    if args.s4:
        print '\t'.join([str(i) for i in [length2, length, ratio, data['window_name']]])
    else:
        print '\t'.join([str(i) for i in [length2, length, ratio]])
        pass
    return


###########################
# choose calc fn (these fns should be put in another class, and selected during the args parsing section)
###########################

results_fn = None

if args.func == 'pi':
    calc_fn = calc_pi
elif args.func == 'mutrate':
    calc_fn = calc_mut_rate
elif args.func == 'plink':
    calc_fn = calc_plink
elif args.func == 'ped_map':
    calc_fn = calc_ped_map
elif args.func == 'visual_genotype':
    calc_fn = calc_visual_genotype
elif args.func == 'visual_fasta':
    calc_fn = calc_visual_fasta
elif args.func == 'freqs':
    calc_fn = calc_freqs
elif args.func == 'sites':
    calc_fn = calc_sites
elif args.func == 'stats':
    calc_fn = calc_summary_stats
elif args.func == 'stats-win':
    calc_fn = calc_summary_stats_windowed
    results_fn = report_summary_stats
elif args.func == 'msfile':
    calc_fn = calc_msfile
elif args.func == 'roh':
    calc_fn = calc_roh
elif args.func == 'dist_matrix':
    calc_fn = calc_dist_matrix
elif args.func == 'mk':
    calc_fn = calc_mk
elif args.func == 'snpcount':
    calc_fn = calc_snp_counts
elif args.func == 'site_freq':
    calc_fn = calc_site_freq
elif args.func == 'fst':
    calc_fn = calc_fst
elif args.func == 'reglen':
    calc_fn = calc_reg_len
elif args.func == 'region-ratio':
    calc_fn = calc_reg_ratio
elif args.func == 'tmrca':
    calc_fn = calc_tmrca
elif args.func == 'tmrca_subsets':
    calc_fn = calc_tmrca_subsets
elif args.func == 'region_s_star':
    calc_fn = calc_region_s_star
elif args.func == '':
    calc_fn = None
    pass


###########################
# run fn on whole genome
###########################

if ('window_file' not in args or args.window_file == None) and ('ms_win_params' not in args or args.ms_win_params == None):
    data = {}
    data['snps'] = read_snps_gen(read_snp)
    data['chimp_snps'] = read_snps_gen(read_chimp)
    if args.regions != None:
        data['region_length'] = args.regions.total_length() * args.base_sample_probability
    else:
        # this isn't great, because it doesn't take into account unmapped regions!
        data['region_length'] = args.baselookup.chr_lens['genome']
        pass
    if args.second_region != None:
        data['second_region_length'] = args.second_region.total_length()
        pass
    data['is_window'] = False
    data['region_tag'] = args.data_name

    calc_fn(data)
    sys.exit(0)
    pass

###########################
# run fn on windows
###########################

file_chr = ''
if 'varfile' in args: next_snp = read_snp()
if 'chimp_snps' in args and args.chimp_snps != None: next_snp_chimp = read_chimp()
list_of_snps = list()
list_of_snps_chimp = list()
prev_chr = 'chr1'
prev_start = 0
prev_end = 0

def next_window():

    if args.ms:
        for window_chr in args.ms_sim_labels:
            window_start = 0
            window_end = window_start + args.ms_win_params[0]
            while window_end <= args.ms_reglen:
                if args.limit != None and args.snp_count > args.limit:
                    return
                yield [window_chr, window_start, window_end, None, None]
                window_start += args.ms_win_params[1]
                window_end = window_start + args.ms_win_params[0]
                pass
            pass
        pass
    elif args.merge_windows_by_fourth_column:
        current_window_key = None
        for line in args.window_file:
            print 'debug window', line
            if len(line.strip().split()) == 0: continue
            [window_chr, window_start, window_end, window_key] = line.strip().split()[:4]
            if window_key != current_window_key:
                if current_window_key != None:
                    if args.limit != None and args.snp_count > args.limit:
                        return
                    yield [current_window_chr, current_window_start, current_window_end, current_window_key, window_hash]
                    pass
                current_window_chr = window_chr
                current_window_start = int(window_start)
                current_window_end = int(window_end)
                current_window_key = window_key
                window_hash = '%s,%s' % (window_start, window_end)
            else:
                current_window_end = int(window_end)
                window_hash += ' %s,%s' % (window_start, window_end)
                continue
                pass
            pass
        # one last time (for the last line)
        yield [current_window_chr, current_window_start, current_window_end, current_window_key, window_hash]
        pass
    else:
        for line in args.window_file:
            if len(line.strip().split()) == 0: continue
            [window_chr, window_start, window_end] = line.strip().split()[:3]
            window_start = int(window_start)
            window_end = int(window_end)
            if args.limit != None and args.snp_count > args.limit:
                return
            yield [window_chr, window_start, window_end, None, None]
            pass
        pass
    return

def site_lte_site(c1, p1, c2, p2):
    if args.ms:
        return (c1 == c2 and p1 <= p2) or c1 < c2
    return (c1 == c2 and p1 <= p2) or args.baselookup.chr_nums[c1] < args.baselookup.chr_nums[c2]

def site_in_region(c, p, rc, rs, re):
    return c == rc and p > rs and p <= re

# calc for a particular window
persistent_data = {}
for [window_chr, window_start, window_end, window_key, window_hash ] in next_window():
    
    if args.debug_2: print 'starting loop:', window_chr, "%d-%d" % (window_start, window_end), snp_str(list_of_snps)
    data = {}

    # make sure we're not going backwards (not needed for ms files)
    if not args.ms and (args.baselookup.chr_nums[window_chr] < args.baselookup.chr_nums[prev_chr] \
       or window_start < prev_start \
       or window_end < prev_end):
        print "current window lies before a prior window!"
        print "current", window_chr, window_start, window_end
        print "prior", prior_chr, prior_start, prior_end
        sys.exit(-1)
        pass
    
    #print args.ms, args.base_sample_probability_normal
    #print args.ms, args.base_sample_probability_normal != None
    if args.regions != None:
        data['region_length'] = args.regions.amount_in_region(window_chr, window_start, window_end) * args.base_sample_probability
    elif args.ms and 'ms_region_length' in args:
        data['region_length'] = args.ms_region_length
    else:
        data['region_length'] = window_end - window_start
        pass


    # remove snps that lie before the window
    list_of_snps[:] = [snp for snp in list_of_snps if snp['chr'] == window_chr and snp['pos'] > window_start]
    list_of_snps_chimp[:] = [snp for snp in list_of_snps_chimp if snp['chr'] == window_chr and snp['pos'] > window_start]

    if args.debug_2: print 'removed snps:', window_chr, "%d-%d" % (window_start, window_end), snp_str(list_of_snps)

    if 'varfile' in args:
        # add new snps, make sure they're in the right window and chr, skip any snps that come before this window
        if args.debug_2: print "about to add snps; current snp:", next_snp
        while next_snp and site_lte_site(next_snp['chr'], next_snp['pos'], window_chr, window_start):
            if args.debug_3: print "  skipping snp:", window_chr, next_snp
            next_snp = read_snp()
            pass
        while next_snp and site_in_region(next_snp['chr'], next_snp['pos'], window_chr, window_start, window_end):
            if args.debug_2: print "  adding snp:", next_snp
            list_of_snps.append(next_snp)
            next_snp = read_snp()
            pass
        if args.debug_2: print 'added new snps:', window_chr, "%d-%d" % (window_start, window_end), snp_str(list_of_snps)
        pass

    # add new primate snps
    if 'chimp_snps' in args and args.chimp_snps != None:
        if args.debug_2: print "about to add primate snps; current snp:", next_snp_chimp
        while next_snp_chimp \
                  and ((next_snp_chimp['chr'] == window_chr and next_snp_chimp['pos'] <= window_start) \
                       or (args.baselookup.chr_nums[next_snp_chimp['chr']] < args.baselookup.chr_nums[window_chr])):
            if args.debug_3: print "  skipping primate snp:", window_chr, next_snp_chimp
            next_snp_chimp = read_chimp()
            pass
        while next_snp_chimp \
                  and next_snp_chimp['chr'] == window_chr \
                  and next_snp_chimp['pos'] > window_start \
                  and next_snp_chimp['pos'] <= window_end:
            if args.debug_2: print "  adding primate snp:", next_snp_chimp
            list_of_snps_chimp.append(next_snp_chimp)
            next_snp_chimp = read_chimp()
            pass
        if args.debug_2: print 'added new primate snps:', window_chr, "%d-%d" % (window_start, window_end), snp_str(list_of_snps_chimp)
        pass

    if 'varfile' in args: data['snps'] = list_of_snps
    if 'chimp_snps' in args and args.chimp_snps != None: data['chimp_snps'] = list_of_snps_chimp
    if args.second_region != None:
        data['second_region_length'] = args.second_region.amount_in_region(window_chr, window_start, window_end)
        pass
    data['is_window'] = True
    data['window_chr'] = window_chr
    data['window_start'] = window_start
    data['window_end'] = window_end
    data['region_tag'] = window_key if window_key != None else ' '.join([str(i) for i in (data['window_chr'], data['window_start'], data['window_end'])])
    data['region_tag2'] = window_hash if window_hash != None else ','.join([str(i) for i in (data['window_start'], data['window_end'])])
    data['persist'] = persistent_data

    calc_fn(data)

    pass

if results_fn != None:
    results_fn(data)
    pass
