#!/usr/bin/env python

"""
Contains argument subparsers
"""

import textwrap
import logging
import sys


# Create logger:
logger = logging.getLogger(f'hybpiper.assemble.{__name__}')

title = textwrap.dedent(
    fr"""
                                                     T
                                                        T
                                         C  G
 _    _            _       _____      T        G        A
| |  | |          | |     |  _  \  A              A  A
| |__| | __    __ | |___  | |_| |  _   _____   _____   _____
|  __  | \ \  / / |  _  \ |  ___/ | | |  _  \ |  _  | |  _  \
| |  | |  \ \/ /  | |_| | | |     | | | |_| | |  __/  | |  --
|_|  |_|   \  /   |_____/ |_|     |_| |  ___/ |_____| |_|
           / /                        | |
          /_/                         |_|

    """
)


sys.stderr.write(title)
sys.stderr.flush()


def add_assemble_parser(subparsers):
    """
    Parser for the main assembly stage of HybPiper i.e. assemble.py.

    :param argparse._SubParsersAction subparsers: subparsers object to add parser(s) to
    :return None: no return value specified; default is None
    """

    parser_assemble = subparsers.add_parser('assemble', help='Assemble gene, intron, and supercontig sequences')
    parser_assemble.add_argument('--readfiles', '-r',
                                 nargs='+',
                                 help='One or more read files to start the pipeline. If exactly two are specified, '
                                      'will assume it is paired Illumina reads.',
                                 required=True)
    group_1 = parser_assemble.add_mutually_exclusive_group(required=True)
    group_1.add_argument('--targetfile_dna', '-t_dna',
                         dest='targetfile_dna',
                         default=False,
                         help='FASTA file containing DNA target sequences for each gene. If there are multiple '
                              'targets for a gene, the id must be of the form: >Taxon-geneName')
    group_1.add_argument('--targetfile_aa', '-t_aa',
                         dest='targetfile_aa',
                         default=False,
                         help='FASTA file containing amino-acid target sequences for each gene. If there are multiple '
                              'targets for a gene, the id must be of the form: >Taxon-geneName')
    group_2 = parser_assemble.add_mutually_exclusive_group()
    group_2.add_argument('--bwa',
                         dest='bwa',
                         action='store_true',
                         help='Use BWA to search reads for hits to target. Requires BWA and a target file that is '
                              'nucleotides!',
                         default=False)
    group_2.add_argument('--diamond',
                         dest='diamond',
                         action='store_true',
                         help='Use DIAMOND instead of BLASTx.',
                         default=False)
    parser_assemble.add_argument('--diamond_sensitivity',
                                 choices=['mid-sensitive', 'sensitive', 'more-sensitive', 'very-sensitive',
                                          'ultra-sensitive'],
                                 help='Use the provided sensitivity for DIAMOND searches.',
                                 default=False)
    parser_assemble.add_argument('--start_from',
                                 choices=['map_reads', 'distribute_reads', 'assemble_reads', 'exonerate_contigs'],
                                 help='Start the pipeline from the given step. Note that this relies on the presence '
                                      'of output files for previous steps, produced by a previous run attempt. '
                                      'Default is map_reads',
                                 dest='start_from',
                                 default='map_reads')
    parser_assemble.add_argument('--cpu',
                                 type=int,
                                 default=0,
                                 help='Limit the number of CPUs. Default is to use all cores available.')
    parser_assemble.add_argument('--evalue',
                                 type=float,
                                 default=1e-4,
                                 help='e-value threshold for blastx hits, default: %(default)s')
    parser_assemble.add_argument('--max_target_seqs',
                                 type=int,
                                 default=10,
                                 help='Max target seqs to save in BLASTx search, default: %(default)s')
    parser_assemble.add_argument('--cov_cutoff',
                                 default=8,
                                 help='Coverage cutoff for SPAdes. Default is: %(default)s')
    parser_assemble.add_argument('--single_cell_assembly',
                                 action='store_true',
                                 dest='spades_single_cell',
                                 default=False,
                                 help='Run SPAdes assemblies in MDA (single-cell) mode. Default is False')
    parser_assemble.add_argument('--kvals', nargs='+',
                                 help='Values of k for SPAdes assemblies. SPAdes needs to be compiled to handle '
                                      'larger k-values! Default is auto-detection by SPAdes.',
                                 default=None)
    parser_assemble.add_argument('--thresh',
                                 type=int,
                                 help='Percent identity threshold for retaining Exonerate hits. Default is 55, '
                                      'but increase this if you are worried about contaminant sequences.', default=55)
    parser_assemble.add_argument('--paralog_min_length_percentage',
                                 default=0.75,
                                 type=float,
                                 help='Minimum length percentage of a contig Exonerate hit vs reference protein '
                                      'length for a paralog warning and sequence to be generated. Default is %('
                                      'default)s')
    parser_assemble.add_argument('--depth_multiplier',
                                 help='Assign a long paralog as the "main" sequence if it has a coverage depth '
                                      '<depth_multiplier> times all other long paralogs. Set to zero to not use '
                                      'depth. Default is %(default)s',
                                 default=10,
                                 type=int)
    parser_assemble.add_argument('--prefix',
                                 help='Directory name for pipeline output, default is to use the FASTQ file name.',
                                 default=None)
    parser_assemble.add_argument('--timeout_assemble',
                                 help='Kill long-running processes if they take longer than X percent of average.',
                                 default=0,
                                 type=int)
    parser_assemble.add_argument('--timeout_exonerate_contigs',
                                 help='Kill long-running processes if they take longer than X seconds. Default is %('
                                      'default)s',
                                 default=120,
                                 type=int)
    parser_assemble.add_argument('--target',
                                 help='Use the target file sequence with this taxon name in Exonerate searches for '
                                      'each gene. Other targets for that gene will be used only for read sorting. Can '
                                      'be a tab-delimited file (one <gene>\\t<taxon_name> per line) or a single taxon '
                                      'name.',
                                 default=None)
    parser_assemble.add_argument('--exclude',
                                 help='Do not use any sequence with the specified taxon name string in Exonerate '
                                      'searches. Sequenced from this taxon will still be used for read sorting.',
                                 default=None)
    parser_assemble.add_argument('--unpaired',
                                 help='Include a single FASTQ file with unpaired reads along with two paired read '
                                      'files',
                                 default=False)
    parser_assemble.add_argument('--no_stitched_contig', dest='no_stitched_contig', action='store_true',
                                 help='Do not create any stitched contigs. The longest single Exonerate hit will be '
                                      'used.',
                                 default=False)
    parser_assemble.add_argument('--bbmap_memory',
                                 default=1,
                                 type=int,
                                 help='GB memory (RAM ) to use for bbmap.sh with exonerate_hits.py. Default is %('
                                      'default)s.')
    parser_assemble.add_argument('--bbmap_subfilter',
                                 default=7,
                                 type=int,
                                 help='Ban alignments with more than this many substitutions. Default is %(default)s.')
    parser_assemble.add_argument('--bbmap_threads',
                                 default=2,
                                 type=int,
                                 help='Number of threads to use for BBmap when searching for chimeric stitched contig. '
                                      'Default is %(default)s.')
    parser_assemble.add_argument('--chimeric_stitched_contig_edit_distance',
                                 help='Minimum number of differences between one read of a read pair vs the '
                                      'stitched contig reference for a read pair to be flagged as discordant.',
                                 default=5,
                                 type=int)
    parser_assemble.add_argument('--chimeric_stitched_contig_discordant_reads_cutoff',
                                 help='Minimum number of discordant reads pairs required to flag a stitched contig as '
                                      'a potential chimera of contigs from multiple paralogs',
                                 default=5,
                                 type=int)
    parser_assemble.add_argument('--merged',
                                 help='For assembly with both merged and unmerged (interleaved) reads.',
                                 action='store_true',
                                 default=False)
    parser_assemble.add_argument('--run_intronerate',
                                 help='Run intronerate to recover fasta files for supercontigs with introns (if '
                                      'present), and introns-only.',
                                 action='store_true',
                                 dest='intronerate',
                                 default=False)
    parser_assemble.add_argument('--keep_intermediate_files',
                                 help='Keep all intermediate files and logs, which can be useful for '
                                      'debugging. Default action is to delete them, which greatly reduces the total '
                                      'file number).',
                                 action='store_true',
                                 dest='keep_intermediate_files',
                                 default=False)
    parser_assemble.add_argument('--no_padding_supercontigs',
                                 help='If Intronerate is run, and a supercontig is created by concatenating multiple '
                                      'SPAdes contigs, do not add 10 "N" characters between contig joins. By default, '
                                      'Ns will be added.',
                                 action='store_true',
                                 dest='no_padding_supercontigs',
                                 default=False)
    parser_assemble.add_argument('--verbose_logging',
                                 help='If supplied, enable verbose login. NOTE: this can increase the size of the log '
                                      'files by an order of magnitude.',
                                 action='store_true',
                                 dest='verbose_logging',
                                 default=False)

    # Set defaults for subparser <parser_assemble>:
    parser_assemble.set_defaults(blast=True)

    return parser_assemble


def add_stats_parser(subparsers):
    """
    Parser for hybpiper_stats, which now includes running get_seq_lengths

    :param argparse._SubParsersAction subparsers:
    :return None: no return value specified; default is None
    """

    parser_stats = subparsers.add_parser('stats', help='Gather statistics about the HybPiper run(s)')
    group_1 = parser_stats.add_mutually_exclusive_group(required=True)
    group_1.add_argument('--targetfile_dna', '-t_dna',
                         dest='targetfile_dna',
                         default=False,
                         help='FASTA file containing DNA target sequences for each gene. If there are multiple '
                              'targets for a gene, the id must be of the form: >Taxon-geneName')
    group_1.add_argument('--targetfile_aa', '-t_aa',
                         dest='targetfile_aa',
                         default=False,
                         help='FASTA file containing amino-acid target sequences for each gene. If there are multiple '
                              'targets for a gene, the id must be of the form: >Taxon-geneName')
    parser_stats.add_argument("sequence_type",
                              help="Sequence type (gene or supercontig) to recover lengths for",
                              choices=["gene", "GENE", "supercontig", "SUPERCONTIG"])
    parser_stats.add_argument('namelist',
                              help="Text file with names of HybPiper output directories, one per line")
    parser_stats.add_argument("--seq_lengths_filename",
                              help="File name for the sequence lengths *.tsv file. Default is <seq_lengths.tsv>.",
                              default='seq_lengths')
    parser_stats.add_argument("--stats_filename",
                              help="File name for the stats *.tsv file. Default is= <hybpiper_stats.tsv>",
                              default='hybpiper_stats')

    return parser_stats


def add_retrieve_sequences_parser(subparsers):
    """
    Parser for retrieve_sequences

    :param argparse._SubParsersAction subparsers: subparsers object to add parser(s) to
    :return None: no return value specified; default is None
    """

    parser_retrieve_sequences = subparsers.add_parser('retrieve_sequences',
                                                      help='Retrieve sequences generated from multiple runs of '
                                                           'HybPiper')
    group_1 = parser_retrieve_sequences.add_mutually_exclusive_group(required=True)
    group_1.add_argument('--targetfile_dna', '-t_dna',
                         dest='targetfile_dna',
                         default=False,
                         help='FASTA file containing DNA target sequences for each gene. If there are multiple '
                              'targets for a gene, the id must be of the form: >Taxon-geneName')
    group_1.add_argument('--targetfile_aa', '-t_aa',
                         dest='targetfile_aa',
                         default=False,
                         help='FASTA file containing amino-acid target sequences for each gene. If there are multiple '
                              'targets for a gene, the id must be of the form: >Taxon-geneName')
    parser_retrieve_sequences.add_argument('--sample_names',
                                           help='Directory containing Hybpiper output OR a file containing HybPiper '
                                                'output names, one per line',
                                           default=None)
    parser_retrieve_sequences.add_argument('--single_sample_name',
                                           help='A single sample name to recover sequences for',
                                           default=None)
    parser_retrieve_sequences.add_argument('sequence_type',
                                           help='Type of sequence to extract',
                                           choices=["dna", "aa", "intron", "supercontig"])
    parser_retrieve_sequences.add_argument("--hybpiper_dir",
                                           help='Specify directory containing HybPiper output')
    parser_retrieve_sequences.add_argument("--fasta_dir",
                                           help='Specify directory for output FASTA files')
    parser_retrieve_sequences.add_argument('--skip_chimeric_genes',
                                           action='store_true',
                                           dest='skip_chimeric',
                                           help='Do not recover sequences for putative chimeric genes',
                                           default=False)
    parser_retrieve_sequences.add_argument('--stats_file',
                                           help='Stats file produced by "hybpiper stats", required for selective '
                                                'filtering of retrieved sequences')
    parser_retrieve_sequences.add_argument('--filter_by', action='append', nargs=3,
                                           help='Provide three space-separated arguments: 1) column of the stats_file '
                                                'to filter by, 2) greater or less than symbol (> or <), '
                                                '3) a threshold - either an integer (raw number of genes) or float ('
                                                'percentage of genes in analysis).')

    return parser_retrieve_sequences


def add_paralog_retriever_parser(subparsers):
    """
    Parser for paralog_retriever

    :param argparse._SubParsersAction subparsers: subparsers object to add parser(s) to
    :return None: no return value specified; default is None
    """

    parser_paralog_retriever = subparsers.add_parser('paralog_retriever', help='Retrieve paralog sequences for a '
                                                                               'given gene, for all samples')
    parser_paralog_retriever.add_argument('namelist',
                                          help='Text file containing list of HybPiper output directories, '
                                               'one per line.')
    group_1 = parser_paralog_retriever.add_mutually_exclusive_group(required=True)
    group_1.add_argument('--targetfile_dna', '-t_dna',
                         dest='targetfile_dna',
                         default=False,
                         help='FASTA file containing DNA target sequences for each gene. Used to extract unique gene '
                              'names for paralog recovery. If there are multiple targets for a gene, the id must be '
                              'of the form: >Taxon-geneName')
    group_1.add_argument('--targetfile_aa', '-t_aa',
                         dest='targetfile_aa',
                         default=False,
                         help='FASTA file containing amino-acid target sequences for each gene. Used to extract '
                              'unique gene names for paralog recovery. If there are multiple targets for a gene, '
                              'the id must be of the form: >Taxon-geneName')
    parser_paralog_retriever.add_argument('--fasta_dir_all',
                                          help='Specify directory for output FASTA files (ALL). Default is '
                                               '"paralogs_all".',
                                          default='paralogs_all')
    parser_paralog_retriever.add_argument('--fasta_dir_no_chimeras',
                                          help='Specify directory for output FASTA files (no putative chimeric '
                                               'sequences). Default is "paralogs_no_chimeras".',
                                          default='paralogs_no_chimeras')
    parser_paralog_retriever.add_argument('--paralog_report_filename',
                                          help='Specify the filename for the paralog *.tsv report table',
                                          default='paralog_report')
    parser_paralog_retriever.add_argument('--paralogs_above_threshold_report_filename',
                                          help='Specify the filename for the *.txt list of genes with paralogs in '
                                               '<paralogs_list_threshold_percentage> number of samples',
                                          default='paralogs_above_threshold_report')
    parser_paralog_retriever.add_argument('--paralogs_list_threshold_percentage',
                                          help='Percent of total number of samples and genes that must have paralog '
                                               'warnings to be reported in the <genes_with_paralogs.txt> report file. '
                                               'The default is 0.0, meaning that all genes and samples with at least '
                                               'one paralog warning will be reported',
                                          type=float,
                                          default=0.0)
    parser_paralog_retriever.add_argument('--heatmap_filename',
                                          help='Filename for the output heatmap, saved by default as a *.png file. '
                                               'Defaults to "paralog_heatmap"',
                                          default='paralog_heatmap')
    parser_paralog_retriever.add_argument('--figure_length',
                                          type=int,
                                          help='Length dimension (in inches) for the output heatmap file. Default is '
                                               'automatically calculated based on the number of genes',
                                          default=None)
    parser_paralog_retriever.add_argument('--figure_height',
                                          type=int,
                                          help='Height dimension (in inches) for the output heatmap file. Default is '
                                               'automatically calculated based on the number of samples',
                                          default=None)
    parser_paralog_retriever.add_argument('--sample_text_size',
                                          type=int,
                                          help='Size (in points) for the sample text labels in the output heatmap '
                                               'file. Default is automatically calculated based on the number of '
                                               'samples',
                                          default=None)
    parser_paralog_retriever.add_argument('--gene_text_size',
                                          type=int,
                                          help='Size (in points) for the gene text labels in the output heatmap file. '
                                               'Default is automatically calculated based on the number of genes',
                                          default=None)
    parser_paralog_retriever.add_argument('--heatmap_filetype',
                                          choices=['png', 'pdf', 'eps', 'tiff', 'svg'],
                                          help='File type to save the output heatmap image as. Default is png',
                                          default='png')
    parser_paralog_retriever.add_argument('--heatmap_dpi',
                                          type=int,
                                          help='Dots per inch (DPI) for the output heatmap image. Default is 300',
                                          default='300')

    return parser_paralog_retriever


def add_gene_recovery_heatmap_parser(subparsers):
    """
    Parser for gene_recovery_heatmap

    :param argparse._SubParsersAction subparsers: subparsers object to add parser(s) to
    :return None: no return value specified; default is None
    """

    parser_gene_recovery_heatmap = subparsers.add_parser('recovery_heatmap',
                                                         help='Create a gene recovery heatmap for the HybPiper run')
    parser_gene_recovery_heatmap.add_argument('seq_lengths_file',
                                              help="Filename for the seq_lengths file (output of the 'hybpiper "
                                                   "stats' command)")
    parser_gene_recovery_heatmap.add_argument('--heatmap_filename',
                                              help='Filename for the output heatmap, saved by default as a *.png file. '
                                                   'Defaults to "heatmap"',
                                              default='heatmap')
    parser_gene_recovery_heatmap.add_argument('--figure_length',
                                              type=int,
                                              help='Length dimension (in inches) for the output heatmap file. '
                                                   'Default is automatically calculated based on the number of '
                                                   'genes',
                                              default=None)
    parser_gene_recovery_heatmap.add_argument('--figure_height',
                                              type=int,
                                              help='Height dimension (in inches) for the output heatmap file. '
                                                   'Default is automatically calculated based on the number of '
                                                   'samples',
                                              default=None)
    parser_gene_recovery_heatmap.add_argument('--sample_text_size',
                                              type=int,
                                              help='Size (in points) for the sample text labels in the output heatmap '
                                                   'file. Default is automatically calculated based on the '
                                                   'number of samples',
                                              default=None)
    parser_gene_recovery_heatmap.add_argument('--gene_text_size',
                                              type=int,
                                              help='Size (in points) for the gene text labels in the output heatmap '
                                                   'file. Default is automatically calculated based on the '
                                                   'number of genes',
                                              default=None)
    parser_gene_recovery_heatmap.add_argument('--heatmap_filetype',
                                              choices=['png', 'pdf', 'eps', 'tiff', 'svg'],
                                              help='File type to save the output heatmap image as. Default is *.png',
                                              default='png')
    parser_gene_recovery_heatmap.add_argument('--heatmap_dpi',
                                              type=int,
                                              help='Dot per inch (DPI) for the output heatmap image. Default is 300',
                                              default='300')

    return parser_gene_recovery_heatmap


def add_check_dependencies_parser(subparsers):
    """
    Parser for check_dependencies

    :param argparse._SubParsersAction subparsers: subparsers object to add parser(s) to
    :return None: no return value specified; default is None
    """

    parser_check_dependencies = subparsers.add_parser('check_dependencies',
                                                      help='Run a check for all pipeline dependencies and exit')

    # Set defaults for subparser <check_dependencies>:
    parser_check_dependencies.set_defaults(logger=None)

    return parser_check_dependencies


def add_check_targetfile_parser(subparsers):
    """
    Parser for check_targetfile

    :param argparse._SubParsersAction subparsers: subparsers object to add parser(s) to
    :return None: no return value specified; default is None
    """

    parser_check_target_file = subparsers.add_parser('check_targetfile',
                                                     help='Check the target file for sequences with low-complexity '
                                                          'regions, then exit')
    group_1 = parser_check_target_file.add_mutually_exclusive_group(required=True)
    group_1.add_argument('--targetfile_dna', '-t_dna',
                         dest='targetfile_dna',
                         default=False,
                         help='FASTA file containing DNA target sequences for each gene. If there are multiple '
                              'targets for a gene, the id must be of the form: >Taxon-geneName')
    group_1.add_argument('--targetfile_aa', '-t_aa',
                         dest='targetfile_aa',
                         default=False,
                         help='FASTA file containing amino-acid target sequences for each gene. If there are multiple '
                              'targets for a gene, the id must be of the form: >Taxon-geneName')
    parser_check_target_file.add_argument('--sliding_window_size',
                                          type=int,
                                          default=None,
                                          help='Number of characters (single-letter DNA or amino-acid codes) to '
                                               'include in the sliding window for low-complexity check')
    parser_check_target_file.add_argument('--complexity_minimum_threshold',
                                          type=float,
                                          default=None,
                                          help='Minimum threshold value. Beneath this value, the sequence in the '
                                               'sliding window is flagged as low-complexity, and the corresponding '
                                               'target file sequence is reported as having low-complexity regions ')

    # Set defaults for subparser <check_target_file>:
    parser_check_target_file.set_defaults(logger=None)

    return parser_check_target_file

