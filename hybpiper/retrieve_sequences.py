#!/usr/bin/env python

"""
This script will get the sequences generated from multiple runs of the 'hybpiper assemble' command.
Specify either a directory with all the HybPiper output directories or a file containing sample names of interest.
It retrieves all the gene names from the target file used in the run of the pipeline.

You must specify whether you want the protein (aa), nucleotide (dna) sequences.

You can also specify 'intron' to retrieve the intron sequences, or 'supercontig' to get intron and exon sequences.

Will output unaligned fasta files, one per gene.
"""

import os
import sys
import argparse
from Bio import SeqIO
import logging


# Create a custom logger

# Log to Terminal (stderr):
console_handler = logging.StreamHandler(sys.stderr)
console_handler.setLevel(logging.INFO)

# Setup logger:
logger = logging.getLogger(f'hybpiper.{__name__}')

# Add handlers to the logger
logger.addHandler(console_handler)
logger.setLevel(logging.DEBUG)  # Default level is 'WARNING'


def get_chimeric_genes_for_sample(sample_directory_name):
    """
    Returns a list of putative chimeric gene sequences for a given sample

    :param str sample_directory_name: directory name for the sample
    :return list chimeric_genes_to_skip: a list of putative chimeric gene sequences for the sample
    """

    chimeric_genes_to_skip = []
    try:
        with open(f'{sample_directory_name}/'
                  f'{sample_directory_name}_genes_derived_from_putative_chimeric_stitched_contig.csv') as chimeric:
            lines = chimeric.readlines()
            for line in lines:
                chimeric_genes_to_skip.append(line.split(',')[1])
    except FileNotFoundError:  # This file should be written in assemble.py even if it's empty
        logger.info(f'No chimeric stitched contig summary file found for gene sample {sample_directory_name}!')
        raise

    return chimeric_genes_to_skip


def recover_sequences_from_all_samples(seq_dir, filename, target_genes, sample_names, hybpiper_dir=None,
                                       fasta_dir=None, skip_chimeric=False):
    """
    Recovers sequences (dna, amino acid, supercontig or intron) for all genes from all samples

    :param str seq_dir: directory to recover sequence from
    :param str filename: file name component used to reconstruct path to file
    :param list target_genes: list of unique gene names in the target file
    :param str sample_names: directory of samples, or text file with list of sample names
    :param None or str hybpiper_dir: if provided, a path to the directory containing HybPiper output
    :param None or str fasta_dir: directory name for output files, default is current directory
    :param bool skip_chimeric: if True, skip putative chimeric genes
    :return None:
    """

    if os.path.isdir(sample_names):
        sampledir = sample_names
        sample_names = [x for x in os.listdir(sampledir) if os.path.isdir(os.path.join(sampledir, x)) and not
                        x.startswith('.')]
    else:
        sample_names = [x.rstrip() for x in open(sample_names)]
        if hybpiper_dir:
            sampledir = hybpiper_dir
        else:
            sampledir = '.'

    if fasta_dir:
        fasta_dir = fasta_dir
        if not os.path.isdir(fasta_dir):
            os.mkdir(fasta_dir)
    else:
        fasta_dir = '.'

    logger.info(f'Retrieving {len(target_genes)} genes from {len(sample_names)} samples')
    for gene in target_genes:
        numSeqs = 0

        # Construct names for intron and supercontig output files:
        if seq_dir in ['intron', 'supercontig']:
            outfilename = f'{gene}_{filename}.fasta'
        else:
            outfilename = f'{gene}.{seq_dir}'

        with open(os.path.join(fasta_dir, outfilename), 'w') as outfile:
            for sample in sample_names:

                # Recover a list of putative chimeric genes for the sample, and skip gene if in list:
                if skip_chimeric:
                    chimeric_genes_to_skip = get_chimeric_genes_for_sample(sample)
                    # print(f'chimeric_genes_to_skip is: {chimeric_genes_to_skip}')
                    if gene in chimeric_genes_to_skip:
                        logger.info(f'Skipping putative chimeric stitched contig sequence for {gene}, sample {sample}')
                        continue

                # Get path to the gene/intron/supercontig sequence:
                if seq_dir == 'intron':
                    sample_path = os.path.join(sampledir, sample, gene, sample, 'sequences', seq_dir,
                                               f'{gene}_{filename}.fasta')
                else:
                    sample_path = os.path.join(sampledir, sample, gene, sample, 'sequences', seq_dir,
                                               f'{gene}.{seq_dir}')
                try:
                    seq = next(SeqIO.parse(sample_path, 'fasta'))
                    # print(seq)
                    SeqIO.write(seq, outfile, 'fasta')
                    numSeqs += 1
                # except FileNotFoundError or StopIteration:  # BioPython 1.80 returns StopIteration error?
                except FileNotFoundError:
                    pass
        logger.info(f'Found {numSeqs} sequences for gene {gene}.')


def recover_sequences_from_one_sample(seq_dir,
                                      filename,
                                      target_genes,
                                      single_sample_name,
                                      fasta_dir=None,
                                      skip_chimeric=False):
    """
    Recovers sequences (dna, amino acid, supercontig or intron) for all genes from one sample

    :param str seq_dir: directory to recover sequence from
    :param str filename: file name component used to reconstruct path to file (None, intron or supercontig)
    :param list target_genes: list of unique gene names in the target file
    :param str single_sample_name: directory of a single sample
    :param None or str fasta_dir: directory name for output files, default is current directory
    :param bool skip_chimeric: if True, skip putative chimeric genes
    :return None:
    """

    if not os.path.isdir(single_sample_name):
        sys.exit(f'Can not find a directory for sample {single_sample_name}, exiting...')

    # Create a user-supplied directory if provided, or write to the current directory if not:
    if fasta_dir:
        fasta_dir = fasta_dir
        if not os.path.isdir(fasta_dir):
            os.mkdir(fasta_dir)
    else:
        fasta_dir = '.'
    logger.info(f'Retrieving {len(target_genes)} genes from sample {single_sample_name}...')
    sequences_to_write = []

    # Construct names for intron and supercontig output files:
    if seq_dir in ['intron', 'supercontig']:
        outfilename = f'{filename}.fasta'
    else:
        outfilename = f'{seq_dir}.fasta'
    for gene in target_genes:
        numSeqs = 0
        # Recover a list of putative chimeric genes for the sample, and skip gene if in list:
        if skip_chimeric:
            chimeric_genes_to_skip = get_chimeric_genes_for_sample(single_sample_name)
            # print(f'chimeric_genes_to_skip is: {chimeric_genes_to_skip}')
            if gene in chimeric_genes_to_skip:
                logger.info(f'Skipping putative chimeric stitched contig sequence for {gene}, sample'
                            f' {single_sample_name}')
                continue

        # Get path to the gene/intron/supercontig sequence:
        if seq_dir == 'intron':
            sample_path = os.path.join(single_sample_name, gene, single_sample_name, 'sequences', seq_dir,
                                       f'{gene}_{filename}.fasta')
        else:
            sample_path = os.path.join(single_sample_name, gene, single_sample_name, 'sequences', seq_dir,
                                       f'{gene}.{seq_dir}')
        try:
            seq = next(SeqIO.parse(sample_path, 'fasta'))
            seq.id = f'{seq.id}-{gene}'
            sequences_to_write.append(seq)
            numSeqs += 1
        # except FileNotFoundError or StopIteration:  # BioPython 1.80 returns StopIteration error?
        except FileNotFoundError:
            pass
        logger.info(f'Found {numSeqs} sequences for gene {gene}.')

    with open(os.path.join(fasta_dir, outfilename), 'w') as outfile:
        SeqIO.write(sequences_to_write, outfile, 'fasta')


def standalone():
    """
    Used when this module is run as a stand-alone script. Parses command line arguments and runs function main().:
    """

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    group_1 = parser.add_mutually_exclusive_group(required=True)
    group_1.add_argument('--targetfile_dna', '-t_dna', dest='targetfile_dna',
                         help='FASTA file containing DNA target sequences for each gene. If there are multiple '
                              'targets for a gene, the id must be of the form: >Taxon-geneName')
    group_1.add_argument('--targetfile_aa', '-t_aa', dest='targetfile_aa',
                         help='FASTA file containing amino-acid target sequences for each gene. If there are multiple '
                              'targets for a gene, the id must be of the form: >Taxon-geneName')
    parser.add_argument('--sample_names',
                        help='Directory containing Hybpiper output OR a file containing HybPiper output names, '
                             'one per line',
                        default=None)
    parser.add_argument('--single_sample_name',
                        help='A single sample name to recover sequences for', default=None)
    parser.add_argument('sequence_type',
                        help='Type of sequence to extract',
                        choices=['dna', 'aa', 'intron', 'supercontig'])
    parser.add_argument('--hybpiper_dir', help='Specify directory containing HybPiper output',
                        default=None)
    parser.add_argument('--fasta_dir', help='Specify directory for output FASTA files',
                        default=None)
    parser.add_argument('--skip_chimeric_genes',
                        action='store_true',
                        dest='skip_chimeric',
                        help='Do not recover sequences for putative chimeric genes',
                        default=False)

    parser.set_defaults(targetfile_dna=False, targetfile_aa=False)

    args = parser.parse_args()
    main(args)


def main(args):
    """
    Entry point for the assemble.py module.

    :param argparse.Namespace args:
    """

    # Set target file name:
    if args.targetfile_dna:
        targetfile = args.targetfile_dna
    elif args.targetfile_aa:
        targetfile = args.targetfile_aa

    # Check some command line parameters:
    if not args.sample_names and not args.single_sample_name:
        sys.exit(f'Please supply either the --sample_names or --single_sample_name flag and corresponding arguments!')
    if args.sample_names and args.single_sample_name:
        sys.exit(f'Please supply either the --sample_names or --single_sample_name flag and corresponding arguments!')

    # Set sequence directory name and file names:
    if args.sequence_type == 'dna':
        seq_dir = "FNA"
        filename = None
    elif args.sequence_type == 'aa':
        seq_dir = "FAA"
        filename = None
    elif args.sequence_type == 'intron':
        seq_dir = 'intron'
        filename = 'introns'
    elif args.sequence_type == 'supercontig':
        seq_dir = 'intron'
        filename = 'supercontig'

    # Use gene names parsed from a target file.
    target_genes = list(set([x.id.split('-')[-1] for x in SeqIO.parse(targetfile, 'fasta')]))

    # Recover sequences from all samples:
    if args.sample_names:
        recover_sequences_from_all_samples(seq_dir,
                                           filename,
                                           target_genes,
                                           args.sample_names,
                                           args.hybpiper_dir,
                                           args.fasta_dir,
                                           args.skip_chimeric)
    elif args.single_sample_name:
        recover_sequences_from_one_sample(seq_dir,
                                          filename,
                                          target_genes,
                                          args.single_sample_name,
                                          args.fasta_dir,
                                          args.skip_chimeric)


if __name__ == "__main__":
    standalone()
