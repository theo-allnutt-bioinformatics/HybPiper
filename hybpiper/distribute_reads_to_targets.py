#!/usr/bin/env python

"""
This script is part of a pipeline to extract phylogenetically-useful sequences from
Illumina data using the targeted (liquid-phase) sequence enrichment approach.

After a BLASTx search of the raw reads against the target sequences, the reads need to be
sorted according to the successful hits. This script takes the BLASTx output (tabular)
and the raw read files, and distributes the reads into FASTA files ready for assembly.

If there are multiple BLAST results (for example, one for each read direction),
concatenate them prior to sorting. # CJJ not still true?
"""

import os
import errno
import argparse
import logging
from hybpiper.distribute_reads_to_targets_bwa import distribute_reads
import progressbar


# Create logger:
logger = logging.getLogger(f'hybpiper.assemble.{__name__}')


def mkdir_p(path):
    """
    Creates a directory corresponding the the given path, if it doesn't already exist.

    :param str path: path of directory to create
    :return:
    """
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def blast2dict(blastfilename):
    f1=open(blastfilename,'r')

    #print("blast2dict")
    
    data={}
    genes=[]
    
    for i in f1:
    
        k=i.split("\t")
        #if k.endswith("/1") or k.endswith("/1"):
            #k=k[:-2]
            
        gene=k[1].split("-")[-1]
        
        if gene not in genes:
        
            genes.append(gene)
            data[gene]=set()
            data[gene].add(k[0])
            #p0=sp.Popen("mkdir %s/%s" %(seqname,gene),shell=True).wait()
            
        else:
            
            data[gene].add(k[0])
            
    f1.close()
    
    return data

def read_sorting(blastfilename):

    read_hit_dict=blast2dict(blastfilename)

    return read_hit_dict


def write_paired_seqs(target, ID1, Seq1, Qual1, ID2, Seq2, Qual2, single=True, merged=False):
    """
    Writes interleaved fasta files, and also interleaved fastq files if merged=True, to the corresponding gene directory

    :param str target: gene name e.g. gene001
    :param str ID1: fasta/fastq header for R1
    :param str Seq1: fasta/fastq sequence for R1
    :param str Qual1: fastq quality scores for R1
    :param str ID2: fasta/fastq head for R2
    :param str Seq2: fasta/fastq sequence for R2
    :param str Qual2: fastq quality scores for R2
    :param bool single: # CJJ hardcoded as True - remove?
    :param bool merged: If True, write fastq seqs as well as fasta
    :return::
    """

    mkdir_p(target)
    if single:  # If True, write paired reads in interleaved format
        if merged:
            outfile = open(os.path.join(target, f'{target}_interleaved.fastq'), 'a')
            outfile.write(f'@{ID1}\n{Seq1}\n+\n{Qual1}\n')
            outfile.write(f'@{ID2}\n{Seq2}\n+\n{Qual2}\n')
            outfile.close()

        outfile = open(os.path.join(target, f'{target}_interleaved.fasta'), 'a')
        outfile.write(f'>{ID1}\n{Seq1}\n')
        outfile.write(f'>{ID2}\n{Seq2}\n')
        outfile.close()


def write_single_seqs(target, ID1, Seq1):
    """
    Writes a fasta file of single-end/unpaired reads to the corresponding gene directory

    :param str target: gene name e.g. gene001
    :param str ID1: fasta/fastq header for R1
    :param str Seq1: fasta/fastq sequence for R1
    :return:
    """

    mkdir_p(target)
    outfile = open(os.path.join(target, f'{target}_unpaired.fasta'), 'a')
    outfile.write(f'>{ID1}\n{Seq1}\n')
    outfile.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('blast_filename', help='Name of blast file from BLASTX mapping of reads')
    parser.add_argument('readfiles', nargs='+', help="List of readfiles")
    parser.add_argument("--merged", help="If provided, write fastq files for bbmerge", action="store_true",
                        default=False)
    args = parser.parse_args()

    logging.info(f'{"[NOTE]:":10} Running script distribute_reads_to_targets.py with {args}')
    readfiles = args.readfiles
    read_hit_dict = blast2dict(args.blast_filename)
    logging.info(f'{"[NOTE]":10} Unique reads with hits: {len(read_hit_dict)}')
    distribute_reads(readfiles, read_hit_dict, merged=args.merged)


if __name__ == '__main__':
    main()
