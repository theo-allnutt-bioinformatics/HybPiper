#!/usr/bin/env python

"""
Takes the seq_lengths.txt file (output from running 'hybpiper get_seq_lengths') as input.

For each sample and for each gene, calculates the percentage length recovered. This percentage is calculated as a
fraction of the mean length for representative gene sequences in the target/bait file provided.

Generates a heatmap of percentage length recovery for each sample and each gene.

"""

import logging
import sys
import argparse
import os
import datetime

# Import non-standard-library modules:

try:
    import pandas as pd
except ImportError:
    sys.exit(f"Required Python package 'pandas' not found. Is it installed for the Python used to run this script?")

try:
    import seaborn as sns
except ImportError:
    sys.exit(f"Required Python package 'seaborn' not found. Is it installed for the Python used to run this script?")

try:
    import matplotlib.pyplot as plt
except ImportError:
    sys.exit(f"Required Python package 'matplotlib' not found. Is it installed for the Python used to run this script?")

########################################################################################################################
########################################################################################################################
# Get current working directory and host name:

cwd = os.getcwd()


# Configure logger:
def setup_logger(name, log_file, console_level=logging.INFO, file_level=logging.DEBUG,
                 logger_object_level=logging.DEBUG):
    """
    Function to create a logger instance.

    By default, logs level DEBUG and above to file.
    By default, logs level INFO and above to stderr and file.

    :param string name: name for the logger instance
    :param string log_file: filename for log file
    :param string console_level: level for logging to console
    :param string file_level: level for logging to file
    :param string logger_object_level: level for logger object
    :return: a logger object
    """

    # Get date and time string for log filename:
    date_and_time = datetime.datetime.now().strftime("%Y-%m-%d-%H_%M_%S")

    # Log to file:
    file_handler = logging.FileHandler(f'{log_file}_{date_and_time}.log', mode='w')
    file_handler.setLevel(file_level)
    file_format = logging.Formatter('%(asctime)s - %(filename)s - %(name)s - %(funcName)s - %(levelname)s - %('
                                    'message)s')
    file_handler.setFormatter(file_format)

    # Log to Terminal (stdout):
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(console_level)
    console_format = logging.Formatter('%(message)s')
    console_handler.setFormatter(console_format)

    # Setup logger:
    logger_object = logging.getLogger(name)
    logger_object.setLevel(logger_object_level)  # Default level is 'WARNING'

    # Add handlers to the logger
    logger_object.addHandler(console_handler)
    logger_object.addHandler(file_handler)

    return logger_object


# Create logger(s):
logger = setup_logger(__name__, 'gene_recovery_heatmap')


########################################################################################################################
########################################################################################################################
# Define functions:

def get_figure_dimensions(df):
    """
    Takes a dataframe and returns figure length (inches), figure height (inches), sample_text_size and gene_id_text_size
    values based on the number of samples and genes in the seq_lengths.txt input file provided.

    :param pandas.core.frame.DataFrame df: pandas dataframe of seq_lengths.txt after filtering and pivot
    :return float fig_length, figure_height, sample_text_size, gene_id_text_size:
    """

    num_samples = len(df.index)
    num_genes = len(df.columns)

    logger.info(f'Number of samples in input lengths file is: {num_samples}')
    logger.info(f'Number of genes in input lengths file is: {num_genes}')

    # Set some dimensions for a given number of samples:
    if num_samples <= 10:
        sample_text_size = 10
        figure_height = 10/2.54
    elif 10 < num_samples <= 20:
        sample_text_size = 8
        figure_height = 10/2.54
    elif 20 < num_samples <= 50:
        sample_text_size = 8
        figure_height = 15/2.54
    elif 50 < num_samples <= 100:
        sample_text_size = 6
        figure_height = 15/2.54
    elif 100 < num_samples <= 200:
        sample_text_size = 4
        figure_height = 21/2.54
    elif 200 < num_samples <= 400:
        sample_text_size = 3
        figure_height = 21/2.54
    elif num_samples > 400:
        sample_text_size = 3
        figure_height = 21/2.54

    # Set some dimensions for a given number of genes (i.e. number of unique genes in target file):
    if num_genes <= 10:
        gene_id_text_size = 10
        fig_length = 10.7/2.54
    elif 10 < num_genes <= 20:
        gene_id_text_size = 10
        fig_length = 15.7/2.54
    elif 20 < num_genes <= 50:
        gene_id_text_size = 8
        fig_length = 20.7/2.54
    elif 50 < num_genes <= 100:
        gene_id_text_size = 6
        fig_length = 20.7/2.54
    elif 100 < num_genes <= 200:
        gene_id_text_size = 4
        fig_length = 29.7/2.54
    elif 200 < num_genes <= 400:
        gene_id_text_size = 12
        fig_length = 254/2.54
    elif num_genes > 400:
        gene_id_text_size = 3
        fig_length = 29.7/2.54

    logger.info(f'fig_length: {fig_length}, figure_height: {figure_height}, sample_text_size: {sample_text_size}, '
                f'gene_id_text_size: {gene_id_text_size}')

    return fig_length, figure_height, sample_text_size, gene_id_text_size


########################################################################################################################
########################################################################################################################
# Run script:

def standalone():
    """
    Used when this module is run as a stand-alone script. Parses command line arguments and runs function main().:
    """

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-seq_lengths_file',
                        help="filename for the seq_lengths file (output of 'hybpiper get_seq_lengths'). If not "
                             "provided, the file 'seq_lengths.txt' will be searched for by default", default=None)

    args = parser.parse_args()
    main(args)


def main(args):
    """
    Entry point for the hybpiper.py module.

    :param argparse.Namespace args:
    """

    logger.info(f'Running {__name__} with: {args}')

    if args.seq_lengths_file and not os.path.exists(args.seq_lengths_file):
        logger.info(f'Can not find file "{args.seq_lengths_file}". Is it in the current working directory?')
    elif args.seq_lengths_file and os.path.exists(args.seq_lengths_file):
        sample_filename = args.seq_lengths_file
    else:
        sample_filename = "seq_lengths.txt"

    # Read in the sequence length file:
    df = pd.read_csv(sample_filename, delimiter='\t', )

    # For each sample, divide each gene length by the MeanLength value for that gene:
    df.loc[:, df.columns[1]:] = df.loc[:, df.columns[1]:].div(df.iloc[0][df.columns[1]:])

    # For each length ratio, if the value is greater than 1, assign it to 1:
    df.where(df.loc[:, df.columns[1]:] < 1, 1, inplace=True)

    # Drop the gene MeanLengths row:
    df.drop(labels=0, axis=0, inplace=True)

    # Melt the dataframe to make it suitable for the df.pivot method:
    df = df.melt(id_vars=['Species'], var_name='gene_id', value_name='percentage_recovery')

    # Change percentage values to numeric:
    df["percentage_recovery"] = df["percentage_recovery"].apply(pd.to_numeric)

    # Pivot the dataframe for input into the seaborn heatmap function:
    df = df.pivot(index='Species', columns='gene_id', values='percentage_recovery')

    # Get figure dimension and label text size based on number of samples and genes:
    fig_length, figure_height, sample_text_size, gene_id_text_size = get_figure_dimensions(df)

    # Create heatmap:
    sns.set(rc={'figure.figsize': (fig_length, figure_height)})
    sns.set_style('ticks')  # options are: white, dark, whitegrid, darkgrid, ticks
    cmap = 'bone_r'  # sets colour scheme
    heatmap = sns.heatmap(df, vmin=0, cmap=cmap)
    heatmap.tick_params(axis='x', labelsize=gene_id_text_size)
    heatmap.tick_params(axis='y', labelsize=sample_text_size)
    heatmap.set_yticklabels(heatmap.get_yticklabels(), rotation=0)
    heatmap.set_xlabel("Gene name", fontsize=14, fontweight='bold', labelpad=20)
    heatmap.set_ylabel("Sample name", fontsize=14, fontweight='bold', labelpad=20)
    plt.title("Percentage length recovery for each gene, relative to mean of baitfile references", fontsize=20, y=1.05)
    plt.tight_layout()

    # Save heatmap as png file:
    logger.info(f'Saving heatmap as file "heatmap.png"')
    plt.savefig('heatmap.png', dpi=300)


########################################################################################################################
########################################################################################################################

if __name__ == '__main__':
    standalone()

########################################################################################################################
########################################################################################################################