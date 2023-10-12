#!/usr/bin/env python3

import argparse
import csv
import os
import sys
import pandas as pd
import scipy.sparse as sparse
import scipy.io as sio


def main():
    '''
    Filter and convert Baysor transcripts Parquet file to sparse counts matrix
    '''
    # Parse input arguments.
    args = parse_args()

    # Check for existence of input file.
    if not os.path.exists(args.baysor):
        print(f'The specified Baysor output {args.baysor} does not exist!')
        sys.exit(0)

    # Check if output folder already exists.
    if os.path.exists(args.out):
        print(f'The specified output folder {args.out} already exists!')
        sys.exit(0)

    # Read 5 columns from transcripts Parquet file
    transcripts_df = pd.read_csv(
        args.baysor,
        usecols=['gene', 'cell', 'assignment_confidence']
    )
    # Remove transcripts below user-specified cutoff
    mask = transcripts_df['assignment_confidence'] >= args.conf_cutoff
    transcripts_df = transcripts_df[mask]

    matrix = transcripts_df.groupby(['gene', 'cell']).size().unstack()
    cells = transcripts_df['cell'].unique()[1:]
    features = transcripts_df['gene'].unique()

    # Call a helper function to create Seurat and Scanpy compatible MTX output
    write_sparse_mtx(args, matrix, cells, features)


# --------------------------
# Helper functions


def parse_args():
    """Parses command-line options for main()."""
    summary = (
        'Map Xenium transcripts to Baysor segmentation result. '
        'Generate Seurat/Scanpy-compatible feature-cell matrix.'
    )
    parser = argparse.ArgumentParser(description=summary)

    # Require args
    required_named = parser.add_argument_group('required named arguments')
    msg = 'The path to the *segmentation.csv file produced by Baysor.'
    required_named.add_argument('-baysor', required=True, help=msg)
    msg = 'The name of output folder in which feature-cell matrix is written.'
    required_named.add_argument('-out', required=True, help=msg)

    # Optional args
    msg = 'Ignore transcripts with assignment confidence below this threshold.'
    parser.add_argument('-conf_cutoff', default='0.9', type=float, help=msg)
    msg = (
        'Reporting interval. Will print message to stdout whenever the '
        'specified # of transcripts is processed.'
    )
    parser.add_argument('-rep_int', default='100000', type=int, help=msg)

    # Parse args or die
    try:
        opts = parser.parse_args()
    except:
        sys.exit(0)

    return opts


def write_sparse_mtx(args, matrix, cells, features):
    """Write feature-cell matrix in Seurat/Scanpy-compatible MTX format"""

    # Create the matrix folder.
    os.mkdir(args.out)

    # Convert matrix to scipy's COO sparse matrix.
    sparse_mat = sparse.coo_matrix(matrix.values)

    # Write matrix in MTX format.
    sio.mmwrite(args.out + "/matrix.mtx", sparse_mat)

    # Write cells as barcodes.tsv. File name is chosen to ensure
    # compatibility with Seurat/Scanpy.
    with open(args.out + "/barcodes.tsv", 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
        for cell in cells:
            writer.writerow(["cell_" + str(cell)])

    # Write features as features.tsv. Write 3 columns to ensure
    # compatibility with Seurat/Scanpy.
    with open(args.out + "/features.tsv", 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
        for feat in features:
            feat = str(feat)
            row = [feat]*2
            if feat.startswith(("NegControlProbe_", "antisense_")):
                row += ["Negative Control Probe"]
            elif feat.startswith("NegControlCodeword_"):
                row += ["Negative Control Codeword"]
            elif feat.startswith("BLANK_"):
                row += ["Blank Codeword"]
            else:
                row += ["Gene Expression"]
            writer.writerow(row)


if __name__ == "__main__":
    main()
