#!/usr/bin/env python3

import argparse
import sys
import pandas as pd


def main():
    # Parse input arguments.
    args = parse_args()
    df = pd.read_csv(args.transcript)

    # Filter transcripts. Ignore negative controls
    filtered_df = df[
        df["qv"].ge(args.min_qv) &
        df["x_location"].between(args.min_x, args.max_x) &
        df["x_location"].between(args.min_x, args.max_x) &
        ~df["feature_name"].str.startswith(
            ("NegControlProbe_", "antisense_", "NegControlCodeword_", "BLANK_")
        )
    ]
    # Change cell_id of cell-free transcripts from -1 to 0
    filtered_df["cell_id"].replace(-1, 0, inplace=True)

    # Output filtered transcripts to CSV
    filename = '_'.join([
        f'X{args.min_x}-{args.max_x}',
        f'Y{args.min_y}-{args.max_y}',
        'filtered_transcripts.csv'
    ])
    filtered_df.to_csv(filename, index=False, encoding='utf-8')


# --------------------------
# Helper functions


def parse_args():
    """
    Parses command-line options for main().
    """
    summary = (
        'Filter transcripts from transcripts.csv based on Q-Score threshold '
        'and upper bounds on x and y coordinates. Remove negative controls.'
    )
    parser = argparse.ArgumentParser(description=summary)
    # Required Args
    required_named = parser.add_argument_group('required named arguments')
    ts_msg = 'The path to the transcripts.csv file produced by Xenium.'
    required_named.add_argument('-transcript', required=True, help=ts_msg)
    # Optional Args
    qs_msg = 'The minimum Q-Score to pass filtering. (default: 20.0)'
    min_msg = (
        'Only keep transcripts whose x/y-coordinate is greater than specified '
        'limit. If no limit is specified, the minimum value will be 0.0'
    )
    max_msg = (
        'Only keep transcripts whose x/y-coordinate is less than specified '
        'limit. If no limit is specified, the default value will retain all '
        'transcripts since Xenium slide is <24000 microns in x and y.'
    )
    parser.add_argument('-min_qv', default='20.0', type=float, help=qs_msg)
    parser.add_argument('-min_x', default='0.0', type=float, help=min_msg)
    parser.add_argument('-max_x', default='24000.0', type=float, help=max_msg)
    parser.add_argument('-min_y', default='0.0', type=float, help=min_msg)
    parser.add_argument('-max_y', default='24000.0', type=float, help=max_msg)

    try:
        opts = parser.parse_args()
    except:
        sys.exit(0)

    return opts


if __name__ == "__main__":
    main()
