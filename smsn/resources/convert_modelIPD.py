import pandas as pd
import sys
import os

import argparse

""" Converts IPDMODEL files into a lightweight hd5 format """


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-i","--input",
                            help='Path of the csv produced by ipdtools',
                            required=True)

    parser.add_argument("-o","--output",
                            help="Path of the output h5 file ",
                            required=True
                        )

    args = parser.parse_args()
    parser.input = os.path.realpath(args.input)
    parser.output = os.path.realpath(args.output)

    df = pd.read_csv(parser.input,sep=";")

    scaffolds = {}

    for scaffold in df["Fasta_ID"].unique():
        tmp_df = df[df["Fasta_ID"] == scaffold].copy()
        tmp_df.to_hdf(parser.output,key=scaffold)
