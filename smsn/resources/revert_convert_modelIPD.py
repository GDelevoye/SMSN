import pandas as pd
import sys
import os

import argparse

""" Converts IPDMODEL files into a lightweight hd5 format """


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-i","--input",
                            help='Path of the HDF file of ipdtools output',
                            required=True)

    parser.add_argument("-o","--output",
                            help="Path of the output .csv file ",
                            required=True
                        )

    args = parser.parse_args()
    parser.input = os.path.realpath(args.input)
    parser.output = os.path.realpath(args.output)

    with pd.HDFStore(parser.input) as tmp:
        scaffolds = list(tmp)
        # print(scaffolds)


    proto_df = []

    for scaffold in scaffolds:
        tmp_df = pd.read_hdf(parser.input,key =scaffold)
        proto_df.append(tmp_df)

    df = pd.concat(proto_df)
    df.to_csv(parser.output,sep=";",index=False)
