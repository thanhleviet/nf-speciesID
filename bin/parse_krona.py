#!/usr/bin/env python3

import pandas as pd
from argparse import ArgumentParser
from StringIO import StringIO


def parse_args():
    parser = ArgumentParser(description="Parsing krona tab output")
    parser.add_argument("-i", "--krona", help="Input file", required=True)
    parser.add_argument("-t", "--taxa", help="Taxa file",
                        default="taxonomy.tab")
    parser.add_argument("--thres", help="Thresold", default=0.15)
    return parser.parse_args()


def main():
    args = parse_args()
    krona = args.krona
    taxa_db = args.taxa
    thres = args.thres
    # Read taxonomy.tab and load only 2 columns taxID and name
    taxa = pd.read_csv(taxa_db, sep="\t", header=0, index_col=[
                       0], usecols=[0, 4], names=["taxID", "name"])
    # Read tab file output by krona
    kr_tsv = pd.read_csv(krona, sep="\t")
    # Count reads by taxaID value
    un = kr_tsv['taxID'].value_counts().to_frame()
    un.index.name = "taxID"
    un = un.rename(columns={"taxID": "counts"})
    # Join taxonomy database with krona
    df = taxa.join(un)
    df = df[df.counts > 0][['name', 'counts']]
    # Estimate percentage for each samples
    df['reads_percentage'] = df.counts / sum(df.counts)
    df = df.sort_values('reads_percentage', ascending=False)  # Sort descending
    # Filter species >= a given threshold
    df = df[df.reads_percentage >= thres]
    if df.shape[0] == 1:
        status = "PASS"
        if df.iloc[0, 2] <= 0.45:
            status = "LIKELY_CONTAMINATED"
        print("{}\t[('{}',{})]".format(status, df.iloc[0, 0], df.iloc[0, 2]))
    else:
        status = "CONTAMINATED"
        species_list = zip(list(df.name), list(
            df.reads_percentage.apply(lambda x: round(x, 2))))
        print("{}\t{}".format(status, species_list))


if __name__ == "__main__":
    main()
