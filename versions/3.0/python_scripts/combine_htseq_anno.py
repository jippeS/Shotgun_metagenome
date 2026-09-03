import gzip
import os
import argparse
import sys
import pandas as pd


class AnnotatePanres:
    def __init__(self, args):
        self.input_htseq = args.input_htseq
        self.input_anno = args.input_anno
        self.input_ortho = args.input_ortho
        self.output_file = args.output

    def parse_htseq(self):
        df_htseq = pd.read_csv(self.input_htseq, sep="\t", skipfooter=5, header=None)
        df_htseq = df_htseq.rename(columns={
            df_htseq.columns[1]: "htseq-count"
        })
        return df_htseq

    def parse_annotate(self):
        df_anno = pd.read_csv(self.input_anno, sep="\t",skiprows=4, skipfooter=3)
        df_anno.columns = df_anno.columns.str.lstrip("#")
        return df_anno

    def parse_ortho(self):
        df_ortho = pd.read_csv(self.input_ortho, sep="\t",skiprows=5, skipfooter=3)
        df_ortho.columns = df_ortho.columns.str.lstrip("#")
        return df_ortho

    def combine(self, df_htseq, df_anno):
        df_combined = df_anno.merge(
            df_htseq,
            left_on="query",
            right_on=df_htseq.columns[0],
            how="left").drop(columns=df_htseq.columns[0])
        return df_combined

    def combined_final(self, df_combined, df_ortho):

        df_combined_final = df_combined.merge(
            df_ortho,
            left_on="query",
            right_on="qseqid",
            how="left"
        ).drop(columns=df_ortho.columns[0])
        return df_combined_final

    def write_output(self, df_combined):
        df_combined.to_csv(self.output_file, sep="\t", index=False)

def argparser():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--input_htseq", type=str, help="Absolute path to the mapstat file")
    parser.add_argument("--input_anno", type=str, help="Absolute path to the mapstat file")
    parser.add_argument("--input_ortho", type=str, help="Absolute path to the mapstat file")
    parser.add_argument("--output", type=str, help="output file")
    args = parser.parse_args()
    return args


def main():
    """
    Executes the program functions in order.
    """
    args = argparser()
    annotate = AnnotatePanres(args)
    df_htseq = annotate.parse_htseq()
    df_anno = annotate.parse_annotate()
    df_ortho = annotate.parse_ortho()
    df_combined = annotate.combine(df_htseq, df_anno)
    df_combined_final = annotate.combined_final(df_combined, df_ortho)
    annotate.write_output(df_combined_final)

if __name__ == '__main__':
    sys.exit(main())
