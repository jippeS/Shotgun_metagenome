import gzip
import os
import argparse
import sys
import pandas as pd


def get_output_file(input_file):
    file_name = ""
    input_list = input_file.split("/")[1:]

    for i in input_list:
        if i != input_list[-1]:
            file_name += f"/{i}"
    else:
        file_name += f"/{input_list[7]}_{input_list[-1]}"
    return file_name


class AnnotatePanres:
    def __init__(self, input_file, annotation_file):
        self.input_file = input_file
        self.input_file = "/export/projects/Transcriptomic/DNA/ELLY/proj1/samples/ELLY01/output/classified/panres/panres_31_07_2025_output.mapstat"
        self.annotation_file = annotation_file
        self.output_file = get_output_file(self.input_file)

    def open_input_file(self):
        # Read the first 6 lines
        with open(self.input_file, 'r') as f:
            first_six_lines = ''.join([next(f) for _ in range(6)])  # store as string
            # Load the rest of the file (starting from the 7th line) into pandas
            df = pd.read_csv(f, delimiter='\t')
        return df

    def open_annotate_file(self):
        with open(self.annotation_file, 'r') as f:
            df = pd.read_csv(f, delimiter='\t')

        #list of all databases to make into columns.
        unique_values = df[df['variable'] == 'database']['value'].unique()

        return df, unique_values


def argparser():
    parser = argparse.ArgumentParser(description="Filter FASTQ sequences with too many consecutive 'G's.")
    parser.add_argument("input_file", type=str, help="Absolute path to the mapstat file")
    parser.add_argument("annotation", type=str, help="Absolute path to the mapstat file")
    args = parser.parse_args()
    return args


def main():
    """
    Executes the program functions in order.
    """
    args = argparser()
    annotate = AnnotatePanres(input_file=args.input_file, annotation_file=args.annotation)
    annotate.open_input_file()
    annotate.open_annotate_file()

if __name__ == '__main__':
    sys.exit(main())
