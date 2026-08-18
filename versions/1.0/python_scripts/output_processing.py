import gzip
import os
import argparse
import sys
import time

class OutputTransformer:
    def __init__(self, input_file, output_file, sample_name):
        self.input_file = input_file
        self.output_file = output_file
        self.sample_name = sample_name



    def filter_fastq(self):
        """
        Opens the gzip files and calls the filter function.
        """
        with open(self.input_file, "r") as infile, open(self.output_file, "w") as outfile:
            filter(infile, outfile)

    def filter(self, infile, outfile):
        """
        Writes the lines to a fastq file if it passes certain checks.
        """
        count = 0
        start = time.time()
        count_remove = 0
        if self.direction == "reverse":
            number = 2
        elif self.direction == "forward":
            number = 1
        else:
            sys.exit(print(f"Something went wrong with the naming: {self.direction}"))

        while True:
            header = self.reform_header(infile.readline(), number)
            sequence = infile.readline()
            plus = infile.readline()
            quality = infile.readline()

            if not header or not sequence or not plus or not quality:
                break
            else:
                if self.is_valid_sequence(sequence):
                    outfile.write(header)
                    outfile.write(sequence)
                    outfile.write(plus)
                    outfile.write(quality)
                    count += 1
                else:
                    removed_file.write(header)
                    removed_file.write(sequence)
                    removed_file.write(plus)
                    removed_file.write(quality)
                    count_remove += 1

            # if count == 100000:
        end = time.time()
        react_time = end - start
        print(f"Sequences removed: {count_remove}")
        print("\n\nIt took the program {} seconds to complete.".format(round(react_time)))

        return count

def argparser():
    parser = argparse.ArgumentParser(description="Filter FASTQ sequences with too many consecutive 'G's.")
    parser.add_argument("input_file", type=str, help="Path to the input tsv")
    parser.add_argument("output_file", type=str, help="Path to the input tsv")
    parser.add_argument("sample", type=str, help="Sample name")
    args = parser.parse_args()
    return args



def main():
    """
    Executes the program functions in order.
    """
    args = argparser()
    filterer = OutputTransformer(input_file=args.input_file, sample_name=args.sample)
    filterer.filter_fastq()

if __name__ == '__main__':
    sys.exit(main())