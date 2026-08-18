import gzip
import os
import argparse
import sys
import time

class FastqFilter:
    def __init__(self, input_file, poly_g, output):
        self.input_file = input_file
        # From an entire path, retrieves the file, and strips away the extension leaving only forward or reverse text.
        self.direction = self.input_file.split("/")[-1].split("_")[1].split(".")[0]
        self.output_dir = output
        self.output_file = self._generate_output_filename(input_file)
        self.min_consecutive_gs = "G" * poly_g

    def _generate_output_filename(self, input_file):
        """
        Generate the outputfile
        """
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        base, ext = os.path.splitext(input_file.split("/")[-1])
        if ext == '.gz':
            base, ext2 = os.path.splitext(base)
            output_file = f"{self.output_dir}{base}_filtered{ext2}{ext}"
            self.remove_seq = f"{self.output_dir}{base}_removed_seqs{ext2}{ext}"
        else:
            output_file = f"{self.output_dir}{base}_filtered.{ext}"
            self.remove_seq = f"{self.output_dir}{base}_removed_seqs.{ext}"
        return output_file

    def is_valid_sequence(self, sequence):
        """
        Check if there isnt a poly G part in the sequence.
        """
        return self.min_consecutive_gs not in sequence

    def filter_fastq(self):
        """
        Opens the gzip files and calls the filter function.
        """
        print(self.input_file)
        print(self.output_file)
        with gzip.open(self.input_file, 'rt') as infile, gzip.open(self.output_file, 'wt') as outfile, gzip.open(self.remove_seq, 'wt') as removed_file:
            self.filter(infile, outfile, removed_file)

    def reform_header(self, header, fw_or_rv):
        """
        Reforms the header to be compatible with the kneaddata paired function.
        """
        reformed_header = header.strip().split(" ")
        reformed_header = f"{reformed_header[0]} {reformed_header[0]}/{fw_or_rv}\n"
        return reformed_header

    def filter(self, infile, outfile, removed_file):
        """
        Writes the lines to a fastq file if it passes certain checks.
        """
        count = 0
        start = time.time()
        count_remove = 0
        if self.direction == "2":
            number = 2
        elif self.direction == "1":
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
                # break

        return count

def argparser():
    parser = argparse.ArgumentParser(description="Filter FASTQ sequences with too many consecutive 'G's.")
    parser.add_argument("input_file", type=str, help="Path to the input FASTQ file (gzipped).")
    parser.add_argument("output", type=str, help="Path to the output FASTQ file (gzipped).")
    parser.add_argument("--poly_g", type=int,
                        help="Minimum number of consecutive 'G's to filter out sequences.", default=30)
    args = parser.parse_args()
    return args



def main():
    """
    Executes the program functions in order.
    """
    args = argparser()
    filterer = FastqFilter(input_file=args.input_file, poly_g=args.poly_g, output=args.output)
    filterer.filter_fastq()

if __name__ == '__main__':
    sys.exit(main())