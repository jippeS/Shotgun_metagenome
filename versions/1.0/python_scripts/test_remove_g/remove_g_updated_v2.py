#!/usr/bin/env python3

import os
import sys
import argparse
import multiprocessing as mp
import tempfile
import subprocess
import time


# --------------------------------------------------
# Header reform (same logic as original)
# --------------------------------------------------

def reform_header(header, fw_or_rv):
    header = header.strip().split(" ")
    return f"{header[0]} {header[0]}/{fw_or_rv}\n"


# --------------------------------------------------
# Filter one chunk
# --------------------------------------------------

def filter_chunk(args):
    chunk_path, poly_g, direction = args

    poly_string = "G" * poly_g

    filtered_path = f"{chunk_path}.filtered"
    removed_path = f"{chunk_path}.removed"

    kept = 0
    removed = 0

    with open(chunk_path, "r") as infile, \
         open(filtered_path, "w") as outfile, \
         open(removed_path, "w") as removed_file:

        while True:
            header = infile.readline()
            seq = infile.readline()
            plus = infile.readline()
            qual = infile.readline()

            if not qual:
                break

            header = reform_header(header, direction)

            if poly_string not in seq:
                outfile.write(header)
                outfile.write(seq)
                outfile.write(plus)
                outfile.write(qual)
                kept += 1
            else:
                removed_file.write(header)
                removed_file.write(seq)
                removed_file.write(plus)
                removed_file.write(qual)
                removed += 1

    return filtered_path, removed_path, kept, removed


# --------------------------------------------------
# Split FASTQ safely
# --------------------------------------------------

def split_fastq(temp_fastq, reads_per_chunk):
    lines_per_chunk = reads_per_chunk * 4
    chunk_files = []
    chunk_index = 0

    with open(temp_fastq, "r") as infile:

        while True:
            lines = []
            for _ in range(lines_per_chunk):
                line = infile.readline()
                if not line:
                    break
                lines.append(line)

            if not lines:
                break

            chunk_path = f"{temp_fastq}.chunk{chunk_index}"

            with open(chunk_path, "w") as chunk:
                chunk.writelines(lines)

            chunk_files.append(chunk_path)
            chunk_index += 1

    return chunk_files


# --------------------------------------------------
# Generate output filenames (same logic as original)
# --------------------------------------------------

def generate_output_paths(input_file, output_dir):

    if not output_dir.endswith("/"):
        output_dir += "/"

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    filename = os.path.basename(input_file)

    base, ext = os.path.splitext(filename)

    if ext == ".gz":
        base, ext2 = os.path.splitext(base)
        output_file = f"{output_dir}{base}_filtered{ext2}{ext}"
        removed_file = f"{output_dir}{base}_removed_seqs{ext2}{ext}"
    else:
        output_file = f"{output_dir}{base}_filtered{ext}"
        removed_file = f"{output_dir}{base}_removed_seqs{ext}"

    return output_file, removed_file


# --------------------------------------------------
# Detect direction from filename (_1 or _2)
# --------------------------------------------------

def detect_direction(input_file):
    filename = os.path.basename(input_file)

    parts = filename.split("_")
    if len(parts) < 2:
        sys.exit("Filename does not contain _1 or _2")

    direction = parts[1].split(".")[0]

    if direction == "1":
        return 1
    elif direction == "2":
        return 2
    else:
        sys.exit(f"Could not determine direction from filename: {filename}")


# --------------------------------------------------
# Main parallel driver
# --------------------------------------------------

def parallel_filter_fastq(input_gz, output_dir, poly_g, threads, reads_per_chunk):

    start_total = time.time()

    if threads is None:
        threads = mp.cpu_count()

    direction = detect_direction(input_gz)

    output_gz, removed_gz = generate_output_paths(input_gz, output_dir)

    with tempfile.TemporaryDirectory() as tmpdir:

        temp_fastq = os.path.join(tmpdir, "decompressed.fastq")

        # Decompress
        print("Decompressing with pigz...")
        with open(temp_fastq, "wb") as outfile:
            subprocess.run(
                ["pigz", "-dc", input_gz],
                stdout=outfile,
                check=True
            )

        # Split
        print("Splitting FASTQ into chunks...")
        chunks = split_fastq(temp_fastq, reads_per_chunk)

        print(f"Created {len(chunks)} chunks.")
        print(f"Filtering using {threads} CPU cores...")

        tasks = [(chunk, poly_g, direction) for chunk in chunks]

        with mp.Pool(threads) as pool:
            results = pool.map(filter_chunk, tasks)

        # Merge filtered + removed
        merged_filtered = os.path.join(tmpdir, "merged_filtered.fastq")
        merged_removed = os.path.join(tmpdir, "merged_removed.fastq")

        total_kept = 0
        total_removed = 0

        with open(merged_filtered, "wb") as out_filt, \
             open(merged_removed, "wb") as out_rem:

            for filt, rem, kept, removed in results:
                total_kept += kept
                total_removed += removed

                with open(filt, "rb") as f:
                    out_filt.write(f.read())

                with open(rem, "rb") as f:
                    out_rem.write(f.read())

        # Compress filtered
        print("Compressing filtered output...")
        with open(output_gz, "wb") as outfile:
            subprocess.run(
                ["pigz", "-p", str(threads), "-c", merged_filtered],
                stdout=outfile,
                check=True
            )

        # Compress removed
        print("Compressing removed sequences...")
        with open(removed_gz, "wb") as outfile:
            subprocess.run(
                ["pigz", "-p", str(threads), "-c", merged_removed],
                stdout=outfile,
                check=True
            )

    end_total = time.time()

    print("\nFinished.")
    print(f"Sequences kept: {total_kept}")
    print(f"Sequences removed: {total_removed}")
    print(f"Runtime: {round(end_total - start_total, 2)} seconds")


# --------------------------------------------------
# CLI
# --------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Parallel FASTQ poly-G filter (same output format as original)"
    )

    parser.add_argument("input", help="Input FASTQ.gz file")
    parser.add_argument("output_dir", help="Output directory")

    parser.add_argument("--poly_g", type=int, default=30)
    parser.add_argument("--threads", type=int, default=None)
    parser.add_argument("--reads_per_chunk", type=int, default=1_000_000)

    return parser.parse_args()


def main():
    args = parse_args()

    parallel_filter_fastq(
        input_gz=args.input,
        output_dir=args.output_dir,
        poly_g=args.poly_g,
        threads=args.threads,
        reads_per_chunk=args.reads_per_chunk
    )


if __name__ == "__main__":
    sys.exit(main())