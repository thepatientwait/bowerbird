# TODO: add check for sra_runs which predate the current pipeline run

# init
# ----

import subprocess
import polars as pl
import gzip

import tempfile

import argparse

import psutil

import time

import threading

# logging
# -------
from log_config import configure_logger
import os

logger_name = os.path.splitext(os.path.basename(__file__))[0]
logger      = configure_logger(logger_name)


# functions
# ---------

def symlink_otu_table(
    input_otu,
    output_otu
    ):
    input_otu       = os.path.abspath(input_otu)
    output_otu      = os.path.abspath(output_otu)

    os.symlink(input_otu, output_otu)


def remove_empty_otu_tables(
    input_otus
    ):
    good_otus = []
    bad_otus = []

    for input_otu in input_otus:
        dirpath = os.path.dirname(input_otu)
        if os.path.exists(os.path.join(dirpath, 'no_reads_found')):
            bad_otus.append(input_otu)
        else:
            good_otus.append(input_otu)

    logger.info(f"Removed empty OTU tables: {bad_otus}")

    return good_otus

def read_output(stream, buffer):
    for line in iter(stream.readline, ''):
        buffer.append(line)
        print(line, end='')

# @profile
def run_singlem_command(otu_tables, sam_accession, output_otu, metapackage):
    with tempfile.NamedTemporaryFile() as f:
        f.write('\n'.join(otu_tables).encode())
        f.flush()
        input_otu_list = f.name

        cmd = f'singlem summarise ' \
            f'--input-gzip-archive-otu-table-list {input_otu_list} ' \
            f'--collapse-to-sample-name {sam_accession} ' \
            f'--output-archive-otu-table >(gzip >{output_otu}) ' \
            f'--metapackage {metapackage}'

        logger.info(f"Attempting SingleM command: {cmd}")

        process = subprocess.Popen(
            ["bash", '-o', 'pipefail', "-c", cmd],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        # Track memory usage of the process and its children
        memory_usage = []

        stdout_buffer = []
        stderr_buffer = []

        stdout_thread = threading.Thread(target=read_output, args=(process.stdout, stdout_buffer))
        stderr_thread = threading.Thread(target=read_output, args=(process.stderr, stderr_buffer))

        stdout_thread.start()
        stderr_thread.start()

        while process.poll() is None:
            try:
                proc = psutil.Process(process.pid)
                mem_info = proc.memory_info()
                memory_usage.append(mem_info.rss)
                for child in proc.children(recursive=True):
                    mem_info = child.memory_info()
                    memory_usage.append(mem_info.rss)
            except psutil.NoSuchProcess:
                pass
            time.sleep(0.1)

        stdout_thread.join()
        stderr_thread.join()

        stdout = ''.join(stdout_buffer)
        stderr = ''.join(stderr_buffer)

        print(stdout + stderr)

        max_memory_usage = max(memory_usage) if memory_usage else 0
        logger.info(f"Max memory usage: {max_memory_usage / (1024 * 1024)} MB")

        return process.returncode
    

def merge_otu_tables(
    input_otus, 
    sam_accession, 
    output_otu, 
    metapackage
    ):

    # First attempt with input_otus
    # input_otus = remove_empty_otu_tables(input_otus)
    output = run_singlem_command(input_otus, sam_accession, output_otu, metapackage)

    # if "ValueError: Length mismatch:" in output.stdout + output.stderr:
    #     good_otus = []
    #     for otu_table in input_otus:
    #         if os.path.getsize(otu_table) == 0:
    #             file_name = os.path.basename(otu_table)
    #             logger.info(f"OTU table for {file_name} is empty. Removing...")
    #             # os.remove(otu_table)
    #         else:
    #             with gzip.open(otu_table) as file:
    #                 if (pl.read_json(file.read())['otus'].dtype == pl.List(pl.Null)):
    #                     file_name = os.path.basename(otu_table)
    #                     logger.info(f"No reads found in {file_name}. Updating no reads found file...")
    #                     output_dir = os.path.dirname(otu_table)
    #                     print(output_dir)
    #                     # with open(os.path.join(output_dir, "no_reads_found"), "a") as f: pass
    #                 else:
    #                     good_otus.append(otu_table)

    #     # Retry with good_otus
    #     if good_otus:
    #         output = run_singlem_command(good_otus, sam_accession, output_otu, metapackage)

    # if output.returncode != 0:
    #     raise Exception("Error in SingleM command.")


# main
def main():
    parser = argparse.ArgumentParser(description="Assign taxonomy to OTU table.")
    parser.add_argument('--input-otus', required=True, help="Path to the input OTU tables list.")
    parser.add_argument('--output-otu', required=True, help="Path to the output taxonomic profile.")

    metapackage = '/mnt/weka/pkg/cmr/n9985158/singlem-metapackage/S4.3.0.GTDB_r220.metapackage_20240523.smpkg.zb'

    args = parser.parse_args()

    logger.info(f"Merging OTU tables for {"test"}...")

    with open(args.input_otus) as f:
        input_otus = f.readlines()
        input_otus = [x.strip() for x in input_otus]
    
    if len(input_otus) == 1:
        logger.info(f"Only one OTU table found for {"test"}...")
        symlink_otu_table(input_otus[0], args.output_otu)
        logger.info(f"OTU table for {"test"} symlinked to {args.output_otu}")
    else:
        logger.info(f"Multiple OTU tables found for {"test"}...")
        merge_otu_tables(input_otus, "test", args.output_otu, metapackage)

        logger.info(f"Merged OTU table for {"test"} written to {args.output_otu}")

    logger.info(f"Done!")

if __name__ == '__main__':
    main()