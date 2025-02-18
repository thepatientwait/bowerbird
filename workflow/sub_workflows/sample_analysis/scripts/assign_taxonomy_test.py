# init
# ----

import subprocess
import polars as pl
import argparse
import os
import re
import gzip
import psutil

import time

import threading

# logging
# -------
from log_config import configure_logger

logger_name = os.path.splitext(os.path.basename(__file__))[0]
logger = configure_logger(logger_name)

# functions
# ---------
def read_output(stream, buffer):
    for line in iter(stream.readline, ''):
        buffer.append(line)
        print(line, end='')

        
def assign_taxonomy(input_otu, tax_profile, metapackage, threads=1):
    """
    Assigns taxonomy to OTU table based on sample accession.
    """

    cmd = f'singlem renew ' \
        f'--input-archive-otu-table <(gunzip -c {input_otu}) ' \
        f'--taxonomic-profile >(gzip >{tax_profile}) ' \
        f'--metapackage {metapackage} ' \
        f'--threads {threads}'
    
    logger.info(f"Attempting SingleM command: {cmd}")
    # output = subprocess.run(
    #     ["bash",'-o','pipefail',"-c", cmd],
    #     capture_output=True, text=True)

    # print(output.stdout + output.stderr)

    # if "had too many unassigned OTUs" in output.stdout + output.stderr:

    #     logger.info("No reads found in the output.")
    #     logger.info("Updating no reads found file...")
    #     with open(os.path.join(os.path.dirname(input_otu), 'no_reads_found'), 'a') as f: pass

    #     logger.info("Writing empty taxonomic profile...")
    #     with gzip.open(tax_profile, 'wt') as f:
    #         f.write('sample\tcoverage\ttaxonomy\n') # type: ignore

    # elif output.returncode != 0:
    #     raise(Exception("Error in SingleM command."))

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


def check_for_empty_otu_table(
    input_otu,
    tax_profile
    ):
    dirpath = os.path.dirname(input_otu)
    file = os.path.join(dirpath, 'no_reads_found')

    if os.path.exists(file):
        logger.info("No reads found in the output.")
        logger.info("Writing empty taxonomic profile...")
        with gzip.open(tax_profile, 'wt') as f:
            f.write('sample\tcoverage\ttaxonomy\n') # type: ignore

        return True          
    
    else:
        return False  
       

# main
def main():
    parser = argparse.ArgumentParser(description="Assign taxonomy to OTU table.")
    parser.add_argument('--input-otu', required=True, help="Path to the input OTU table.")
    parser.add_argument('--tax-profile', required=True, help="Path to the output taxonomic profile.")

    metapackage = '/mnt/weka/pkg/cmr/n9985158/singlem-metapackage/S4.3.0.GTDB_r220.metapackage_20240523.smpkg.zb'

    args = parser.parse_args()

    logger.info(f"Assigning taxonomy to OTU table {args.input_otu}...")
    is_empty = check_for_empty_otu_table(args.input_otu, args.tax_profile)
    if is_empty:
        logger.info(f"Taxonomic profile written to {args.tax_profile}.")
        logger.info("Done!")
    else:
        assign_taxonomy(args.input_otu, args.tax_profile, metapackage, 50)
        logger.info(f"Taxonomic profile written to {args.tax_profile}.")
        logger.info("Done!")

if __name__ == '__main__':
    main()