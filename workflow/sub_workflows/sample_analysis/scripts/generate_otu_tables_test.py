# init
# ----

import extern
import os
import argparse
from memory_profiler import profile

import subprocess

# logging
# -------
from log_config import configure_logger
import os

logger_name = os.path.splitext(os.path.basename(__file__))[0]
logger      = configure_logger(logger_name)

# functions
# ---------

def analyse(run_accession, run_otu_table, metapackage, threads):
    some_ena_analysed = False
    paired_result = None
    unpaired_result = None
    if os.path.exists(f'{run_accession}_1.fastq.gz'):
        if os.path.exists(f'{run_accession}_2.fastq.gz'):
            sequence_input_arg = f'--forward {run_accession}_1.fastq.gz --reverse {run_accession}_2.fastq.gz'
            some_ena_analysed = True
            paired_result = run_singlem(sequence_input_arg, run_otu_table, metapackage, threads)
        else:
            raise Exception("Found a forward read file but no reverse read file from ENA")
    
    if os.path.exists(f'{run_accession}.fastq.gz'):
        sequence_input_arg = f'--forward {run_accession}.fastq.gz'
        some_ena_analysed = True
        unpaired_result = run_singlem(sequence_input_arg, run_otu_table, metapackage, threads)

    if not some_ena_analysed:
        raise Exception("Unexpected form of ENA download")


def run_singlem(sequence_input_arg, run_otu_table, metapackage, threads):
    cmd = f'singlem pipe {sequence_input_arg} ' \
        f'--no-assign-taxonomy ' \
        f'--archive-otu-table {run_otu_table} ' \
        f'--metapackage {metapackage} ' \
        f'--threads {threads}'

    logger.info(f"Attempting SingleM command: {cmd}")
    output = subprocess.run(
        ["bash",'-o','pipefail',"-c", cmd],
        capture_output=True, text=True)

    print(output.stdout + output.stderr)

    if output.returncode != 0:
        raise(Exception("Error in SingleM command."))
        
    if "No reads found" in output.stdout + output.stderr:
        print("No reads found in the output.")
        output_dir = os.path.dirname(run_otu_table)
        with open(os.path.join(output_dir, "no_reads_found"), "a") as f: pass

# main
def main():
    parser = argparse.ArgumentParser(description="Generate OTU table.")
    parser.add_argument('--dir', required=False, help="Path to the run directory.", default=os.getcwd())
    parser.add_argument('--accession', required=True, help="Run accession.")
    parser.add_argument('--output_otu_table', required=True, help="Path to the output OTU table.")

    args = parser.parse_args()

    metapackage = '/mnt/weka/pkg/cmr/n9985158/singlem-metapackage/S4.3.0.GTDB_r220.metapackage_20240523.smpkg.zb'
    threads = 1

    os.chdir(args.dir)

    logger.info(f"Analysing run...")
    
    analyse(args.accession, args.output_otu_table, metapackage, threads)

    logger.info("Done!")

if __name__ == '__main__':
    main()