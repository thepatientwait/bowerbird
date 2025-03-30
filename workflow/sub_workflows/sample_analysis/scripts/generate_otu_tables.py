# init
# ----

import subprocess
import gzip

from utils import remove_previous_log_files, check_for_previous_output, fail_if_max_attempts

# logging
# -------
from log_config import configure_logger
import os

logger_name = os.path.splitext(os.path.basename(__file__))[0]
logger      = configure_logger(logger_name)

# snakefile parameters
# --------------------
RUN_DIR                 = snakemake.input['run_dir']

METAPACKAGE             = snakemake.params['metapackage']
RUN_OTU_TABLE           = snakemake.params['run_otu_table']

ATTEMPT                 = snakemake.resources['attempt']
TOTAL_ATTEMPTS          = snakemake.resources['total_attempts']

MODE                    = snakemake.wildcards['mode']
if MODE == 'plastic':
    MODE = 'singlem'

RUN_ACCESSION           = snakemake.wildcards['run_accession']


THREADS                 = snakemake.threads


# functions
# ---------

def analyse(
    run_accession, 
    run_otu_table, 
    mode, 
    metapackage, 
    threads
    ):

    some_ena_analysed = False

    if os.path.exists(f'{run_accession}_1.fastq.gz'):
        if os.path.exists(f'{run_accession}_2.fastq.gz'):
            sequence_input_arg = f'--forward {run_accession}_1.fastq.gz --reverse {run_accession}_2.fastq.gz'
            run_singlem(sequence_input_arg, run_otu_table, mode, metapackage, threads)
        else:
            raise Exception("Found a forward read file but no reverse read file from ENA")
    
    if os.path.exists(f'{run_accession}.fastq.gz'):
        sequence_input_arg = f'--forward {run_accession}.fastq.gz'
        run_singlem(sequence_input_arg, run_otu_table, mode, metapackage, threads)

def run_singlem(
    sequence_input_arg, 
    run_otu_table, 
    mode, 
    metapackage, 
    threads
    ):
    
    cmd = f'{mode} pipe {sequence_input_arg} ' \
        f'--no-assign-taxonomy ' \
        f'--archive-otu-table >(gzip >{run_otu_table}) ' \
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
        logger.info("No reads found in the input.")
        output_dir = os.path.dirname(run_otu_table)
        with open(os.path.join(output_dir, "log"), "w") as f:
            f.write(f"No reads found in the input.")
        with open(os.path.join(output_dir, "failed"), "w") as f:
            f.write(f"{RUN_ACCESSION}\n")
        logger.info("Logged failed run accession.")
        logger.info("Exiting without error.")
        exit(0)
        

# main
def main():

    logger.info(f"Analysing run {RUN_ACCESSION}...")

    # housekeeping
    OUTPUT_DIR = os.path.dirname(RUN_OTU_TABLE)

    # remove previous log files
    remove_previous_log_files(OUTPUT_DIR)
    # check if the output file already exists
    check_for_previous_output(RUN_OTU_TABLE)
    # fail if max attempts reached
    fail_if_max_attempts(ATTEMPT, TOTAL_ATTEMPTS, OUTPUT_DIR, RUN_ACCESSION)

    run_path = os.path.join(RUN_DIR, RUN_ACCESSION)
    analyse(run_path, RUN_OTU_TABLE, MODE, METAPACKAGE, THREADS)
    logger.info("Done!")

if __name__ == '__main__':
    main()