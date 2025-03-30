# init
# ----

import subprocess
import gzip

from utils import remove_previous_log_files, check_for_previous_output, check_for_missing_otu_tables, fail_if_max_attempts

# logging
# -------
from log_config import configure_logger
import os

logger_name = os.path.splitext(os.path.basename(__file__))[0]
logger      = configure_logger(logger_name)

# snakefile parameters
# --------------------
SINGLETON               = snakemake.input['singleton']

INPUT_DONE              = snakemake.input['input_done']
INPUT_OTU               = INPUT_DONE.replace('.done', '.tsv.gz')

METAPACKAGE             = snakemake.params['metapackage']
TAX_PROFILE             = snakemake.params['tax_profile']

MODE                    = snakemake.wildcards['mode']
if MODE == 'plastic':
    MODE = 'singlem'

ATTEMPT                 = snakemake.resources['attempt']
TOTAL_ATTEMPTS          = snakemake.resources['total_attempts']

THREADS                 = snakemake.threads

if 'sam_accession' in snakemake.wildcards.keys():
    ACCESSION           = snakemake.wildcards['sam_accession']
    ACC_TYPE            = 'sam'
else:
    ACCESSION           = snakemake.wildcards['run_accession']
    ACC_TYPE            = 'run'


# functions
# ---------
def assign_taxonomy(
        input_otu, 
        tax_profile, 
        mode,
        metapackage, 
        threads=1):
    """
    Assigns taxonomy to OTU table based on sample accession.
    """

    cmd = f'{mode} renew ' \
        f'--input-archive-otu-table <(gunzip -c {input_otu}) ' \
        f'--taxonomic-profile >(gzip >{tax_profile}) ' \
        f'--metapackage {metapackage} ' \
        f'--threads {threads}'
        
    logger.info(f"Attempting SingleM command: {cmd}")
    output = subprocess.run(
        ["bash",'-o','pipefail',"-c", cmd],
        capture_output=True, text=True)

    print(output.stdout + output.stderr)

    if "had too many unassigned OTUs" in output.stdout + output.stderr:

        logger.info("No reads found in the output.")
        logger.info("Updating no reads found file...")
        with open(os.path.join(os.path.dirname(input_otu), 'no_reads_found'), 'a') as f: pass

        logger.info("Writing empty taxonomic profile...")
        with gzip.open(tax_profile, 'wt') as f:
            f.write('sample\tcoverage\ttaxonomy\n') # type: ignore

    elif output.returncode != 0:
        raise(Exception("Error in SingleM command."))


def symlink_run_taxonomy(
    input_tax_profile, 
    output_tax_profile
    ):

    input_tax_profile       = os.path.abspath(input_tax_profile)
    output_tax_profile      = os.path.abspath(output_tax_profile)

    os.symlink(input_tax_profile, output_tax_profile)

# main
def main():
    
    logger.info(f"Assigning taxonomy for {ACCESSION}...")

    # housekeeping
    OUTPUT_DIR = os.path.dirname(TAX_PROFILE)

    # remove previous log files
    remove_previous_log_files(OUTPUT_DIR)
    # check if the output file already exists
    check_for_previous_output(TAX_PROFILE)
    # check for missing OTU tables
    check_for_missing_otu_tables(INPUT_OTU, OUTPUT_DIR, ACCESSION)
    # fail if max attempts reached
    fail_if_max_attempts(ATTEMPT, TOTAL_ATTEMPTS, OUTPUT_DIR, ACCESSION)


    if SINGLETON:
        logger.info(f"Only one run found for {ACCESSION}...")
        symlink_run_taxonomy(INPUT_OTU, TAX_PROFILE)
        logger.info(f"Taxonomic profile {INPUT_OTU} symlinked to {TAX_PROFILE}.")
    else:
        if ACC_TYPE == 'sam':
            logger.info(f"Multiple runs found for {ACCESSION}.")

        logger.info(f"Assigning taxonomy to OTU table {INPUT_OTU}...")
        assign_taxonomy(INPUT_OTU, TAX_PROFILE, MODE, METAPACKAGE, THREADS)
        logger.info(f"Taxonomic profile for {ACCESSION} written to {TAX_PROFILE}.")

    logger.info("Done!")

if __name__ == '__main__':
    main()