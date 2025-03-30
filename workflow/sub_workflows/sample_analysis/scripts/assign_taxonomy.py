# init
# ----

import subprocess
import gzip

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
INPUT_OTU               = INPUT_DONE.replace('.done', '.json.gz')

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
    
def check_for_missing_otu_table(
    input_otu,
    output_dir,
    accession
    ):
    
    if not os.path.exists(input_otu):
        logger.info(f"OTU table for {input_otu} is missing.")
        with open(os.path.join(output_dir, "log"), "w") as f:
            f.write(f"OTU table for {input_otu} is missing.")
        with open(os.path.join(output_dir, "failed"), "w") as f:
            f.write(f"{accession}\n")
        logger.info("Logged failed accession.")
        logger.info("Exiting without error.")
        exit(0)


def fail_if_max_attempts(
    attempt,
    total_attempts,
    output_dir,
    accession
    ):    

    if attempt <= total_attempts:
        return

    # If max attempts are reached, log failure and exit
    logger.info("No more attempts! OTU table generation failed.")
    with open(os.path.join(output_dir, "log"), "w") as f:
        f.write(f"No more attempts! OTU table generation failed.")
    with open(os.path.join(output_dir, "failed"), "w") as f:
        f.write(f"{accession}\n")
    logger.info("Logged failed run accession.")
    logger.info("Exiting without error.")
    exit(0)

# main
def main():
    
    logger.info(f"Assigning taxonomy for {ACCESSION}...")

    is_empty = check_for_empty_otu_table(INPUT_OTU, TAX_PROFILE)

    # check for missing OTU tables
    check_for_missing_otu_table(INPUT_OTU, os.path.dirname(TAX_PROFILE), ACCESSION)
    # fail if max attempts reached
    fail_if_max_attempts(ATTEMPT, TOTAL_ATTEMPTS, os.path.dirname(TAX_PROFILE), ACCESSION)

    if is_empty:
        logger.info(f"Taxonomic profile written to {TAX_PROFILE}.")
        logger.info("Done!")
    elif SINGLETON:
        logger.info(f"Only one run found for {ACCESSION}...")
        symlink_run_taxonomy(SINGLETON[0], TAX_PROFILE)
        logger.info(f"Taxonomic profile {SINGLETON} symlinked to {TAX_PROFILE}.")
    else:
        if ACC_TYPE == 'sam':
            logger.info(f"Multiple runs found for {ACCESSION}.")

        logger.info(f"Assigning taxonomy to OTU table {INPUT_OTU}...")
        assign_taxonomy(INPUT_OTU, TAX_PROFILE, MODE, METAPACKAGE, THREADS)
        logger.info(f"Taxonomic profile for {ACCESSION} written to {TAX_PROFILE}.")

    logger.info("Done!")

if __name__ == '__main__':
    main()