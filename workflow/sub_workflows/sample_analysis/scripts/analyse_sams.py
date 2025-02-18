# init
# ----

import subprocess
import polars as pl
import gzip
import tempfile

# logging
# -------
from log_config import configure_logger
import os

logger_name = os.path.splitext(os.path.basename(__file__))[0]
logger      = configure_logger(logger_name)

# snakefile parameters
# --------------------
INPUT_OTUS              = snakemake.input['otu_tables']
SINGLETON               = snakemake.input['singleton']

METAPACKAGE             = snakemake.params['metapackage']

THREADS                 = snakemake.threads

SAM_ACCESSION           = snakemake.wildcards['sam_accession']

SAM_OTU_TABLE           = snakemake.output['otu_table']
SAM_TAX_PROFILE         = snakemake.output['tax_profile']

# functions
# ---------

def symlink_otu_table(
    input_otu,
    sam_otu_table
    ):
    input_otu       = os.path.abspath(input_otu)
    sam_otu_table   = os.path.abspath(sam_otu_table)

    os.symlink(input_otu, sam_otu_table)


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

    logger.info(f"Excluding empty OTU tables: {bad_otus}")

    return good_otus


def merge_otu_tables_command(otu_tables, sam_accession, sam_otu_table, metapackage):
    with tempfile.NamedTemporaryFile() as f:
        f.write('\n'.join(otu_tables).encode())
        f.flush()
        input_otu_list = f.name

        cmd = f'singlem summarise ' \
            f'--input-gzip-archive-otu-table-list {input_otu_list} ' \
            f'--collapse-to-sample-name {sam_accession} ' \
            f'--output-archive-otu-table >(gzip >{sam_otu_table}) ' \
            f'--metapackage {metapackage}'

        logger.info(f"Attempting SingleM command: {cmd}")

        output = subprocess.run(
            ["bash",'-o','pipefail',"-c", cmd],
            capture_output=True, text=True)

        print(output.stdout + output.stderr)
        return output
    

def merge_otu_tables(
    input_otus, 
    sam_accession, 
    sam_otu_table, 
    metapackage
    ):

    # First attempt with input_otus
    input_otus = remove_empty_otu_tables(input_otus)
    output = merge_otu_tables_command(input_otus, sam_accession, sam_otu_table, metapackage)

    if "ValueError: Length mismatch:" in output.stdout + output.stderr:
        good_otus = []
        empty_otus = 0
        for otu_table in input_otus:
            if os.path.getsize(otu_table) == 0:
                file_name = os.path.basename(otu_table)
                logger.info(f"OTU table for {file_name} is empty. Removing...")
                os.remove(otu_table)
                empty_otus += 1
            else:
                with gzip.open(otu_table) as file:
                    if (pl.read_json(file.read())['otus'].dtype == pl.List(pl.Null)):
                        file_name = os.path.basename(otu_table)
                        logger.info(f"No reads found in {file_name}. Updating no reads found file...")
                        output_dir = os.path.dirname(otu_table)
                        with open(os.path.join(output_dir, "no_reads_found"), "a") as f: pass
                    else:
                        good_otus.append(otu_table)
        
        if empty_otus > 0:
            raise Exception(f"There were {empty_otus} empty OTU tables.")

        # Retry with good_otus
        if good_otus:
            output = merge_otu_tables_command(good_otus, sam_accession, sam_otu_table, metapackage)

    if output.returncode != 0:
        raise Exception("Error in SingleM command.")


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


# main
def main():

    logger.info(f"Merging OTU tables for {SAM_ACCESSION}...")

    # check if all otu tables exist (mostly here for if the rule retries)
    for input_otu in INPUT_OTUS:
        if not os.path.exists(input_otu):
            raise Exception(f"OTU table {input_otu} not found.")
    
    if len(INPUT_OTUS) == 1:
        logger.info(f"Only one OTU table found for {SAM_ACCESSION}...")
        symlink_otu_table(INPUT_OTUS[0], SAM_OTU_TABLE)
        logger.info(f"OTU table for {SAM_ACCESSION} symlinked to {SAM_OTU_TABLE}")
    else:
        logger.info(f"Multiple OTU tables found for {SAM_ACCESSION}...")
        merge_otu_tables(INPUT_OTUS, SAM_ACCESSION, SAM_OTU_TABLE, METAPACKAGE)
        logger.info(f"Merged OTU table for {SAM_ACCESSION} written to {SAM_OTU_TABLE}")

        

    logger.info(f"Assigning taxonomy for {SAM_ACCESSION}...")

    is_empty = check_for_empty_otu_table(SAM_OTU_TABLE, SAM_TAX_PROFILE)

    if is_empty:
        logger.info(f"Taxonomic profile written to {SAM_TAX_PROFILE}.")
        logger.info("Done!")
    elif SINGLETON:
        logger.info(f"Only one run found for {SAM_ACCESSION}...")
        symlink_run_taxonomy(SINGLETON[0], SAM_TAX_PROFILE)
        logger.info(f"Taxonomic profile {SINGLETON} symlinked to {SAM_TAX_PROFILE}.")
    else:
        logger.info(f"Multiple runs found for {SAM_ACCESSION}.")

        logger.info(f"Assigning taxonomy to OTU table {SAM_OTU_TABLE}...")
        assign_taxonomy(SAM_OTU_TABLE, SAM_TAX_PROFILE, METAPACKAGE, THREADS)
        logger.info(f"Taxonomic profile for {SAM_ACCESSION} written to {SAM_TAX_PROFILE}.")

    logger.info("Done!")

if __name__ == '__main__':
    main()