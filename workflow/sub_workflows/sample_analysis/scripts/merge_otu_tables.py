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

METAPACKAGE             = snakemake.params['metapackage']

SAM_ACCESSION           = snakemake.wildcards['sam_accession']

OUTPUT_OTU              = snakemake.output['otu_table']

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

    logger.info(f"Excluding empty OTU tables: {bad_otus}")

    return good_otus


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

        output = subprocess.run(
            ["bash",'-o','pipefail',"-c", cmd],
            capture_output=True, text=True)

        print(output.stdout + output.stderr)
        return output
    

def merge_otu_tables(
    input_otus, 
    sam_accession, 
    output_otu, 
    metapackage
    ):

    # First attempt with input_otus
    input_otus = remove_empty_otu_tables(input_otus)
    output = run_singlem_command(input_otus, sam_accession, output_otu, metapackage)

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
            output = run_singlem_command(good_otus, sam_accession, output_otu, metapackage)

    if output.returncode != 0:
        raise Exception("Error in SingleM command.")


# main
def main():

    logger.info(f"Merging OTU tables for {SAM_ACCESSION}...")

    # check if all otu tables exist (mostly here for if the rule retries)
    for input_otu in INPUT_OTUS:
        if not os.path.exists(input_otu):
            raise Exception(f"OTU table {input_otu} not found.")
    
    if len(INPUT_OTUS) == 1:
        logger.info(f"Only one OTU table found for {SAM_ACCESSION}...")
        symlink_otu_table(INPUT_OTUS[0], OUTPUT_OTU)
        logger.info(f"OTU table for {SAM_ACCESSION} symlinked to {OUTPUT_OTU}")
    else:
        logger.info(f"Multiple OTU tables found for {SAM_ACCESSION}...")
        merge_otu_tables(INPUT_OTUS, SAM_ACCESSION, OUTPUT_OTU, METAPACKAGE)
        logger.info(f"Merged OTU table for {SAM_ACCESSION} written to {OUTPUT_OTU}")

    logger.info(f"Done!")

if __name__ == '__main__':
    main()