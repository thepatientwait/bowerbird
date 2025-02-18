# init
# ----

import subprocess

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

RUN_ACCESSION           = snakemake.wildcards['run_accession']
THREADS                 = snakemake.threads

RUN_OTU_TABLE           = snakemake.output['run_otu_table']
RUN_TAX_PROFILE         = snakemake.output['run_tax_profile']

# functions
# ---------

def analyse(run_accession, run_otu_table, run_tax_profile, metapackage, threads):
    some_ena_analysed = False

    if os.path.exists(f'{run_accession}_1.fastq.gz'):
        if os.path.exists(f'{run_accession}_2.fastq.gz'):
            sequence_input_arg = f'--forward {run_accession}_1.fastq.gz --reverse {run_accession}_2.fastq.gz'
            some_ena_analysed = True
            run_singlem(sequence_input_arg, run_otu_table, run_tax_profile, metapackage, threads)
        else:
            raise Exception("Found a forward read file but no reverse read file from ENA")
    
    if os.path.exists(f'{run_accession}.fastq.gz'):
        sequence_input_arg = f'--forward {run_accession}.fastq.gz'
        some_ena_analysed = True
        run_singlem(sequence_input_arg, run_otu_table, run_tax_profile, metapackage, threads)

    if not some_ena_analysed:
        raise Exception("Unexpected form of ENA download")

def run_singlem(sequence_input_arg, run_otu_table, run_tax_profile, metapackage, threads):
    cmd = f'singlem pipe {sequence_input_arg} ' \
        f'--taxonomic-profile >(gzip >{run_tax_profile}) ' \
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
        
    # this error handling was from when otu table and tax profile were generated by seperate scripts
    # may not work as expected
    if "No reads found" in output.stdout + output.stderr:
        logger.info("No reads found in the output.")
        output_dir = os.path.dirname(run_otu_table)
        with open(os.path.join(output_dir, "no_reads_found"), "a") as f: pass

    if "had too many unassigned OTUs" in output.stdout + output.stderr:

        logger.info("No reads found in the output.")
        logger.info("Updating no reads found file...")
        with open(os.path.join(os.path.dirname(input_otu), 'no_reads_found'), 'a') as f: pass

        logger.info("Writing empty taxonomic profile...")
        with gzip.open(tax_profile, 'wt') as f:
            f.write('sample\tcoverage\ttaxonomy\n') # type: ignore

# main
def main():

    logger.info(f"Analysing run {RUN_ACCESSION}...")

    run_path = os.path.join(RUN_DIR, RUN_ACCESSION)
    
    analyse(run_path, RUN_OTU_TABLE, RUN_TAX_PROFILE, METAPACKAGE, THREADS)

    logger.info("Done!")

if __name__ == '__main__':
    main()