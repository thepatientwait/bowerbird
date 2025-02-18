# init
# ----

import extern
import polars as pl
import gzip
import glob

# logging
# -------
from log_config import configure_logger
import os

logger_name = os.path.splitext(os.path.basename(__file__))[0]
logger      = configure_logger(logger_name)

# snakefile parameters
# --------------------
RUN_ACCESSION           = snakemake.wildcards['run_accession']

THREADS                 = snakemake.threads

OUTPUT_DIR              = snakemake.output['output_dir']

# functions
# ---------

def download_run(
    run_accession, 
    output_dir,
    threads
    ):

    cmd = f'kingfisher get ' \
        f'-r {run_accession} ' \
        f'--output-directory {output_dir} ' \
        f'--download-threads {threads} ' \
        f'--extraction-threads {threads} ' \
        f'--output-format-possibilities fastq.gz --hide-download-progress ' \
        f'-m ena-ascp ena-ftp prefetch aws-http aws-cp'

    logger.info(f"Downloading reads for {run_accession}...")
    try:
        logger.info(f"Attempting SingleM command: {cmd}")
        extern.run(cmd)
    except extern.ExternCalledProcessError as e:
        logger.error(e)
        raise(e)


def check_files(
    run_accession, 
    output_dir
    ):

    if os.path.exists(f'{output_dir}/{run_accession}_1.fastq.gz'):
        if os.path.exists(f'{output_dir}/{run_accession}_2.fastq.gz'):
            return True, f"Found paired read files {run_accession}_1.fastq.gz and {run_accession}_2.fastq.gz"
        else:
            os.remove(f'{output_dir}/{run_accession}_1.fastq.gz')
            raise Exception("Found a forward read file but no reverse read file from ENA")
    else:
        return False, f"Found single read file {run_accession}.fastq.gz"


def check_interleaved(
    fastq_reads, 
    START_CHECK_PAIRS=5, 
    END_CHECK_PAIRS=5
    ):
    """
    Return True if:
     - reads have duplicate metadata in consecutive reads
     - total_count is even
    """

    start_list = []
    end_list = []
    read_count = 0

    with gzip.open(fastq_reads, "rt") as f:

        for line_count, line in enumerate(f):
            if line_count % 4 != 0:
                continue

            read_count += 1
            if len(start_list) < START_CHECK_PAIRS * 2:
                start_list.append([line.strip(), read_count])

            end_list.append([line.strip(), read_count])
            if len(end_list) > END_CHECK_PAIRS * 2:
                end_list.pop(0)

    try:
        line_count # type: ignore
    except UnboundLocalError:
        return False, "Empty file"

    READS_COLUMNS = {
        "read": str,
        "number": int,
        }

    reads = pl.DataFrame(start_list + end_list, orient="row", schema=READS_COLUMNS)


    READ_METADATA_REGEX = r"[^ ]*(.*)"

    if read_count % 2 == 1:
        return False, "Odd readcount"

    processed = (
        reads
        .select(
            pl.col("number").add(1).floordiv(2).alias("index"),
            pl.col("number").mod(2).alias("rank"),
            pl.col("read").str.extract(READ_METADATA_REGEX)
            )
        .pivot(values="read", index="index", on="rank", aggregate_function=None)
        .filter(pl.col("0") != pl.col("1"))
    )

    if processed.height > 0:
        return False, f"Consecutive reads do not match ({processed.height}/{reads.height // 2})"
    else:
        return True, "Duplicate read names in consecutive reads with even readcount"

# main
def main():

    logger.info(f"Downloading run {RUN_ACCESSION}...")

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    download_run(RUN_ACCESSION, OUTPUT_DIR, THREADS)

    # remove files '*aria2' that are created by kingfisher
    for f in glob.glob(os.path.join(OUTPUT_DIR, '*.aria2')):
        os.remove(f)

    two_files, reason = check_files(RUN_ACCESSION, OUTPUT_DIR)

    if two_files:
        logger.info(f"Run {RUN_ACCESSION} has separate paired read files.")
    else:
        logger.info(f"Run {RUN_ACCESSION} has single file.")

        interleaved, reason = check_interleaved(f'{OUTPUT_DIR}/{RUN_ACCESSION}.fastq.gz')

        logger.info(f"Run {RUN_ACCESSION} has {'non-' if not interleaved else ''}interleaved reads: {reason}.")

        with open(OUTPUT_DIR + "/interleaved", "w") as f:
            f.write(f"{interleaved}\t{reason}\n")


    logger.info(f"Downloaded reads for {RUN_ACCESSION}.")

if __name__ == '__main__':
    main()