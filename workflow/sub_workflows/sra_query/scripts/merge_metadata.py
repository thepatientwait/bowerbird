# init
# ----
import polars as pl
# logging
# -------
from log_config import configure_logger
import os

logger_name = os.path.splitext(os.path.basename(__file__))[0]
logger      = configure_logger(logger_name)

# snakefile parameters
# --------------------
NEW_SAMPLES                 = snakemake.input['new_samples']

OLD_DATE                    = snakemake.params['old_date']

TARGET                      = snakemake.wildcards['target']
NEW_DATE                    = snakemake.wildcards['new_date']

JSON                        = snakemake.output['json']

# functions
# ---------

# main
def main():
    logger.info("Merging metadata...")
    # read new samples
    new_samples = pl.read_ndjson(NEW_SAMPLES)

    if OLD_DATE:
        logger.info(f"Merging with old metadata from {OLD_DATE}")
        OLD_SAMPLES = JSON.replace(NEW_DATE, OLD_DATE)
        try:
            old_samples = pl.read_ndjson(OLD_SAMPLES)
            new_samples = pl.concat([old_samples, new_samples]).unique('acc')
        except FileNotFoundError:
            logger.error(f"No old metadata found for {OLD_SAMPLES}")
            logger.error("""
Check dates of existing metadata and the `outputs/last_run` log.
If these don't match something didn't work last run.

Update the `outputs/last_run` log to match the existing metadata files.
Or remove the existing metadata files + `outputs/last_run` to trigger a fresh start.
            """)
            raise FileNotFoundError(f"No old metadata found for {OLD_SAMPLES}")

    else:
        logger.info("No old metadata found. Writing new metadata...")

    # write
    new_samples.write_ndjson(JSON)
    logger.info(f"Updated metadata: {JSON}")

    if OLD_DATE:
        logger.info(f"Removing old metadata file: {OLD_SAMPLES}")
        os.remove(OLD_SAMPLES)


if __name__ == '__main__':
    main()