# init
# ----
from google.cloud import bigquery

import polars as pl

# logging
# -------
from log_config import configure_logger
import os

logger_name = os.path.splitext(os.path.basename(__file__))[0]
logger      = configure_logger(logger_name)

# snakefile parameters
# --------------------
QUERY                   = snakemake.input['query']
JSON                    = snakemake.output['json']

if QUERY == 'sql_queries/metadata_query.sql':
    NEW_DATE                = snakemake.params['new_date']
    OLD_DATE                = snakemake.params['old_date']
    TARGET                  = snakemake.params['target']



# define sql queries
# ------------------

with open(QUERY, 'r') as file:
    sql = file.read()

if QUERY == 'sql_queries/metadata_query.sql':

    start_date_filter = ""

    if NEW_DATE is not None:
        start_date_filter = f"""
            AND releasedate < '{NEW_DATE}'"""

    if OLD_DATE is not None:
        start_date_filter += f"""
            AND releasedate >= '{OLD_DATE}'"""    

    if TARGET == 'shotgun':
        target_filter = """
        AND platform = 'ILLUMINA'
        AND (mbases > 1000
            OR (libraryselection = 'RANDOM'
                AND mbases > 100))
        AND mbases <= 200000"""
    elif TARGET == 'amplicon':
        target_filter = """
        AND platform = 'ILLUMINA'
        AND libraryselection = 'PCR'
        AND mbases <= 6000"""
    elif TARGET == 'nanopore':
        target_filter = """
        AND platform = 'OXFORD_NANOPORE'"""
    else:
        raise ValueError("Unknown target: {}".format(TARGET))

    sql = sql.format(start_date_filter=start_date_filter, target_filter=target_filter)

# run query
# ---------
client      = bigquery.Client()

logger.info(f"Running query...")
query_job   = client.query(sql)
results     = query_job.result().to_dataframe()

# write results to ndjson file
logger.info(f"Writing results to json file: {JSON}...")
pl.from_pandas(results).write_ndjson(JSON)