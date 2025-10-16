# Overview

This projects main goal is to regularly update the Sandpiper database.
It uses `snakemake` to manage the workflow and `conda` to ensure a consistent environment.

The workflow is divided into 2 sub-workflows to optimise the speed of snakemake's DAG creation.

## SRA Query
This sub-workflow queries the SRA database to get a list of all the new and updated SRA runs since the last update.
It directly queries the SRA database using Google BigQuery, then downloads the new metadata as a new-line delimited JSON (ndjson) file.

## Sample Analysis
This sub-workflow downloads and processes the new SRA runs.
It also combines runs from the same sample to create a single OTU table per sample.
Both individual runs and combined samples are stored in the Sandpiper database.

1. Download the SRA runs using `kingfisher`.
2. Create OTU tables using `singlem pipe`.
3. Combine OTU tables from runs belonging to the same sample using `singlem summarise`.
4. Assign taxonomy to OTUs using `singlem renew`.

## Main Workflow
The main workflow simply runs the SRA Query then chunks the new runs into manageable batches to be processed by the Sample Analysis workflow.

# Features

Has been designed with extensibility in mind. It currently is just designed to target shotgun metagenome samples, but has been designed to be easily extended to other types of samples (e.g. 16S rRNA gene amplicon sequencing, long-read metagenomics) and metapackages (e.g. `lyrebird`, other more niche targets).

# TODO

[ ] Fix symlink bug.
- Currently samples with just 1 run get a symlink to the outputs of that run.
- The original code used relative symlinks, but the main output dir was moved, breaking thousands of symlinks.
- These all need to be replaced with absolute symlink and the relevant code updated.

[ ] Link with sandpiper database creation workflow.
- This will probably require some changes to the database schema to accommodate multiple types of samples.

[ ] Add tests.

[ ] Create cron job.
