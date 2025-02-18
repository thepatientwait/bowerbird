#!/usr/bin/env python3

###############################################################################
#
#    Copyright (C) 2020 Ben Woodcroft
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

__author__ = "Ben Woodcroft"
__copyright__ = "Copyright 2020"
__credits__ = ["Ben Woodcroft"]
__license__ = "GPL3"
__maintainer__ = "Ben Woodcroft"
__email__ = "benjwoodcroft near gmail.com"
__status__ = "Development"

import argparse
import logging
import sys
import os

import extern

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser()

    parent_parser.add_argument('--date',required=True, help='date of run e.g. 20111130')
    
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")
    parent_parser.add_argument('--target', help='amplicon or shotgun', choices=['amplicon','shotgun'], required=True)
    parent_parser.add_argument('--bucket-name', help='bucket name', default='singlem-sra-metadata-gathering')

    parent_parser.add_argument('--start-date', help='start date for query', default=None, required=False)
    parent_parser.add_argument('--dry-run', help='do not run the query', action="store_true")

    args = parent_parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    logging.info("Starting query_sra_for_metadata.py")
    
    # define metagenome taxonomy table
    sql_pretable = """
    -- Loop counter
    DECLARE counter int64 default 1;

    -- Create intermediate table with metagenome taxonomy entries
    CREATE OR REPLACE TABLE 
        metadata.emp AS
    WITH 
        cte AS (
            SELECT 
                1 AS xlevel, tax_id, parent_id, sci_name
            FROM 
                `nih-sra-datastore.sra_tax_analysis_tool.taxonomy` e
            WHERE 
                e.sci_name = 'metagenomes'
            )
    SELECT *
    FROM cte;

    WHILE EXISTS (
        SELECT c.*
        FROM 
            metadata.emp p
        INNER JOIN 
            `nih-sra-datastore.sra_tax_analysis_tool.taxonomy` c 
            ON p.tax_id = c.parent_id
        WHERE 
            p.xlevel = counter
        )
    DO
    -- Insert next level
    INSERT INTO 
        metadata.emp ( xlevel, tax_id, parent_id, sci_name )
    SELECT 
        counter + 1 AS xlevel, c.tax_id, c.parent_id, c.sci_name
    FROM 
        metadata.emp p
    INNER JOIN 
        `nih-sra-datastore.sra_tax_analysis_tool.taxonomy` c ON p.tax_id = c.parent_id
    WHERE 
        p.xlevel = counter;

    SET counter = counter + 1;

    -- Loop safely
    IF counter > 10
    THEN
        BREAK;
    END IF;
    END WHILE;
    """

    sql = """
    EXPORT DATA OPTIONS(
        uri='__SQL_JSON_RESULTS_URI__',
        format='JSON') AS

    SELECT
        -- count(*)
        acc,
        mbases,
        mbytes,
        assay_type,
        releasedate,
        organism,
        libraryselection,
        librarysource,
        biosample,
        bioproject,
        sra_study,
        avgspotlen,
        attributes
    FROM
        `nih-sra-datastore.sra.metadata`
    WHERE
        acc IN (
        SELECT
            acc
        FROM
            `nih-sra-datastore.sra.metadata`
        WHERE
            (librarysource = 'METAGENOMIC'
                OR organism IN (SELECT sci_name FROM metadata.emp))
            AND librarysource NOT IN ('VIRAL RNA', 'METATRANSCRIPTOMIC', 'TRANSCRIPTOMIC')
            AND platform = 'ILLUMINA'
            AND consent = 'public'"""

    # Add date filter
    if args.start_date:
        sql += f"""
            AND releasedate >= '{args.start_date}'"""

    # Add target filter
    if args.target == 'shotgun':
        sql += """
            AND (mbases > 1000
                OR (libraryselection = 'RANDOM'
                    AND mbases > 100))
            AND mbases <= 200000)
    ORDER BY releasedate ASC
    LIMIT 10
            ;"""

    elif args.target == 'amplicon':
        sql += """
            AND libraryselection = 'PCR'
            AND mbases <= 2000)
            ;"""

    else: raise ValueError("Unknown target: {}".format(args.target))

    # Now do the full query
    destination_path = f'gs://{args.bucket_name}/{args.target}_sra_{args.date}/*'
    # TODO: combine with initial query
    sql = sql.replace('__SQL_JSON_RESULTS_URI__',destination_path)

    

    if args.dry_run:
        logging.info("Dry run: not running query")
        logging.info(f"Would have written to {destination_path}")
        logging.info("Would have run the following SQL:")
        logging.info(sql)
        logging.info("NOTE: this did not test the SQL, it may not work as expected.")
        exit(0)

    ## CREATE METADATA DATASET
    logging.info("Creating metadata dataset...")
    extern.run('bq mk -f -d metadata')

    ## CREATE TEMPORARY TABLE
    logging.info("Creating temporary taxonomy table...")
    extern.run('bq query --use_legacy_sql=false', stdin=sql_pretable)

    ## FULL QUERY
    logging.info(f"Querying and writing {args.target} results to {destination_path}...")
    extern.run('bq query --use_legacy_sql=false', stdin=sql)

    # TODO: combine with initial query
    # Now write the taxonomy table
    destination_path = 'gs://{}/sra_taxonomy_table_{}/*'.format(args.bucket_name, args.date)
    logging.info("Querying and writing taxonomy table results to {}...".format(destination_path))

    sql = "EXPORT DATA OPTIONS(uri='{}', format='JSON') AS select * from metadata.emp".format(destination_path)

    extern.run('bq query --use_legacy_sql=false', stdin=sql)
    logging.info("Finished querying for the taxonomy table.")