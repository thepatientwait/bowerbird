import os
import gzip
import logging
import sys

def get_logger():
    """
    Dynamically get the logger of the calling script.
    """
    # Get the name of the calling module
    script_name = os.path.splitext(os.path.basename(sys.argv[0]))[0]
    script_name = script_name.split('.')[-1]  # Remove the file extension

    # Return a logger with the calling module's name
    return logging.getLogger(script_name)

def remove_previous_log_files(
    output_dir
    ):
    """
    Remove any previous log files in the output directory.
    """
    logger = get_logger()
    logger.info(f"Removing previous log files in {output_dir}.")
    for file in os.listdir(output_dir):
        if file == 'log' or file == 'failed':
            file_path = os.path.join(output_dir, file)
            logger.info(f"Removing {file_path}.")
            os.remove(file_path)


def check_for_previous_output(output):
    """
    Check if the run_otu_table already exists and is valid.
    If it does, exit the script.
    """
    logger = get_logger()
    if os.path.exists(output) or os.path.islink(output):
        if os.path.islink(output):
            target = os.readlink(output)
            if not os.path.exists(target):
                logger.info(f"Output {output} is a broken symlink. Removing it.")
                os.remove(output)
                return

        if os.path.getsize(output) == 0:
            logger.info(f"Output {output} exists but is empty. Removing it.")
            os.remove(output)
            return

        try:
            with gzip.open(output, 'rt') as f:
                f.read(1)
        except (OSError, gzip.BadGzipFile):
            logger.info(f"Output {output} exists but is corrupted. Removing it.")
            os.remove(output)
            return

        logger.info(f"Output {output} already exists and is valid.")
        logger.info("Exiting without error.")
        exit(0)

def check_for_missing_otu_tables(
    input_otus, 
    output_dir, 
    accession
    ):
    logger = get_logger()

    error_message = None

    if isinstance(input_otus, list) and len(input_otus) > 1:
        n_missing = 0
        for input_otu in input_otus:
            if not os.path.exists(input_otu):
                n_missing += 1

        if n_missing > 0:
            error_message = f"Missing {n_missing} out of {len(input_otus)} OTU tables for {accession}."

    elif isinstance(input_otus, str) or (isinstance(input_otus, list) and len(input_otus) == 1):
        input_otu = input_otus
        if not os.path.exists(input_otu):
            error_message = f"OTU table for {input_otu} is missing."

    if error_message:
        logger.info(error_message)
        with open(os.path.join(output_dir, "log"), "w") as f:
            f.write(error_message)
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

    logger = get_logger()
    if attempt <= total_attempts:
        return

    logger.info("No more attempts! OTU table generation failed.")
    with open(os.path.join(output_dir, "log"), "w") as f:
        f.write(f"No more attempts! OTU table generation failed.")
    with open(os.path.join(output_dir, "failed"), "w") as f:
        f.write(f"{accession}\n")
    logger.info("Logged failed run accession.")
    logger.info("Exiting without error.")
    exit(0)