# log_config.py
import logging

def configure_logger(name):
    logger = logging.getLogger(name)

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s [%(name)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    return logger