import logging
from logging import Logger

def get_logger(name:str) -> Logger:
    """
    Create a Logger instance with the specified name, and two handlers for INFO
    and ERROR messages with different formatters

    Args:
        name (str): Name of the logger, e.g. the module from where it is called

    Returns:
        Logger
    """
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)

    # Handler for INFO level messages
    info_handler = logging.StreamHandler()
    info_handler.setLevel(logging.INFO)
    info_format = logging.Formatter('%(message)s')
    info_handler.setFormatter(info_format)
    # Create filter to only handle INFO messages
    class InfoFilter(logging.Filter):
        def filter(self, record):
            return record.levelname == 'INFO'
    f = InfoFilter()
    info_handler.addFilter(f)

    # Handler for ERROR level messages
    error_handler = logging.StreamHandler()
    error_handler.setLevel(logging.ERROR)
    error_format = logging.Formatter('%(levelname)s - %(name)s - %(message)s')
    error_handler.setFormatter(error_format)

    # Add handlers to logger
    logger.addHandler(info_handler)
    logger.addHandler(error_handler)

    return logger