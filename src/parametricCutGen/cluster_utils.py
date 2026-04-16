import logging
import os
import tomllib
from pyscipopt import Model, SCIP_EVENTTYPE


def parse_logger(name):
    """Set logger information from cluster paramters files."""
    logger = logging.getLogger(name)
    logging_level = os.getenv("LOGGING_LEVEL")
    if logging_level == "debug":
        logger.info(f"Adjusting the default logging level to DEBUG")    
        logger.setLevel(logging.DEBUG)
    elif logging_level == "warning":
        logger.info(f"Adjusting the default logging level to WARNING")
        logger.setLevel(logging.WARNING)
    elif logging_level == "error":
        logger.info(f"Adjusting the default logging level to ERROR")
        logger.setLevel(logging.ERROR)
    else:
        logger.setLevel(logging.INFO)
    logger.info(f"Logging level set to {initial_gen_logger.level}")
    return logger

class validate_paths(dict, logger):
    """Checks that specified dictionary with keys "name_of_path" and values path"""
    def __call__(self):
        paths = self
        paths_are_valid = True
        bad_paths = {}
        for name_of_path in paths:
            if os.path.exists(paths[name_of_path]):
                logger.debug(f"Path {name_of_path} with dir {paths[name_of_path]} exists.")
            else:
                logger.debug(f"Path {name_of_path} with dir {paths[name_of_path]} does not exists.")
                paths_are_valid = False
                bad_paths[path] = paths[path]    
        if paths_are_valid is False:
            logger.error(f"All paths should exists. We have the following bad paths: {bad_paths}")
        return paths

