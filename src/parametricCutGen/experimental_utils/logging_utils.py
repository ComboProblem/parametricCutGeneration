import logging

# Control of logging facilites from igp, spam, parametricCutGen

_alias_logger_names = {"igp_functions": "cutgeneratingfunctionology.igp.functions", "igp_paramateric_real_field": "cutgeneratingfunctionology.igp.parametric", "igp_mfc" :"cutgeneratingfunctionology.igp.minimal_funciton_cell_description"}

 # "pcg_scip_data":  , "pcg_cut_generation_problem": "pcg_cut_score": } # avoid using cut optimization problem because that has an acrynom of cop. ACAB, most likely.

def config_logger(logger_level_default = logging.INFO, **kwds):
    for module_alias in _alias_logger_names:
        try:
            logger = logging.getLogger(_alias_logger_names[module_alias])
            logger.setLevel(kwds[module_alias])
        except KeyError:
            logger = logging.getLogger(_alias_logger_names[module_alias])
            logger.setLevel(logger_level_default)

default_cgf_logging = {"igp_functions": logging.ERROR, "igp_paramateric_real_field": logging.ERROR,  "igp_mfc": logging.ERROR}
default_cgp_logging = {}

config_logger(**default_cgf_logging)

cgp_written_log = logging.FileHandler
