"""

"""

import os
import logging
import traceback


def get_logger(filepath, level=logging.DEBUG):
    # create logger
    logger = logging.getLogger(os.path.basename(filepath))
    logger.setLevel(level)
    # create console handler and set level to debug
    ch = logging.StreamHandler()
    ch.setLevel(level)
    # create formatter
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    # add formatter to ch
    ch.setFormatter(formatter)
    # add ch to logger
    logger.addHandler(ch)

    return logger


def log_error():
    """Nextflow specific function that logs an error upon unexpected failing
    """

    with open(".status", "w") as status_fh:
        status_fh.write("error")


class MainWrapper:

    def __init__(self, f):

        self.f = f

    def __call__(self, *args, **kwargs):

        context = self.f.__globals__

        logger = context.get("logger", None)
        build_versions = context.get("build_versions", None)

        try:
            if build_versions:
                build_versions()
            self.f(*args, **kwargs)
        except:
            if logger:
                logger.error("Module exited unexpectedly with error:"
                             "\\n{}".format(traceback.format_exc()))
            log_error()
