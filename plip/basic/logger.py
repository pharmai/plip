import inspect
import logging


def get_logger():
    """
    Configures a base logger and returns a module-specific sub-logger of the calling module.
    """
    frame = inspect.stack()[1]
    module_name = inspect.getmodule(frame[0]).__name__
    if module_name != '__main__':
        logger = logging.getLogger(module_name)
        if not logger.parent.handlers:
            ch = logging.StreamHandler()
            formatter = logging.Formatter('%(asctime)s [%(levelname)s] [%(filename)s:%(lineno)d] %(name)s: %(message)s')
            ch.setFormatter(formatter)
            logger.parent.addHandler(ch)
    else:
        logger = logging.getLogger('plip')

    return logger
