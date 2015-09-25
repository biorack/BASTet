"""
Module providing functionality for logging based on the python logging module.
The module is intended toease the use of logging while a developer
can still access the standard python logging mechanism if needed.
"""
import logging


class log_helper(object):
    """
    BASTet helper module to ease the use of logging

    Class Variables:

        * ``log_levels`` : Dictionary describing the different available logging levels.
    """
    initialized = False

    log_levels = {'CRITICAL': logging.CRITICAL,
                  'ERROR': logging.ERROR,
                  'WARNING': logging.WARNING,
                  'INFO': logging.INFO,
                  'DEBUG': logging.DEBUG,
                  'NOTSET': logging.NOTSET}

    global_log_level = log_levels['WARNING']

    @classmethod
    def setup_logging(cls, level=log_levels['WARNING']):
        """
        Call this function at the beginning of your code to initiate logging.

        Parameters:

            * ``level`` : The default log level to be used. One of ``log_helper.log_level``.

        """
        log_helper.global_log_level = level
        logging.basicConfig(level=level, format=cls.get_default_format())

    @classmethod
    def get_default_format(cls):
        """
        Get default formatting string.
        """
        return '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

    @classmethod
    def set_log_level(cls, level):
        """
        Set the logging level for all BASTet loggers

        Parameters:

            * ``level`` : The logging levels to be used, one of the values specified in log_helper.log_levels.
        """
        if level not in log_helper.log_levels.values():
            try:
                level = log_helper.log_levels[level]
            except KeyError:
                raise KeyError('Invalid log level given')
        log_helper.global_log_level = level
        ld = logging.Logger.manager.loggerDict
        for k in ld.keys():
            if k.startswith('omsi.'):
                cls.get_logger(k).setLevel(level)

    @classmethod
    def get_logger(cls, module_name):
        """
        Get the logger for a particular module. The module_name
        should always be set to the __name__ variable of the calling module.

        Parameters:

            * ``module_name`` : __name__ of the calling module or None in case the ROOT logger should be used.

        Returns:

            * Python logging.Logger retrieved via logging.getLogger.

        """
        if module_name is not None:
            logobj = logging.getLogger(module_name)
            # Make sure we set the correct logging level if we have not created the logger before
            logobj.setLevel(log_helper.global_log_level)
            return logobj
        else:
            return logging.getLogger()

    @classmethod
    def debug(cls, module_name, message, *args, **kwargs):
        """
        Create a debug log entry. This function is typically called as:

        log_helper.debug(module_name=__name__, message="your message")

        Parameters:

            * ``module_name`` : __name__ of the calling module or None in case the ROOT logger should be used.
            * ``args`` : Additional positional arguments for the python logger.debug function. See the python docs.
            * ``kwargs`` : Additional keyword arguments for the python logger.debug function. See the python docs.

        """
        cls.get_logger(module_name).debug(message, *args, **kwargs)

    @classmethod
    def info(cls, module_name, message, *args, **kwargs):
        """
        Create a info log entry. This function is typically called as:

        log_helper.info(module_name=__name__, message="your message")

        Parameters:

            * ``module_name`` : __name__ of the calling module or None in case the ROOT logger should be used.
            * ``args`` : Additional positional arguments for the python logger.info function. See the python docs.
            * ``kwargs`` : Additional keyword arguments for the python logger.info function. See the python docs.

        """
        cls.get_logger(module_name).info(message, *args, **kwargs)

    @classmethod
    def warning(cls, module_name, message, *args, **kwargs):
        """
        Create a warning log entry. This function is typically called as:

        log_helper.warning(module_name=__name__, message="your message")

        Parameters:

            * ``module_name`` : __name__ of the calling module or None in case the ROOT logger should be used.
            * ``args`` : Additional positional arguments for the python logger.warning function. See the python docs.
            * ``kwargs`` : Additional keyword arguments for the python logger.warning function. See the python docs.

        """
        cls.get_logger(module_name).warning(message, *args, **kwargs)

    @classmethod
    def error(cls, module_name, message, *args, **kwargs):
        """
        Create a error log entry. This function is typically called as:

        log_helper.error(module_name=__name__, message="your message")

        Parameters:

            * ``module_name`` : __name__ of the calling module or None in case the ROOT logger should be used.
            * ``args`` : Additional positional arguments for the python logger.error function. See the python docs.
            * ``kwargs`` : Additional keyword arguments for the python logger.error function. See the python docs.

        """
        cls.get_logger(module_name).error(message, *args, **kwargs)

    @classmethod
    def critical(cls, module_name, message, *args, **kwargs):
        """
        Create a critical log entry. This function is typically called as:

        log_helper.critical(module_name=__name__, message="your message")

        Parameters:

            * ``module_name`` : __name__ of the calling module or None in case the ROOT logger should be used.
            * ``args`` : Additional positional arguments for the python logger.critical function. See the python docs.
            * ``kwargs`` : Additional keyword arguments for the python logger.critical function. See the python docs.

        """
        cls.get_logger(module_name).critical(message, *args, **kwargs)

    @classmethod
    def exception(cls, module_name, message, *args, **kwargs):
        """
        Create a exception log entry. This function is typically called as:

        log_helper.exception(module_name=__name__, message="your message")

        Parameters:

            * ``module_name`` : __name__ of the calling module or None in case the ROOT logger should be used.
            * ``args`` : Additional positional arguments for the python logger.exception function. See the python docs.
            * ``kwargs`` : Additional keyword arguments for the python logger.exception function. See the python docs.

        """
        cls.get_logger(module_name).exception(message, *args, **kwargs)


if not log_helper.initialized:
    log_helper.setup_logging()
    log_helper.set_log_level(log_helper.global_log_level)
    log_helper.info(__name__, 'Initialized logging')
