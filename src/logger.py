import logging

def setup_logger(name):
    # Set up Logger and log starting message:
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    # main handler for logging
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    
    # add the handlers to the logger
    ch.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s', datefmt = "%H:%M:%S"))
    logger.addHandler(ch)

    return logger

# TODO Fix this eventually
def set_log_level(logger_name: str, level: str):
    logger = logging.getLogger(logger_name)
    match level:
        case 'info':
            logger.setLevel(logging.INFO)
        case 'debug':
            logger.setLevel(logging.DEBUG)


# TODO should I do this? Does it matter? Useful in jupyter notebooks to avoid adding duplicate handlers
def clear_logger(logger):
    logger.handlers.clear()