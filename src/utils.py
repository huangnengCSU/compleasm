import os
import tarfile
import time
import urllib.request
import hashlib
from subprocess import check_output, PIPE, Popen
from sys import stderr
import logging
import time
class Error(Exception):
    """
    Module-specific exception
    """

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return self.value

def subprocess_popen(args, stdin=None, stdout=PIPE, stderr=stderr, bufsize=8388608):
    return Popen(args, stdin=stdin, stdout=stdout, stderr=stderr, bufsize=bufsize, universal_newlines=True)

# def init_logger(log_file=None):
#     log_format = logging.Formatter("[%(asctime)s %(levelname)s] %(message)s")
#     logger = logging.getLogger()
#     logger.setLevel(logging.INFO)
#
#     console_handler = logging.StreamHandler()
#     console_handler.setFormatter(log_format)
#     logger.handlers = [console_handler]
#
#     if log_file and log_file != '':
#         file_handler = logging.FileHandler(log_file)
#         file_handler.setFormatter(log_format)
#         logger.addHandler(file_handler)
#     return logger

class MinibuscoLogger():
    def __init__(self, logger=None):
        self.logger = logging.getLogger(logger)
        self.logger.setLevel(logging.DEBUG)

        self.log_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        self.pid = str(os.getpid())
        self.log_path = "logs"
        self.log_name = logger.split('.')[-1]+'.log'

        if not os.path.exists(self.log_path):
            os.mkdir(self.log_path)

        fh = logging.FileHandler(os.path.join(self.log_path, self.log_name), mode='a', encoding='utf-8')
        fh.setLevel(logging.DEBUG)

        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)

        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        fh.setFormatter(formatter)
        ch.setFormatter(formatter)

        self.logger.addHandler(fh)
        self.logger.addHandler(ch)

        fh.close()
        ch.close()

    def getlog(self):
        return self.logger
