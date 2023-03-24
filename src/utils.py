import os
import tarfile
import urllib.request
import hashlib
class Error(Exception):
    """
    Module-specific exception
    """

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return self.value