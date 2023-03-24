import os
import hashlib
import urllib.request
from utils import Error


class URLError(OSError):
    # URLError is a sub-type of OSError, but it doesn't share any of
    # the implementation.  need to override __init__ and __str__.
    # It sets self.args for compatibility with other OSError
    # subclasses, but args doesn't have the typical format with errno in
    # slot 0 and strerror in slot 1.  This may be better than nothing.
    def __init__(self, reason, filename=None):
        self.args = reason,
        self.reason = reason
        if filename is not None:
            self.filename = filename

    def __str__(self):
        return '<urlopen error %s>' % self.reason


def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def download_library(remote_filepath, local_filepath, expected_hash):
    try:
        urllib.request.urlretrieve(remote_filepath, local_filepath)
        observed_hash = md5(local_filepath)
        if observed_hash != expected_hash:
            print("md5 hash is incorrect: {} while {} expected".format(str(observed_hash), str(expected_hash)))
            print("deleting corrupted file {}".format(local_filepath))
            os.remove(local_filepath)
            raise Error("Unable to download necessary files")
        else:
            print("md5 hash is {}".format(observed_hash))
    except URLError:
        print("Cannot reach {}".format(remote_filepath))
        return False
    return True


class Downloader:
    def __init__(self):
        self.base_url = "https://busco-data.ezlab.org/v5/data/"
        self.default_lineage = ["archaea", "bacteria", "eukaryota"]

    def download_lineage(self, lineage=None):
        if lineage is None:
            lineages = self.default_lineage
        else:
            lineages = [lineage]
        for lineage in lineages:
            remote_filepath = self.base_url + "lineages/{}_{}.tar.gz"

