import os
import hashlib
import tarfile
import shutil
import urllib.request
from .utils import Error, MinibuscoLogger

logger = MinibuscoLogger().getlog(__name__)


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


class Downloader:
    def __init__(self, download_dir=None):
        self.base_url = "https://busco-data.ezlab.org/v5/data/"
        self.default_lineage = ["archaea_odb10", "bacteria_odb10", "eukaryota_odb10"]
        if download_dir is None:
            self.download_dir = "downloads"
        else:
            self.download_dir = download_dir
        if not os.path.exists(self.download_dir):
            os.mkdir(self.download_dir)
        self.placement_dir = os.path.join(self.download_dir, "placement_files")
        if not os.path.exists(self.placement_dir):
            os.mkdir(self.placement_dir)
        self.lineage_description, self.placement_description = self.download_file_version_document()
        logger.info("Get file version description done.")
        self.download_placement()  # download placement files
        for lineage in self.default_lineage:
            try:
                if self.check_lineage(lineage):
                    pass
                    # self.lineage_description[lineage].append(os.path.join(self.download_dir, lineage))
                else:
                    self.download_lineage(lineage)
            except KeyError:
                raise Error("invalid lineage name: {}".format(lineage))

    @staticmethod
    def download_single_file(remote_filepath, local_filepath, expected_hash):
        try:
            urllib.request.urlretrieve(remote_filepath, local_filepath)
            observed_hash = md5(local_filepath)
            if observed_hash != expected_hash:
                logger.info("md5 hash is incorrect: {} while {} expected".format(str(observed_hash), str(expected_hash)))
                logger.info("deleting corrupted file {}".format(local_filepath))
                # os.remove(local_filepath)
                logger.error("Unable to download necessary files")
                raise Error("Unable to download necessary files")
            else:
                logger.info("Success download from {}".format(remote_filepath))
        except URLError:
            logger.error("Cannot reach {}".format(remote_filepath))
            return False
        return True

    def download_file_version_document(self):
        file_version_url = self.base_url + "file_versions.tsv"
        file_version_download_path = os.path.join(self.download_dir, "file_versions.tsv")
        hash_url = self.base_url + "file_versions.tsv.hash"
        hash_download_path = os.path.join(self.download_dir, "file_versions.tsv.hash")

        if not os.path.exists(hash_download_path):
            # download hash file
            try:
                urllib.request.urlretrieve(hash_url, hash_download_path)
            except URLError:
                logger.error("Cannot reach {}".format(hash_url))
                raise Error("Unable to download necessary files")
        expected_file_version_hash = ""
        with open(hash_download_path, 'r') as fin:
            expected_file_version_hash = fin.readline().strip()

        if os.path.exists(file_version_download_path):
            download_success = True
        else:
            # download file version
            download_success = self.download_single_file(file_version_url, file_version_download_path,
                                                         expected_file_version_hash)
        lineages_description_dict = {}
        placement_description_dict = {}
        if download_success:
            with open(file_version_download_path, 'r') as fin:
                for line in fin:
                    strain, date, hash_value, category, info = line.strip().split()
                    if info == "lineages":
                        lineages_description_dict[strain] = [date, hash_value, category]
                    elif info == "placement_files":
                        placement_description_dict[strain] = [date, hash_value, category]
            return lineages_description_dict, placement_description_dict
        else:
            return None, None

    def check_lineage(self, lineage):
        try:
            if os.path.exists(os.path.join(self.download_dir, lineage)):
                self.lineage_description[lineage].append(os.path.join(self.download_dir, lineage))
                return True
            else:
                return False
        except KeyError:
            raise Error("invalid lineage name: {}".format(lineage))

    def download_lineage(self, lineage=None):
        if lineage is None:
            lineages = self.default_lineage
        else:
            lineages = [lineage]
        for lineage in lineages:
            if self.check_lineage(lineage):
                continue
            date, expected_hash = self.lineage_description[lineage][0:2]  # [date, hash_value, category]
            remote_url = self.base_url + "lineages/{}.{}.tar.gz".format(lineage, date)
            download_path = os.path.join(self.download_dir, "{}.{}.tar.gz".format(lineage, date))
            download_success = self.download_single_file(remote_url, download_path, expected_hash)
            if download_success:
                tar = tarfile.open(download_path)
                tar.extractall(self.download_dir)
                tar.close()
                logger.info("Lineage file extraction path: {}/{}".format(self.download_dir, lineage))
                local_lineage_dir = os.path.join(self.download_dir, lineage)
                self.lineage_description[lineage].append(local_lineage_dir)

    def download_placement(self):
        for strain in self.placement_description.keys():
            date, expected_hash, category = self.placement_description[strain]
            if strain.startswith("supermatrix"):
                prefix, aln, version, sufix = strain.split(".")
                download_file_name = "{}.{}.{}.{}.{}.tar.gz".format(prefix, aln, version, date, sufix)
            else:
                prefix, version, sufix = strain.split(".")
                download_file_name = "{}.{}.{}.{}.tar.gz".format(prefix, version, date, sufix)
            if os.path.exists(os.path.join(self.placement_dir, download_file_name)):
                self.placement_description[strain].append(
                    os.path.join(self.placement_dir, download_file_name.replace(".tar.gz", "")))
                continue
            else:
                remote_url = self.base_url + "placement_files/{}".format(download_file_name)
                download_path = os.path.join(self.placement_dir, download_file_name)
                download_success = self.download_single_file(remote_url, download_path, expected_hash)
                if download_success:
                    tar = tarfile.open(download_path)
                    tar.extractall(self.placement_dir)
                    tar.close()
                    logger.info("Placement file extraction path: {}/{}".format(self.placement_dir,
                                                                         download_file_name.replace(".tar.gz", "")))
                    self.placement_description[strain].append(
                        os.path.join(self.placement_dir, download_file_name.replace(".tar.gz", "")))


if __name__ == "__main__":
    d = Downloader("downloads")
    d.download_lineage("thiotrichales_odb10")
    d.download_lineage("saccharomycetes_odb10")
    d.download_lineage("diptera_odb10")
    # d.download_lineage("rhodospirillales_odb10")
    # d.download_lineage("cheoctovirus_odb10")
    # d.download_lineage("plasmodium_odb10")
    # for k in d.lineage_description:
    #     print(k)
    #     print(d.lineage_description[k])

    for k in d.placement_description:
        print(k)
        print(d.placement_description[k])
