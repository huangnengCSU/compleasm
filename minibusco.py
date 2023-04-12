#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Neng Huang
Email: neng@ds.dfci.harvard.edu
Date: 2023 Apr 11
"""

import os
import argparse
import hashlib
import sys
import tarfile
import urllib.request


### utils
class Error(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return self.value


### DownloadLineage.py

class URLError(OSError):
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
    def __init__(self, download_dir=None, download_lineage=True, download_placement=True):
        self.base_url = "https://busco-data.ezlab.org/v5/data/"
        self.default_lineage = ["eukaryota_odb10"]
        if download_dir is None:
            self.download_dir = "mb_downloads"
        else:
            self.download_dir = download_dir

        self.placement_dir = os.path.join(self.download_dir, "placement_files")

        if not os.path.exists(self.download_dir):
            os.mkdir(self.download_dir)

        self.lineage_description, self.placement_description = self.download_file_version_document()

        if download_placement:
            if not os.path.exists(self.placement_dir):
                os.mkdir(self.placement_dir)
            self.download_placement()  # download placement files

        if download_lineage:
            for lineage in self.default_lineage:
                self.download_lineage(lineage)

    def download_single_file(self, remote_filepath, local_filepath, expected_hash):
        try:
            urllib.request.urlretrieve(remote_filepath, local_filepath)
            observed_hash = md5(local_filepath)
            if observed_hash != expected_hash:
                print("md5 hash is incorrect: {} while {} expected".format(str(observed_hash), str(expected_hash)))
                print("deleting corrupted file {}".format(local_filepath))
                # os.remove(local_filepath)
                print("Unable to download necessary files")
                raise Error("Unable to download necessary files")
            else:
                print("Success download from {}".format(remote_filepath))
        except URLError:
            print("Cannot reach {}".format(remote_filepath))
            return False
        return True

    def download_file_version_document(self):
        file_version_url = self.base_url + "file_versions.tsv"
        file_version_download_path = os.path.join(self.download_dir, "file_versions.tsv")
        hash_url = self.base_url + "file_versions.tsv.hash"
        hash_download_path = os.path.join(self.download_dir, "file_versions.tsv.hash")

        if os.path.exists(file_version_download_path + ".tmp"):
            sys.exit("file_versions.tsv.tmp exists, another process is downloading, please run again later.")

        if not os.path.exists(file_version_download_path + ".done"):
            open(file_version_download_path + ".tmp", 'w').close()

            # download hash file
            try:
                urllib.request.urlretrieve(hash_url, hash_download_path)
            except URLError:
                print("Cannot reach {}".format(hash_url))
                os.remove(file_version_download_path + ".tmp")
                raise Error("Unable to download necessary file {}".format("file_versions.tsv.hash"))

            with open(hash_download_path, 'r') as fin:
                expected_file_version_hash = fin.readline().strip()

            download_success = self.download_single_file(file_version_url,
                                                         file_version_download_path,
                                                         expected_file_version_hash)
            if not download_success:
                os.remove(file_version_download_path + ".tmp")
                raise Error("Unable to download necessary file {}".format("file_versions.tsv"))
            else:
                open(file_version_download_path + ".done", 'w').close()
                os.remove(file_version_download_path + ".tmp")

        lineages_description_dict = {}
        placement_description_dict = {}
        with open(file_version_download_path, 'r') as fin:
            for line in fin:
                strain, date, hash_value, category, info = line.strip().split()
                if info == "lineages":
                    lineages_description_dict[strain] = [date, hash_value, category]
                elif info == "placement_files":
                    placement_description_dict[strain] = [date, hash_value, category]
        return lineages_description_dict, placement_description_dict

    def download_lineage(self, lineage):
        if not lineage.endswith("_odb10"):
            lineage = lineage + "_odb10"
        if os.path.exists(os.path.join(self.download_dir, lineage) + ".tmp"):
            sys.exit("{}.tmp exists, another process is downloading, please run again later.".format(lineage))

        if not os.path.exists(os.path.join(self.download_dir, lineage) + ".done"):
            open(os.path.join(self.download_dir, lineage) + ".tmp", 'w').close()

            try:
                date, expected_hash = self.lineage_description[lineage][0:2]  # [date, hash_value, category]
            except KeyError:
                os.remove(os.path.join(self.download_dir, lineage) + ".tmp")
                raise Error("invalid lineage name: {}".format(lineage))

            remote_url = self.base_url + "lineages/{}.{}.tar.gz".format(lineage, date)
            download_path = os.path.join(self.download_dir, "{}.{}.tar.gz".format(lineage, date))
            download_success = self.download_single_file(remote_url, download_path, expected_hash)

            if not download_success:
                os.remove(os.path.join(self.download_dir, lineage) + ".tmp")
                raise Error("Unable to download necessary file {}".format("{}.{}.tar.gz".format(lineage, date)))

            if download_success:
                tar = tarfile.open(download_path)
                # tar.extractall(self.download_dir)
                try:
                    tar.extractall(self.download_dir, members=[tar.getmember('{}/refseq_db.faa.gz'.format(lineage)),
                                                               tar.getmember('{}/links_to_ODB10.txt'.format(lineage))])
                except KeyError:
                    tar.extractall(self.download_dir, members=[tar.getmember('{}/refseq_db.faa.gz'.format(lineage))])
                except:
                    os.remove(os.path.join(self.download_dir, lineage) + ".tmp")
                    raise Error("Unable to extract file: {} or {} from {}".format("refseq_db.faa.gz",
                                                                                  "links_to_ODB10.txt",
                                                                                  download_path))

                tar.close()
                print("Lineage file extraction path: {}/{}".format(self.download_dir, lineage))
                local_lineage_dir = os.path.join(self.download_dir, lineage)
                self.lineage_description[lineage].append(local_lineage_dir)
                open(os.path.join(self.download_dir, lineage) + ".done", 'w').close()
                os.remove(os.path.join(self.download_dir, lineage) + ".tmp")

    def download_placement(self):
        if os.path.exists(self.placement_dir + ".tmp"):
            sys.exit("placement_files.tmp exists, another process is downloading, please run again later.")

        if not os.path.exists(self.placement_dir + ".done"):
            open(self.placement_dir + ".tmp", 'w').close()

            for strain in self.placement_description.keys():
                date, expected_hash, category = self.placement_description[strain]
                if strain.startswith("supermatrix"):
                    prefix, aln, version, sufix = strain.split(".")
                    download_file_name = "{}.{}.{}.{}.{}.tar.gz".format(prefix, aln, version, date, sufix)
                else:
                    prefix, version, sufix = strain.split(".")
                    download_file_name = "{}.{}.{}.{}.tar.gz".format(prefix, version, date, sufix)

                if "eukaryota" not in download_file_name:
                    continue

                remote_url = self.base_url + "placement_files/{}".format(download_file_name)
                download_path = os.path.join(self.placement_dir, download_file_name)
                download_success = self.download_single_file(remote_url, download_path, expected_hash)

                if not download_success:
                    os.remove(self.placement_dir + ".tmp")
                    raise Error("Unable to download necessary file {}".format(download_file_name))

                if download_success:
                    tar = tarfile.open(download_path)
                    tar.extractall(self.placement_dir)
                    tar.close()
                    print("Placement file extraction path: {}/{}".format(self.placement_dir,
                                                                         download_file_name.replace(".tar.gz", "")))
                    self.placement_description[strain].append(
                        os.path.join(self.placement_dir, download_file_name.replace(".tar.gz", "")))
            open(self.placement_dir + ".done", 'w').close()
            os.remove(self.placement_dir + ".tmp")


def download(args):
    if args.destination is not None:
        downloader = Downloader(args.destination)
    else:
        downloader = Downloader()
    downloader.download_lineage(args.lineage)


def list_lineages(args):
    if not args.local and not args.remote:
        sys.exit("\n Usage error: Please specify whether to list local or remote lineages."
                 "\n e.g. minibusco.py list --remote or minibusco.py list --local --library_path /path/to/lineage_folder\n")
    if args.local:
        if args.library_path is None:
            sys.exit("\n Usage error: Please specify the folder path to stored lineages."
                     "\n e.g. minibusco list --local --library_path /path/to/lineages_folder\n")
        else:
            print("Local available lineages:")
            for file in os.path.listdir(args.library_path):
                if file.endswith("_odb10.done"):
                    print(file.replace("_odb10.done", ""))
    if args.remote:
        if args.library_path is not None:
            downloader = Downloader(args.library_path, download_lineage=False, download_placement=False)
        else:
            downloader = Downloader(download_lineage=False, download_placement=False)
        print("Remote available lineages:")
        for lineage in downloader.lineage_description.keys():
            print(lineage)


### main.py
def main():
    parser = argparse.ArgumentParser(description="MiniBusco")
    subparser = parser.add_subparsers(dest="command", help="Minibusco modules help", required=True)

    ### sub-command: download
    download_parser = subparser.add_parser("download", help="Download BUSCO lineage")
    download_parser = subparser.add_parser("download", help="Download specified BUSCO lineage")
    download_parser.add_argument("-l", "--lineage", type=str,
                                 help="Specify the name of the BUSCO lineage to be downloaded. (e.g. eukaryota, primates, saccharomycetes etc.)",
                                 required=True)
    download_parser.add_argument("-d", "--destination", type=str, help="Folder path to download folder", default=None)
    download_parser.set_defaults(func=download)

    ### sub-command: list
    list_parser = subparser.add_parser("list", help="List local or remote BUSCO lineages")
    list_parser.add_argument("--remote", action="store_true", help="List remote BUSCO lineages")
    list_parser.add_argument("--local", action="store_true", help="List local BUSCO lineages")
    list_parser.add_argument("--library_path", type=str, help="Folder path to stored lineages. ", default=None)
    list_parser.set_defaults(func=list_lineages)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
