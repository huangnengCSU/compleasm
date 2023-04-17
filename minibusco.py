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
import subprocess
import shlex
import shutil
import re
import json
from enum import Enum
from collections import defaultdict
import pandas as pd
import time


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
                    if "{}/refseq_db.faa.gz".format(lineage) not in tar.getnames():
                        os.remove(os.path.join(self.download_dir, lineage) + ".tmp")
                        os.remove(download_path)
                        sys.exit(
                            "No refseq_db.faa.gz in lineage {}, this lineage cannot be used in minibusco! Lineage file has been deleted.".format(
                                lineage))
                    tar.extractall(self.download_dir, members=[tar.getmember('{}/refseq_db.faa.gz'.format(lineage))])
                except:
                    os.remove(os.path.join(self.download_dir, lineage) + ".tmp")
                    os.remove(download_path)
                    raise Error("Unable to extract file: {} or {} from {}".format("refseq_db.faa.gz",
                                                                                  "links_to_ODB10.txt",
                                                                                  download_path))

                tar.close()
                print("Lineage file extraction path: {}/{}".format(self.download_dir, lineage))
                local_lineage_dir = os.path.join(self.download_dir, lineage)
                self.lineage_description[lineage].append(local_lineage_dir)
                open(os.path.join(self.download_dir, lineage) + ".done", 'w').close()
                os.remove(os.path.join(self.download_dir, lineage) + ".tmp")
        else:
            self.lineage_description[lineage].append(os.path.join(self.download_dir, lineage))

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
        else:
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
                self.placement_description[strain].append(
                    os.path.join(self.placement_dir, download_file_name.replace(".tar.gz", "")))


### miniprot ###
def listfiles(folder):
    for root, folders, files in os.walk(folder):
        for filename in folders + files:
            yield os.path.join(root, filename)


class MiniprotRunner:
    def __init__(self, miniprot_execute_command, nthreads=1, autolineage_mode=False):
        if miniprot_execute_command is None:
            miniprot_execute_command = self.search_miniprot()

        print("miniprot execute command:\n {}".format(miniprot_execute_command))
        self.miniprot_execute_command = miniprot_execute_command
        self.threads = nthreads
        self.autolineage = autolineage_mode

    def search_miniprot(self):
        ## Search for miniprot in the path where "minibusco.py" is located
        print("Searching for miniprot in the path where minibusco.py is located")
        script_path = os.path.dirname(os.path.realpath(__file__))
        for fpath in listfiles(script_path):
            path, file = os.path.split(fpath)
            if file == "miniprot" and os.path.isfile(fpath):
                miniprot_execute_command = fpath
                return miniprot_execute_command
        ## Search miniprot in the current execution path
        print("Searching for miniprot in the current execution path")
        execute_path = os.getcwd()
        for fpath in listfiles(execute_path):
            path, file = os.path.split(fpath)
            if file == "miniprot" and os.path.isfile(fpath):
                miniprot_execute_command = fpath
                return miniprot_execute_command
        ## Search for miniprot in PATH
        print("Searching for miniprot in $PATH")
        env_dict = os.environ
        if "PATH" in env_dict:
            path_list = env_dict["PATH"].split(":")
            for path in path_list:
                for fpath in listfiles(path):
                    path, file = os.path.split(fpath)
                    if file == "miniprot" and os.path.isfile(fpath):
                        miniprot_execute_command = fpath
                        return miniprot_execute_command
        sys.exit(
            "miniprot is not found in the path where minibusco.py is located, the current execution path, or PATH. Please check the installation of miniprot.")

    def run_miniprot(self, assembly_filepath, lineage_filepath, alignment_outdir):
        if not os.path.exists(alignment_outdir):
            os.mkdir(alignment_outdir)
        output_filepath = os.path.join(alignment_outdir, "miniprot_output.gff")

        fout = open(output_filepath, "w")
        if self.autolineage:
            ### if autolineage mode, use --aln option to output Ata_seq ###
            miniprot_process = subprocess.Popen(shlex.split(
                "{} -u -I --outs=0.95 -t {} --aln --gff {} {}".format(self.miniprot_execute_command, self.threads,
                                                                      assembly_filepath, lineage_filepath,
                                                                      output_filepath)), stdout=fout, bufsize=8388608)
        else:
            miniprot_process = subprocess.Popen(shlex.split(
                "{} -u -I --outs=0.95 -t {} --gff {} {}".format(self.miniprot_execute_command, self.threads,
                                                                assembly_filepath, lineage_filepath,
                                                                output_filepath)), stdout=fout, bufsize=8388608)
        miniprot_process.wait()
        fout.close()
        return output_filepath


### auto lineage ###
class AutoLineager:
    def __init__(self, sepp_output_directory, sepp_tmp_directory, library_path, threads, sepp_execute_command=None):
        self.sepp_output_folder = sepp_output_directory
        self.sepp_tmp_folder = sepp_tmp_directory
        self.threads = threads
        self.downloader = Downloader(library_path)
        self.lineage_description = self.downloader.lineage_description
        self.placement_description = self.downloader.placement_description
        self.library_folder = self.downloader.download_dir
        self.placement_file_folder = self.downloader.placement_dir
        self.sepp_execute_command = sepp_execute_command

    def run_sepp(self, marker_genes_filapath):
        # select the best one in ["archaea_odb10", "bacteria_odb10", "eukaryota_odb10"] as search_lineage to run repp
        search_lineage = "eukaryota_odb10"
        sepp_output_folder = self.sepp_output_folder
        tmp_file_folder = self.sepp_tmp_folder
        tree_nwk_path = self.placement_description["tree.{}.nwk".format(search_lineage)][3]
        tree_metadata_path = self.placement_description["tree_metadata.{}.txt".format(search_lineage)][3]
        supermaxtix_path = self.placement_description["supermatrix.aln.{}.faa".format(search_lineage)][3]
        print("tree_nwk_path: {}".format(tree_nwk_path))
        print("tree_metadata_path: {}".format(tree_metadata_path))
        print("supermaxtix_path: {}".format(supermaxtix_path))
        if os.path.exists(sepp_output_folder):
            shutil.rmtree(sepp_output_folder)
        if os.path.exists(tmp_file_folder):
            shutil.rmtree(tmp_file_folder)
        sepp_process = "{} --cpu {} --outdir {} -t {} -r {} -a {} -f {} -F 15 -m amino -p {}".format(
            self.sepp_execute_command, self.threads, sepp_output_folder, tree_nwk_path, tree_metadata_path,
            supermaxtix_path, marker_genes_filapath, tmp_file_folder)
        os.system(sepp_process)
        # sepp = subprocess.Popen(sepp_process)
        # sepp.wait()
        return search_lineage

    # Code from https://gitlab.com/ezlab/busco
    def pick_dataset(self, search_lineage):
        # run_folder = self.run_folder

        # load busco dataset name by id in a dict {taxid:name}
        datasets_mapping = {}

        taxid_busco_file_name = "mapping_taxids-busco_dataset_name.{}.txt".format(search_lineage)
        taxid_busco_file_path = self.placement_description[taxid_busco_file_name][3]
        with open(taxid_busco_file_path) as f:
            for line in f:
                parts = line.strip().split("\t")
                tax_id = parts[0]
                dataset = parts[1].split(",")[0]
                datasets_mapping.update({tax_id: dataset})

        # load the lineage for each taxid in a dict {taxid:reversed_lineage}
        # lineage is 1:2:3:4:5:6 => {6:[6,5,4,3,2,1]}
        lineages = set()
        parents = {}
        taxid_dataset = {}
        for t in datasets_mapping:
            taxid_dataset.update({t: t})

        taxid_lineage_name = "mapping_taxid-lineage.{}.txt".format(search_lineage)
        taxid_lineage_file_path = self.placement_description[taxid_lineage_name][3]
        with open(taxid_lineage_file_path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                lineage = line.strip().split("\t")[4]
                lineages.add(lineage)

                # for each line, e.g. 6\t1:2:3:4:5:6, create/update the lineage for each level
                # 6:[1,2,3,4,5,6], 5:[1,2,3,4,5], 4:[1,2,3,4], etc.
                levels = lineage.split(",")

                for i, t in enumerate(levels):
                    parents.update({t: levels[0: i + 1][::-1]})

        for t in parents:
            for p in parents[t]:  # get the deepest parent, not the root one
                if p in datasets_mapping:
                    taxid_dataset.update({t: p})
                    break
        # load json
        # load "tree" in a string
        # load placements
        # obtain a dict of taxid num of markers
        # figure out which taxid to use by using the highest number of markers and some extra rules

        try:
            with open(os.path.join(self.sepp_output_folder, "output_placement.json")) as json_file:
                data = json.load(json_file)
            tree = data["tree"]
            placements = data["placements"]
        except FileNotFoundError:
            raise Error("Placements failed. Try to rerun increasing the memory or select a lineage manually.")

        node_weight = {}
        n_p = 0
        for placement in placements:
            n_p += 1
            for individual_placement in placement["p"]:
                # find the taxid in tree
                node = individual_placement[0]

                match = re.findall(  # deal with weird character in the json file, see the output yourself.
                    # if this pattern is inconsistant with pplacer version, it may break buscoplacer.
                    "[^0-9][0-9]*:[0-9]*[^0-9]{0,1}[0-9]*[^0-9]{0,2}[0-9]*\[%s\]"
                    % node,
                    tree,
                )
                # extract taxid:
                try:
                    if re.match("^[A-Za-z]", match[0]):
                        taxid = match[0][7:].split(":")[0]
                    else:
                        taxid = match[0][1:].split(":")[0]
                except IndexError as e:
                    raise e
                if taxid_dataset[taxid] in node_weight:
                    node_weight[taxid_dataset[taxid]] += 1
                else:
                    node_weight[taxid_dataset[taxid]] = 1
                break  # Break here to keep only the best match. In my experience, keeping all does not change much.

        # from here, define which placement can be trusted
        max_markers = 0
        choice = []

        # taxid for which no threshold or minimal amount of placement should be considered.
        # If it is the best, go for it.
        no_rules = ["204428"]

        ratio = 2.5
        if search_lineage.split("_")[-2] == "archaea":
            ratio = 1.2
        min_markers = 12

        node_with_max_markers = None
        for n in node_weight:
            if node_weight[n] > max_markers:
                max_markers = node_weight[n]
                node_with_max_markers = n
        if node_with_max_markers in no_rules:
            choice = [node_with_max_markers]
        else:
            for n in node_weight:
                # if the ration between the best and the current one is not enough, keep both
                if node_weight[n] * ratio >= max_markers:
                    choice.append(n)
        if len(choice) > 1:
            # more than one taxid should be considered, pick the common ancestor
            choice = self._get_common_ancestor(choice, parents)
        elif len(choice) == 0:
            if search_lineage.split("_")[-2] == "bacteria":
                choice.append("2")
            elif search_lineage.split("_")[-2] == "archaea":
                choice.append("2157")
            elif search_lineage.split("_")[-2] == "eukaryota":
                choice.append("2759")
        if max_markers < min_markers and not (choice[0] in no_rules):
            if search_lineage.split("_")[-2] == "bacteria":
                key_taxid = "2"
            elif search_lineage.split("_")[-2] == "archaea":
                key_taxid = "2157"
            elif search_lineage.split("_")[-2] == "eukaryota":
                key_taxid = "2759"
            else:
                key_taxid = None  # unexpected. Should throw an exception or use assert.

            lineage = datasets_mapping[taxid_dataset[key_taxid]]

        else:

            lineage = datasets_mapping[taxid_dataset[choice[0]]]
        lineage = "{}_{}".format(lineage, "odb10")
        placed_markers = sum(node_weight.values())
        return [lineage, max_markers, placed_markers]

    @staticmethod
    def _get_common_ancestor(choice, parents):
        # starts with the parents of the first choice
        all_ancestors = set(parents[choice[0]])
        # order will be lost with sets, so keep in a list the lineage of one entry to later pick the deepest ancestor
        ordered_lineage = []
        for c in choice:
            if len(parents[c]) > len(ordered_lineage):
                # probably useless. Init with parents[choice[0] should work
                ordered_lineage = parents[c]
            # keep in set only entries that are in the currently explored lineage
            all_ancestors = all_ancestors.intersection(parents[c])

        # go through the ordered list of the deepest linage until you found a common ancestor.
        for parent in ordered_lineage:
            if parent in all_ancestors:
                return [parent]

    def Run(self, marker_gene_filepath):
        search_lineage = self.run_sepp(marker_gene_filepath)
        lineage, max_markers, placed_markers = self.pick_dataset(search_lineage)
        return lineage


### Analysis miniprot alignment ###
AminoAcid = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
             "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]


class GeneLabel(Enum):
    Single = 1
    Duplicated = 2
    Fragmented = 3
    Interspaced = 4
    Missing = 5


class MiniprotGffItems:
    def __init__(self):
        self.atn_seq = ""
        self.ata_seq = ""
        self.target_id = ""
        self.contig_id = ""
        self.protein_length = 0
        self.protein_start = 0
        self.protein_end = 0
        self.contig_start = 0
        self.contig_end = 0
        self.strand = ""
        self.score = 0
        self.rank = 0
        self.identity = 0
        self.positive = 0
        self.codons = []
        self.frameshift_events = 0
        self.frameshift_lengths = 0
        self.frame_shifts = []

    def show(self):
        return [self.atn_seq,
                self.ata_seq,
                self.target_id,
                self.contig_id,
                self.protein_length,
                self.protein_start,
                self.protein_end,
                self.contig_start,
                self.contig_end,
                self.strand,
                self.score,
                self.rank,
                self.identity,
                self.positive,
                "|".join(self.codons),
                self.frameshift_events,
                self.frameshift_lengths,
                self.frame_shifts]


def get_region_clusters(regions):
    sorted_regions = sorted(regions, key=lambda x: x[0], reverse=False)
    clusters = []
    for (start, stop) in sorted_regions:
        if not clusters:
            clusters.append([start, stop])
        else:
            last_cluster = clusters[-1]
            if last_cluster[0] <= start <= last_cluster[1]:
                # has overlap
                clusters[-1][0] = min(last_cluster[0], start)
                clusters[-1][1] = max(last_cluster[1], stop)
            else:
                clusters.append([start, stop])
    return clusters


class OutputFormat:
    def __init__(self):
        self.gene_label = None
        self.data_record = None


def find_frameshifts(cs_seq):
    frameshifts = []
    frameshift_events = 0
    frameshift_lengths = 0
    pt = r"[0-9]+[MIDFGNUV]"
    it = re.finditer(pt, cs_seq)
    for m in it:
        if m.group(0).endswith("F") or m.group(0).endswith("G"):
            frameshifts.append(m.group(0))
            frameshift_events += 1
            frameshift_lengths += int(m.group(0)[:-1])
    return frameshifts, frameshift_events, frameshift_lengths


def load_dbinfo(dbinfo_file):
    dbinfo = {}
    with open(dbinfo_file, "r") as f:
        for line in f:
            gene_id, db, link = line.strip().split("\t")
            dbinfo[gene_id] = [link, db]
    return dbinfo


class MiniprotAlignmentParser:
    def __init__(self, run_folder, gff_file, lineage, min_length_percent, min_diff, min_identity, min_complete,
                 min_rise, specified_contigs, autolineage, library_path=None):
        self.autolineage = autolineage
        self.run_folder = run_folder
        if not os.path.exists(run_folder):
            os.makedirs(run_folder)
        if lineage is None:
            self.completeness_output_file = os.path.join(self.run_folder, "summary.txt")
            self.full_table_output_file = os.path.join(self.run_folder, "full_table.tsv")
            self.full_table_busco_format_output_file = os.path.join(self.run_folder, "full_table_busco_format.tsv")
        else:
            if not lineage.endswith("_odb10"):
                lineage = lineage + "_odb10"
            self.run_folder = os.path.join(run_folder, lineage)
            if not os.path.exists(self.run_folder):
                os.makedirs(self.run_folder)
            self.completeness_output_file = os.path.join(run_folder, "summary.txt")
            self.full_table_output_file = os.path.join(self.run_folder, "full_table.tsv")
            self.full_table_busco_format_output_file = os.path.join(self.run_folder, "full_table_busco_format.tsv")
        self.library_path = library_path
        self.gff_file = gff_file
        self.lineage = lineage
        self.min_length_percent = min_length_percent
        self.min_diff = min_diff
        self.min_identity = min_identity
        self.min_complete = min_complete
        self.min_rise = min_rise
        self.specified_contigs = specified_contigs
        self.marker_gene_path = os.path.join(self.run_folder, "gene_marker.fasta")

    @staticmethod
    def parse_miniprot_records(gff_file, autolineage=False):
        items = MiniprotGffItems()
        with open(gff_file, "r") as gff:
            while True:
                line = gff.readline()
                if not line:
                    if items.target_id != "":
                        yield items
                    return
                if line.startswith("##gff-version"):
                    continue
                if line.startswith("##PAF"):
                    if items.target_id != "":
                        yield items
                    items.__init__()
                    fields = line.strip().split("\t")[1:]
                    items.target_id = fields[0]
                    items.protein_length = int(fields[1])
                    items.protein_start = int(fields[2])
                    items.protein_end = int(fields[3])
                    items.strand = fields[4]
                    items.contig_id = fields[5]
                    items.contig_start = int(fields[7])
                    items.contig_end = int(fields[8])
                    if fields[5] == "*":
                        ## Unmapped protein
                        items.score = 0
                        items.frameshift_events = 0
                        items.frameshift_lengths = 0
                        items.frame_shifts = []
                        continue
                    items.score = int(fields[13].strip().split(":")[2])
                    # Score = fields[12].strip().split(":")[2]
                    cg = fields[17].replace("cg:Z:", "")
                    cs = fields[18].replace("cs:Z:", "")
                    frame_shifts, frameshift_events, frameshift_lengths = find_frameshifts(cg)
                    items.frameshift_events = frameshift_events
                    items.frameshift_lengths = frameshift_lengths
                    items.frame_shifts = frame_shifts

                    if autolineage:
                        # This is for parse the ATN/ATA/AAS/AQA sequences in miniprot --aln
                        atn_line = gff.readline()
                        ata_line = gff.readline()
                        aas_line = gff.readline()
                        aqa_line = gff.readline()
                        items.atn_seq = atn_line.strip().split("\t")[1].replace("-", "")
                        ata_seq = ata_line.strip().split("\t")[1]
                        new_ata = []
                        for i in range(len(ata_seq)):
                            if ata_seq[i].upper() not in AminoAcid:
                                continue
                            else:
                                new_ata.append(ata_seq[i])
                        items.ata_seq = "".join(new_ata)
                else:
                    fields = line.strip().split("\t")
                    if fields[2] == "mRNA":
                        info_dict = dict(v.split("=") for v in fields[8].split()[0].split(";"))
                        items.rank = int(info_dict["Rank"])
                        items.identity = float(info_dict["Identity"])
                        items.positive = float(info_dict["Positive"])
                    if fields[2] == "CDS":
                        seq_id = fields[0]
                        codon_start = int(fields[3])  # 1-based
                        codon_end = int(fields[4])  # 1-based
                        codon_strand = fields[6]
                        info_dict = dict(x.split("=") for x in fields[8].split()[0].split(";"))
                        target_id = info_dict["Target"]
                        assert target_id == items.target_id and seq_id == items.contig_id
                        items.codons.append("{}_{}_{}".format(codon_start, codon_end, codon_strand))

    @staticmethod
    def record_1st_gene_label(dataframe, min_identity, min_complete):
        # check records with same tid of the best record
        output = OutputFormat()
        gene_id = dataframe.iloc[0]["Target_id"]
        # dataframe = dataframe[dataframe["Identity"] >= min_identity]
        if dataframe.shape[0] == 0:
            output.gene_label = GeneLabel.Missing
            return output
        elif dataframe.shape[0] == 1:
            if dataframe.iloc[0]["Protein_mapped_rate"] >= min_complete:
                output.gene_label = GeneLabel.Single
                output.data_record = dataframe.iloc[0]
                return output
            else:
                output.gene_label = GeneLabel.Fragmented
                output.data_record = dataframe.iloc[0]
                return output
        else:
            complete_regions = []
            fragmented_regions = []
            for i in range(dataframe.shape[0]):
                if dataframe.iloc[i]["Protein_mapped_rate"] >= min_complete:
                    complete_regions.append(
                        (dataframe.iloc[i]["Contig_id"], dataframe.iloc[i]["Start"], dataframe.iloc[i]["Stop"]))
                else:
                    fragmented_regions.append(
                        (dataframe.iloc[i]["Contig_id"], dataframe.iloc[i]["Start"], dataframe.iloc[i]["Stop"]))
            if len(complete_regions) == 0:
                output.gene_label = GeneLabel.Fragmented
                output.data_record = dataframe.iloc[0]
                return output
            elif len(complete_regions) == 1:
                output.gene_label = GeneLabel.Single
                output.data_record = dataframe.iloc[0]
                return output
            else:
                ctgs = [x[0] for x in complete_regions]
                if len(set(ctgs)) > 1:
                    output.gene_label = GeneLabel.Duplicated
                    output.data_record = dataframe
                    return output
                regions = [(x[1], x[2]) for x in complete_regions]
                clusters = get_region_clusters(regions)
                if len(clusters) == 1:
                    output.gene_label = GeneLabel.Single
                    output.data_record = dataframe.iloc[0]
                    return output
                else:
                    output.gene_label = GeneLabel.Duplicated
                    output.data_record = dataframe
                    return output

    @staticmethod
    def record_1st_2nd_gene_label(dataframe_1st, dataframe_2nd, min_identity, min_complete, min_rise):
        # check top 1st and 2nd records whether they are the same gene
        output = OutputFormat()
        # dataframe_1st = dataframe_1st[dataframe_1st["Identity"] >= min_identity]
        # dataframe_2nd = dataframe_2nd[dataframe_2nd["Identity"] >= min_identity]
        if dataframe_1st.shape[0] >= 1 and dataframe_2nd.shape[0] == 0:
            out = MiniprotAlignmentParser.record_1st_gene_label(dataframe_1st, min_identity, min_complete)
            return out
        if dataframe_1st.shape[0] == 0 and dataframe_2nd.shape[0] >= 1:
            out = MiniprotAlignmentParser.record_1st_gene_label(dataframe_2nd, min_identity, min_complete)
            return out
        if dataframe_1st.shape[0] == 0 and dataframe_2nd.shape[0] == 0:
            output.gene_label = GeneLabel.Missing
            return output
        else:
            gene_id1 = dataframe_1st.iloc[0]["Target_id"]
            gene_id2 = dataframe_2nd.iloc[0]["Target_id"]
            protein_length1 = dataframe_1st.iloc[0]["Protein_length"]
            protein_length2 = dataframe_2nd.iloc[0]["Protein_length"]

            label_length = defaultdict(list)
            out1 = MiniprotAlignmentParser.record_1st_gene_label(dataframe_1st, min_identity, min_complete)
            label_length[out1.gene_label].append(protein_length1)
            out2 = MiniprotAlignmentParser.record_1st_gene_label(dataframe_2nd, min_identity, min_complete)
            label_length[out2.gene_label].append(protein_length2)
            if label_length.keys() == {GeneLabel.Single}:
                output.gene_label = GeneLabel.Single
                output.data_record = dataframe_1st.iloc[0]
                return output
            elif label_length.keys() == {GeneLabel.Fragmented}:
                output.gene_label = GeneLabel.Fragmented
                output.data_record = dataframe_1st.iloc[0]
                return output
            elif label_length.keys() == {GeneLabel.Duplicated}:
                output.gene_label = GeneLabel.Duplicated
                output.data_record = dataframe_1st
                return output
            elif label_length.keys() == {GeneLabel.Single, GeneLabel.Fragmented}:
                if label_length[GeneLabel.Fragmented][0] > label_length[GeneLabel.Single][0] * (1 + min_rise):
                    output.gene_label = GeneLabel.Fragmented
                    if out1.gene_label == GeneLabel.Fragmented:
                        output.data_record = dataframe_1st.iloc[0]
                    elif out2.gene_label == GeneLabel.Fragmented:
                        output.data_record = dataframe_2nd.iloc[0]
                    else:
                        raise ValueError
                    return output
                else:
                    output.gene_label = GeneLabel.Single
                    if out1.gene_label == GeneLabel.Single:
                        output.data_record = dataframe_1st.iloc[0]
                    elif out2.gene_label == GeneLabel.Single:
                        output.data_record = dataframe_2nd.iloc[0]
                    else:
                        raise ValueError
                    return output
            elif label_length.keys() == {GeneLabel.Single, GeneLabel.Duplicated}:
                if label_length[GeneLabel.Duplicated][0] > label_length[GeneLabel.Single][0] * (1 + min_rise):
                    output.gene_label = GeneLabel.Duplicated
                    if out1.gene_label == GeneLabel.Duplicated:
                        output.data_record = dataframe_1st
                    elif out2.gene_label == GeneLabel.Duplicated:
                        output.data_record = dataframe_2nd
                    else:
                        raise ValueError
                    return output
                else:
                    output.gene_label = GeneLabel.Single
                    if out1.gene_label == GeneLabel.Single:
                        output.data_record = dataframe_1st.iloc[0]
                    elif out2.gene_label == GeneLabel.Single:
                        output.data_record = dataframe_2nd.iloc[0]
                    else:
                        raise ValueError
                    return output
            elif label_length.keys() == {GeneLabel.Fragmented, GeneLabel.Duplicated}:
                if label_length[GeneLabel.Fragmented][0] > label_length[GeneLabel.Duplicated][0] * (1 + min_rise):
                    output.gene_label = GeneLabel.Fragmented
                    if out1.gene_label == GeneLabel.Fragmented:
                        output.data_record = dataframe_1st.iloc[0]
                    elif out2.gene_label == GeneLabel.Fragmented:
                        output.data_record = dataframe_2nd.iloc[0]
                    else:
                        raise ValueError
                    return output
                else:
                    output.gene_label = GeneLabel.Duplicated
                    if out1.gene_label == GeneLabel.Duplicated:
                        output.data_record = dataframe_1st
                    elif out2.gene_label == GeneLabel.Duplicated:
                        output.data_record = dataframe_2nd
                    else:
                        raise ValueError
                    return output
            else:
                print("Error: wrong permutation of label_length!")
                raise ValueError

    @staticmethod
    def Ost_eval(dataframe, difficial_rate, min_identity, min_complete, min_rise):
        if dataframe.shape[0] == 0:
            output = OutputFormat()
            output.gene_label = GeneLabel.Missing
            return output
        if dataframe.shape[0] == 1:
            return MiniprotAlignmentParser.record_1st_gene_label(
                dataframe[dataframe["Target_id"] == dataframe.iloc[0]["Target_id"]], min_identity, min_complete)
        record_1st = dataframe.iloc[0]
        record_1st_tid = record_1st["Target_id"]
        record_2nd = dataframe.iloc[1]
        record_2nd_tid = record_2nd["Target_id"]
        if (record_1st["I+L"] - record_2nd["I+L"]) / (record_2nd["I+L"] + 1e-9) >= difficial_rate:
            return MiniprotAlignmentParser.record_1st_gene_label(dataframe[dataframe["Target_id"] == record_1st_tid],
                                                                 min_identity, min_complete)
        else:
            if record_1st_tid == record_2nd_tid:
                return MiniprotAlignmentParser.record_1st_gene_label(
                    dataframe[dataframe["Target_id"] == record_1st_tid], min_identity, min_complete)
            else:
                return MiniprotAlignmentParser.record_1st_2nd_gene_label(
                    dataframe[dataframe["Target_id"] == record_1st_tid],
                    dataframe[dataframe["Target_id"] == record_2nd_tid], min_identity, min_complete, min_rise)

    @staticmethod
    def refine_fragmented(dataframe):
        target_ids = dataframe["Target_id"].unique()
        identity_plus_length = defaultdict(int)
        for tid in target_ids:
            sub_dataframe = dataframe[dataframe["Target_id"] == tid]
            if sub_dataframe.shape[0] == 1:
                identity_plus_length[tid] = sub_dataframe["I+L"].iloc[0]
            else:
                regions = sorted(sub_dataframe[["Protein_Start", "Protein_End", "Contig_id", "Start", "Stop",
                                                "I+L"]].values.tolist(), key=lambda x: x[0])
                clusters = []
                for region in regions:
                    if len(clusters) == 0:
                        clusters.append(region)
                    else:
                        if region[0] > clusters[-1][1]:
                            clusters.append(region)
                        else:
                            if region[1] - region[0] > clusters[-1][1] - clusters[-1][0]:
                                clusters[-1] = region
                            else:
                                pass
                if len(clusters) == 1:
                    identity_plus_length[tid] = clusters[0][5]
                else:
                    new_clusters = defaultdict(list)
                    for region in clusters:
                        if len(new_clusters.keys()) == 0:
                            new_clusters[region[2]].append(region)
                        else:
                            if region[2] not in new_clusters:
                                new_clusters[region[2]].append(region)
                            else:
                                cid = region[2]
                                if region[3] > new_clusters[cid][-1][4]:
                                    new_clusters[cid].append(region)
                                else:
                                    if region[4] - region[3] > new_clusters[cid][-1][4] - new_clusters[cid][-1][3]:
                                        new_clusters[cid][-1] = region
                                    else:
                                        pass
                    for cid in new_clusters.keys():
                        for region in new_clusters[cid]:
                            identity_plus_length[tid] += region[5]
        identity_plus_length = sorted(identity_plus_length.items(), key=lambda x: x[1], reverse=True)
        tid = identity_plus_length[0][0]
        output = OutputFormat()
        sub_dataframe = dataframe[dataframe["Target_id"] == tid]
        if sub_dataframe.shape[0] > 1 and identity_plus_length[0][1] > dataframe.iloc[0]["I+L"]:
            output.gene_label = GeneLabel.Interspaced
            output.data_record = sub_dataframe.iloc[0]
        else:
            output.gene_label = GeneLabel.Fragmented
            output.data_record = sub_dataframe.iloc[0]
        return output

    def Run(self):
        single_genes = []
        duplicate_genes = []
        fragmented_genes = []
        interspaced_genes = []
        missing_genes = []
        records = []
        single_complete_proteins = []
        gff_file = self.gff_file
        try:
            reader = iter(self.parse_miniprot_records(gff_file, self.autolineage))
            for items in reader:
                (Atn_seq, Ata_seq, Target_id, Contig_id, Protein_length, Protein_Start, Protein_End, Start, Stop,
                 Strand, Score, Rank, Identity, Positive, Codons, Frameshift_events, Frameshift_lengths,
                 Frame_shifts) = items.show()
                Target_species = Target_id.split("_")[0]
                records.append([Target_species, Target_id, Contig_id, Protein_length, Protein_Start, Protein_End,
                                Protein_End - Protein_Start, (Protein_End - Protein_Start) / Protein_length, Start,
                                Stop, Stop - Start, Strand, Rank, Identity, Positive,
                                (Protein_End - Protein_Start) / Protein_length * Identity,
                                Frameshift_events, Frameshift_lengths, Score, Atn_seq, Ata_seq, Codons])
        except StopIteration:
            pass
        records_df = pd.DataFrame(records, columns=["Target_species", "Target_id", "Contig_id", "Protein_length",
                                                    "Protein_Start", "Protein_End", "Protein_mapped_length",
                                                    "Protein_mapped_rate", "Start", "Stop", "Genome_mapped_length",
                                                    "Strand", "Rank", "Identity", "Positive", "I+L",
                                                    "Frameshift_events", "Frameshift_lengths", "Score", "Atn_seq",
                                                    "Ata_seq", "Codons"])
        all_species = records_df["Target_species"].unique()
        all_contigs = records_df["Contig_id"].unique()
        if self.specified_contigs is not None:
            if len(set(all_contigs) & set(self.specified_contigs)) == 0:
                raise Exception("No contigs found in the specified contigs!")
        grouped_df = records_df.groupby(["Target_species"])
        full_table_writer = open(self.full_table_output_file, "w")
        full_table_writer.write(
            "Gene\tStatus\tSequence\tGene Start\tGene End\tStrand\tScore\tLength\tIdentity\tFraction\tFrameshift events\tBest gene\tCodons\n")
        dbinfo = None
        if self.library_path is not None and self.lineage is not None:
            dbinfo_path = os.path.join(self.library_path, self.lineage, "links_to_ODB10.txt")
            if os.path.exists(dbinfo_path):
                dbinfo = load_dbinfo(dbinfo_path)
        full_table_busco_format_writer = open(self.full_table_busco_format_output_file, "w")
        if dbinfo is None:
            full_table_busco_format_writer.write(
                "# Busco id\tStatus\tSequence\tGene Start\tGene End\tStrand\tScore\tLength\n")
        else:
            full_table_busco_format_writer.write(
                "# Busco id\tStatus\tSequence\tGene Start\tGene End\tStrand\tScore\tLength\tOrthoDB url\tDescription\n")
        for gene_id in all_species:
            mapped_records = grouped_df.get_group(gene_id)
            mapped_records = mapped_records.sort_values(by=["I+L"], ascending=False)
            if self.specified_contigs is not None:
                mapped_records = mapped_records[mapped_records["Contig_id"].isin(self.specified_contigs)]
            pass_tids = mapped_records[(mapped_records["Protein_mapped_rate"] >= self.min_length_percent) | (
                    mapped_records["Identity"] >= self.min_identity)]["Target_id"].unique()

            # mean std filter
            if len(pass_tids) >= 5:
                lengths_dict = {}
                for tid in pass_tids:
                    lengths_dict[tid] = mapped_records[mapped_records["Target_id"] == tid].iloc[0]["Protein_length"]
                length_df = pd.DataFrame.from_dict(lengths_dict, orient="index", columns=["Protein_length"])
                mean_v = length_df["Protein_length"].mean()
                std_v = length_df["Protein_length"].std()
                lower_bound = mean_v - 1.5 * std_v
                upper_bound = mean_v + 1.5 * std_v
                pass_tids = length_df[(length_df["Protein_length"] >= lower_bound) &
                                      (length_df["Protein_length"] <= upper_bound)].index.tolist()

            if len(pass_tids) > 0:
                mapped_records = mapped_records[mapped_records["Target_id"].isin(pass_tids)]
                output = self.Ost_eval(mapped_records, self.min_diff, self.min_identity, self.min_complete,
                                       self.min_rise)
                if output.gene_label == GeneLabel.Single:
                    single_complete_proteins.append(
                        ">{}\n{}\n".format(output.data_record["Target_id"], output.data_record["Ata_seq"]))

                if output.gene_label == GeneLabel.Fragmented:
                    output = self.refine_fragmented(mapped_records)
            else:
                output = OutputFormat()
                output.gene_label = GeneLabel.Missing

            if output.gene_label == GeneLabel.Missing:
                full_table_writer.write("{}\t{}\n".format(gene_id, output.gene_label.name))
                full_table_busco_format_writer.write("{}\t{}\n".format(gene_id, output.gene_label.name))
            else:
                assert output.data_record.shape[0] >= 1
                if output.data_record.ndim == 1:
                    full_table_writer.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".
                                            format(output.data_record["Target_species"],
                                                   output.gene_label.name,
                                                   output.data_record["Contig_id"],
                                                   output.data_record["Start"],
                                                   output.data_record["Stop"],
                                                   output.data_record["Strand"],
                                                   output.data_record["Score"],
                                                   output.data_record["Protein_mapped_length"],
                                                   output.data_record["Identity"],
                                                   output.data_record["Protein_mapped_rate"],
                                                   output.data_record["Frameshift_events"],
                                                   output.data_record["Target_id"],
                                                   output.data_record["Codons"]))
                    if output.gene_label.name == "Single":
                        status = "Complete"
                    elif output.gene_label.name == "Interspaced":
                        status = "Fragmented"
                    else:
                        status = output.gene_label.name
                    if dbinfo is None:
                        full_table_busco_format_writer.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".
                                                             format(output.data_record["Target_species"],
                                                                    status,
                                                                    output.data_record["Contig_id"],
                                                                    output.data_record["Start"],
                                                                    output.data_record["Stop"],
                                                                    output.data_record["Strand"],
                                                                    output.data_record["Score"],
                                                                    output.data_record["Protein_mapped_length"]))
                    else:
                        try:
                            full_table_busco_format_writer.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".
                                                                 format(output.data_record["Target_species"],
                                                                        status,
                                                                        output.data_record["Contig_id"],
                                                                        output.data_record["Start"],
                                                                        output.data_record["Stop"],
                                                                        output.data_record["Strand"],
                                                                        output.data_record["Score"],
                                                                        output.data_record["Protein_mapped_length"],
                                                                        dbinfo[output.data_record["Target_species"]][0],
                                                                        dbinfo[output.data_record["Target_species"]][
                                                                            1]))
                        except KeyError:
                            full_table_busco_format_writer.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".
                                                                 format(output.data_record["Target_species"],
                                                                        status,
                                                                        output.data_record["Contig_id"],
                                                                        output.data_record["Start"],
                                                                        output.data_record["Stop"],
                                                                        output.data_record["Strand"],
                                                                        output.data_record["Score"],
                                                                        output.data_record["Protein_mapped_length"],
                                                                        "*",
                                                                        "*"))
                else:
                    for dri in range(output.data_record.shape[0]):
                        full_table_writer.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".
                                                format(output.data_record.iloc[dri]["Target_species"],
                                                       output.gene_label.name,
                                                       output.data_record.iloc[dri]["Contig_id"],
                                                       output.data_record.iloc[dri]["Start"],
                                                       output.data_record.iloc[dri]["Stop"],
                                                       output.data_record.iloc[dri]["Strand"],
                                                       output.data_record.iloc[dri]["Score"],
                                                       output.data_record.iloc[dri]["Protein_mapped_length"],
                                                       output.data_record.iloc[dri]["Identity"],
                                                       output.data_record.iloc[dri]["Protein_mapped_rate"],
                                                       output.data_record.iloc[dri]["Frameshift_events"],
                                                       output.data_record.iloc[dri]["Target_id"],
                                                       output.data_record.iloc[dri]["Codons"]))
                        if output.gene_label.name == "Single":
                            status = "Complete"
                        elif output.gene_label.name == "Interspaced":
                            status = "Fragmented"
                        else:
                            status = output.gene_label.name
                        if dbinfo is None:
                            full_table_busco_format_writer.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".
                                                                 format(output.data_record.iloc[dri]["Target_species"],
                                                                        status,
                                                                        output.data_record.iloc[dri]["Contig_id"],
                                                                        output.data_record.iloc[dri]["Start"],
                                                                        output.data_record.iloc[dri]["Stop"],
                                                                        output.data_record.iloc[dri]["Strand"],
                                                                        output.data_record.iloc[dri]["Score"],
                                                                        output.data_record.iloc[dri][
                                                                            "Protein_mapped_length"]))
                        else:
                            try:
                                full_table_busco_format_writer.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".
                                                                     format(
                                    output.data_record.iloc[dri]["Target_species"],
                                    status,
                                    output.data_record.iloc[dri]["Contig_id"],
                                    output.data_record.iloc[dri]["Start"],
                                    output.data_record.iloc[dri]["Stop"],
                                    output.data_record.iloc[dri]["Strand"],
                                    output.data_record.iloc[dri]["Score"],
                                    output.data_record.iloc[dri]["Protein_mapped_length"],
                                    dbinfo[output.data_record.iloc[dri]["Target_species"]][0],
                                    dbinfo[output.data_record.iloc[dri]["Target_species"]][1]))
                            except KeyError:
                                full_table_busco_format_writer.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".
                                                                     format(
                                    output.data_record.iloc[dri]["Target_species"],
                                    status,
                                    output.data_record.iloc[dri]["Contig_id"],
                                    output.data_record.iloc[dri]["Start"],
                                    output.data_record.iloc[dri]["Stop"],
                                    output.data_record.iloc[dri]["Strand"],
                                    output.data_record.iloc[dri]["Score"],
                                    output.data_record.iloc[dri]["Protein_mapped_length"],
                                    "*",
                                    "*"))
            if output.gene_label == GeneLabel.Single:
                single_genes.append(gene_id)
            elif output.gene_label == GeneLabel.Duplicated:
                duplicate_genes.append(gene_id)
            elif output.gene_label == GeneLabel.Fragmented:
                fragmented_genes.append(gene_id)
            elif output.gene_label == GeneLabel.Interspaced:
                interspaced_genes.append(gene_id)
            elif output.gene_label == GeneLabel.Missing:
                missing_genes.append(gene_id)
            else:
                print("Error: output.gene_label!")
                raise ValueError
        full_table_writer.close()
        full_table_busco_format_writer.close()

        total_genes = len(all_species)
        d = total_genes - len(single_genes) - len(duplicate_genes) - len(fragmented_genes) - len(
            interspaced_genes) - len(missing_genes)
        print()
        print("S:{:.2f}%, {}".format(len(single_genes) / total_genes * 100, len(single_genes)))
        print("D:{:.2f}%, {}".format(len(duplicate_genes) / total_genes * 100, len(duplicate_genes)))
        print("F:{:.2f}%, {}".format(len(fragmented_genes) / total_genes * 100, len(fragmented_genes)))
        print("I:{:.2f}%, {}".format(len(interspaced_genes) / total_genes * 100, len(interspaced_genes)))
        print("M:{:.2f}%, {}".format((len(missing_genes) + d) / total_genes * 100, len(missing_genes) + d))
        print("N:{}".format(total_genes))
        print()
        with open(self.completeness_output_file, 'a') as fout:
            if self.lineage is not None:
                fout.write("## lineage: {}\n".format(self.lineage))
            else:
                fout.write("## lineage: xx_xx\n")
            fout.write("S:{:.2f}%, {}\n".format(len(single_genes) / total_genes * 100, len(single_genes)))
            fout.write("D:{:.2f}%, {}\n".format(len(duplicate_genes) / total_genes * 100, len(duplicate_genes)))
            fout.write("F:{:.2f}%, {}\n".format(len(fragmented_genes) / total_genes * 100, len(fragmented_genes)))
            fout.write("I:{:.2f}%, {}\n".format(len(interspaced_genes) / total_genes * 100, len(interspaced_genes)))
            fout.write(
                "M:{:.2f}%, {}\n".format((len(missing_genes) + d) / total_genes * 100, len(missing_genes) + d))
            fout.write("N:{}\n".format(total_genes))

        with open(self.marker_gene_path, "w") as fout:
            for x in single_complete_proteins:
                fout.write(x)


### Minibusco Runner ###
class MinibuscoRunner:
    def __init__(self, assembly_path, output_folder, library_path, lineage, autolineage, nthreads,
                 miniprot_execute_command, sepp_execute_command, min_diff, min_length_percent, min_identity,
                 min_complete, min_rise, specified_contigs):
        if lineage is not None:
            if not lineage.endswith("_odb10"):
                lineage = lineage + "_odb10"
        self.lineage = lineage
        self.autolineage = autolineage
        self.output_folder = output_folder
        self.library_path = library_path
        self.assembly_path = assembly_path
        self.min_diff = min_diff
        self.min_length_percent = min_length_percent
        self.min_identity = min_identity
        self.min_complete = min_complete
        self.min_rise = min_rise
        self.specified_contigs = specified_contigs

        self.miniprot_runner = MiniprotRunner(miniprot_execute_command, nthreads, autolineage)
        self.downloader = Downloader(library_path)

        sepp_output_path = os.path.join(output_folder, "sepp_output")
        sepp_tmp_path = os.path.join(output_folder, "sepp_tmp")
        self.lineage_searcher = AutoLineager(sepp_output_path, sepp_tmp_path, library_path, nthreads,
                                             sepp_execute_command)

    def Run(self):
        begin_time = time.time()
        if self.autolineage:
            lineage = "eukaryota_odb10"
        else:
            lineage = self.lineage
        download_lineage_start_time = time.time()
        self.downloader.download_lineage(lineage)
        download_lineage_end_time = time.time()
        print("lineage: {}".format(lineage))
        lineage_filepath = os.path.join(self.downloader.lineage_description[lineage][3], "refseq_db.faa.gz")
        alignment_output_dir = os.path.join(self.output_folder, lineage)
        if not os.path.exists(alignment_output_dir):
            os.makedirs(alignment_output_dir)

        run_miniprot_start_time = time.time()
        miniprot_output_path = self.miniprot_runner.run_miniprot(self.assembly_path,
                                                                 lineage_filepath,
                                                                 alignment_output_dir)
        run_miniprot_end_time = time.time()
        analysis_miniprot_start_time = time.time()
        miniprot_alignment_parser = MiniprotAlignmentParser(run_folder=self.output_folder,
                                                            gff_file=miniprot_output_path,
                                                            lineage=lineage,
                                                            min_diff=self.min_diff,
                                                            min_length_percent=self.min_length_percent,
                                                            min_identity=self.min_identity,
                                                            min_complete=self.min_complete,
                                                            min_rise=self.min_rise,
                                                            specified_contigs=self.specified_contigs,
                                                            autolineage=self.autolineage,
                                                            library_path=self.library_path)

        if os.path.exists(miniprot_alignment_parser.completeness_output_file):
            os.remove(miniprot_alignment_parser.completeness_output_file)
        miniprot_alignment_parser.Run()
        analysis_miniprot_end_time = time.time()
        if self.autolineage:
            autolineage_start_time = time.time()
            marker_genes_filepath = miniprot_alignment_parser.marker_gene_path
            best_match_lineage = self.lineage_searcher.Run(marker_genes_filepath)
            print("best_match_lineage: {}".format(best_match_lineage))
            autolineage_end_time = time.time()
            if best_match_lineage == lineage:
                end_time = time.time()
                print("## Download lineage: {:.2f}(s)".format(download_lineage_end_time - download_lineage_start_time))
                print("## Run miniprot: {:.2f}(s)".format(run_miniprot_end_time - run_miniprot_start_time))
                print(
                    "## Analyze miniprot: {:.2f}(s)".format(analysis_miniprot_end_time - analysis_miniprot_start_time))
                print("## Autolineage: {:.2f}(s)".format(autolineage_end_time - autolineage_start_time))
                print("## Total runtime: {:.2f}(s)".format(end_time - begin_time))
                return
            self.downloader.download_lineage(best_match_lineage)
            lineage = best_match_lineage
            lineage_filepath = os.path.join(self.downloader.lineage_description[lineage][3], "refseq_db.faa.gz")
            alignment_output_dir = os.path.join(self.output_folder, lineage)
            if not os.path.exists(alignment_output_dir):
                os.makedirs(alignment_output_dir)
            second_run_miniprot_start_time = time.time()
            miniprot_output_path = self.miniprot_runner.run_miniprot(self.assembly_path,
                                                                     lineage_filepath,
                                                                     alignment_output_dir)
            second_run_miniprot_end_time = time.time()
            second_analysis_miniprot_start_time = time.time()
            miniprot_alignment_parser = MiniprotAlignmentParser(run_folder=self.output_folder,
                                                                gff_file=miniprot_output_path,
                                                                lineage=lineage,
                                                                min_diff=self.min_diff,
                                                                min_length_percent=self.min_length_percent,
                                                                min_identity=self.min_identity,
                                                                min_complete=self.min_complete,
                                                                min_rise=self.min_rise,
                                                                specified_contigs=self.specified_contigs,
                                                                autolineage=self.autolineage,
                                                                library_path=self.library_path)
            miniprot_alignment_parser.Run()
            second_analysis_miniprot_end_time = time.time()
        end_time = time.time()
        print("## Download lineage: {:.2f}(s)".format(download_lineage_end_time - download_lineage_start_time))
        print("## Run miniprot: {:.2f}(s)".format(run_miniprot_end_time - run_miniprot_start_time))
        print("## Analyze miniprot: {:.2f}(s)".format(analysis_miniprot_end_time - analysis_miniprot_start_time))
        if self.autolineage:
            print("## Autolineage: {:.2f}(s)".format(autolineage_end_time - autolineage_start_time))
            print("## Second run miniprot: {:.2f}(s)".format(
                second_run_miniprot_end_time - second_run_miniprot_start_time))
            print("## Second analyze miniprot: {:.2f}(s)".format(
                second_analysis_miniprot_end_time - second_analysis_miniprot_start_time))
        print("## Total runtime: {:.2f}(s)".format(end_time - begin_time))


### main function ###

def download(args):
    downloader = Downloader(args.library_path)
    for lineage in args.lineage:
        downloader.download_lineage(lineage)


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
            for file in os.listdir(args.library_path):
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


def miniprot(args):
    mr = MiniprotRunner(args.exec_path, args.threads, autolineage_mode=False)
    mr.run_miniprot(args.assembly, args.protein, args.outdir)


def analysis(args):
    ar = MiniprotAlignmentParser(run_folder=args.output_dir,
                                 gff_file=args.gff,
                                 lineage=None,
                                 min_length_percent=args.min_length_percent,
                                 min_diff=args.min_diff,
                                 min_identity=args.min_identity,
                                 min_complete=args.min_complete,
                                 min_rise=args.min_rise,
                                 specified_contigs=args.specified_contigs,
                                 autolineage=False)
    ar.Run()


def run(args):
    assembly_path = args.assembly_path
    output_folder = args.output_dir
    library_path = args.library_path
    lineage = args.lineage
    autolineage = args.autolineage
    nthreads = args.threads
    miniprot_execute_command = args.miniprot_execute_path
    sepp_execute_command = args.sepp_execute_path
    min_diff = args.min_diff
    min_length_percent = args.min_length_percent
    min_identity = args.min_identity
    min_complete = args.min_complete
    min_rise = args.min_rise
    specified_contigs = args.specified_contigs

    if lineage is None and autolineage is False:
        sys.exit(
            "\n Usage error: Please specify the lineage name with -l. e.g. eukaryota, primates, saccharomycetes etc."
            "\n Or specify --autolineage and --sepp_execute_path to automaticly search the best matching lineage\n")
    if autolineage and sepp_execute_command is None:
        sys.exit("\n Usage error: Please specify the path to SEPP executable file with --sepp_execute_path\n")

    mr = MinibuscoRunner(assembly_path=assembly_path,
                         output_folder=output_folder,
                         library_path=library_path,
                         lineage=lineage,
                         autolineage=autolineage,
                         nthreads=nthreads,
                         miniprot_execute_command=miniprot_execute_command,
                         sepp_execute_command=sepp_execute_command,
                         min_diff=min_diff,
                         min_length_percent=min_length_percent,
                         min_identity=min_identity,
                         min_complete=min_complete,
                         min_rise=min_rise,
                         specified_contigs=specified_contigs)
    mr.Run()


### main.py
def main():
    parser = argparse.ArgumentParser(description="MiniBusco")
    subparser = parser.add_subparsers(dest="command", help="Minibusco modules help", required=True)

    ### sub-command: download
    download_parser = subparser.add_parser("download", help="Download specified BUSCO lineages")
    download_parser.add_argument("-l", "--lineage", type=str, nargs='+',
                                 help="Specify the names of the BUSCO lineages to be downloaded. (e.g. eukaryota, primates, saccharomycetes etc.)",
                                 required=True)
    download_parser.add_argument("--library_path", type=str, default="mb_downloads",
                                 help="The destination folder to store the downloaded lineage files."
                                      "If not specified, a folder named \"mb_downloads\" will be created on the current running path.")
    download_parser.set_defaults(func=download)

    ### sub-command: list
    list_parser = subparser.add_parser("list", help="List local or remote BUSCO lineages")
    list_parser.add_argument("--remote", action="store_true", help="List remote BUSCO lineages")
    list_parser.add_argument("--local", action="store_true", help="List local BUSCO lineages")
    list_parser.add_argument("--library_path", type=str, help="Folder path to stored lineages. ", default=None)
    list_parser.set_defaults(func=list_lineages)

    ### sub-command: run_miniprot
    run_miniprot_parser = subparser.add_parser("run_miniprot", help="Run miniprot alignment")
    run_miniprot_parser.add_argument("-a", "--assembly", type=str, help="Input genome file in FASTA format",
                                     required=True)
    run_miniprot_parser.add_argument("-p", "--protein", type=str, help="Input protein file", required=True)
    run_miniprot_parser.add_argument("-o", "--outdir", type=str, help="Miniprot alignment output directory",
                                     required=True)
    run_miniprot_parser.add_argument("-t", "--threads", type=int, help="Number of threads to use", default=1)
    run_miniprot_parser.add_argument("--exec_path", type=str, help="Path to miniprot executable", default=None)
    run_miniprot_parser.set_defaults(func=miniprot)

    ### sub-command: analysis
    analysis_parser = subparser.add_parser("analysis",
                                           help="Evaluate genome completeness from provided miniprot alignment")
    analysis_parser.add_argument("-g", "--gff", type=str, help="Miniprot output gff file", required=True)
    analysis_parser.add_argument("-o", "--output_dir", type=str, help="Output analysis folder", required=True)
    analysis_parser.add_argument("--specified_contigs", type=str, nargs='+', default=None,
                                 help="Specify the contigs to be evaluated, e.g. chr1 chr2 chr3. If not specified, all contigs will be evaluated.")
    analysis_parser.add_argument("--min_diff", type=float, default=0.2,
                                 help="The thresholds for the best matching and second best matching.")
    analysis_parser.add_argument("--min_identity", type=float, default=0.4,
                                 help="The identity threshold for valid mapping results. [0, 1]")
    analysis_parser.add_argument("--min_length_percent", type=float, default=0.6,
                                 help="The fraction of protein for valid mapping results.")
    analysis_parser.add_argument("--min_complete", type=float, default=0.9,
                                 help="The length threshold for complete gene.")
    analysis_parser.add_argument("--min_rise", type=float, default=0.5,
                                 help="Minimum length threshold to make dupicate take precedence over single or fragmented over single/duplicate.")
    analysis_parser.set_defaults(func=analysis)

    ### sub-command: run
    run_parser = subparser.add_parser("run",
                                      help="Run minibusco including miniprot alignment and completeness evaluation")
    run_parser.add_argument("-a", "--assembly_path", type=str, help="Input genome file in FASTA format.", required=True)
    run_parser.add_argument("-o", "--output_dir", type=str, help="The output folder.", required=True)
    run_parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to use")
    run_parser.add_argument("-l", "--lineage", type=str, default=None,
                            help="Specify the name of the BUSCO lineage to be used. (e.g. eukaryota, primates, saccharomycetes etc.)")
    run_parser.add_argument("--library_path", type=str, default="mb_downloads",
                            help="Folder path to download lineages or already downloaded lineages. "
                                 "If not specified, a folder named \"mb_downloads\" will be created on the current running path by default to store the downloaded lineage files.")
    run_parser.add_argument("--specified_contigs", type=str, nargs='+', default=None,
                            help="Specify the contigs to be evaluated, e.g. chr1 chr2 chr3. If not specified, all contigs will be evaluated.")
    run_parser.add_argument("--miniprot_execute_path", type=str, default=None,
                            help="Path to miniprot executable")
    run_parser.add_argument("--autolineage", action="store_true",
                            help="Automatically search for the best matching lineage without specifying lineage file.")
    run_parser.add_argument("--sepp_execute_path", type=str, default=None,
                            help="Path to sepp executable")
    run_parser.add_argument("--min_diff", type=float, default=0.2,
                            help="The thresholds for the best matching and second best matching.")
    run_parser.add_argument("--min_identity", type=float, default=0.4,
                            help="The identity threshold for valid mapping results.")
    run_parser.add_argument("--min_length_percent", type=float, default=0.6,
                            help="The fraction of protein for valid mapping results.")
    run_parser.add_argument("--min_complete", type=float, default=0.9,
                            help="The length threshold for complete gene.")
    run_parser.add_argument("--min_rise", type=float, default=0.5,
                            help="Minimum length threshold to make dupicate take precedence over single or fragmented over single/duplicate.")
    run_parser.set_defaults(func=run)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
