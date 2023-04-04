import os
import json
import re
import shutil

from .DownloadLineage import Downloader
from .utils import Error

# logger = MinibuscoLogger(__name__).getlog()


class AutoLineager:
    def __init__(self, sepp_output_directory, sepp_tmp_directory, config):
        library_path = config.library_path

        self.sepp_output_folder = sepp_output_directory
        self.sepp_tmp_folder = sepp_tmp_directory

        self.threads = config.threads
        # self.run_folder = config.autolineage_folder
        self.downloader = Downloader(library_path)
        self.lineage_description = self.downloader.lineage_description
        self.placement_description = self.downloader.placement_description
        self.library_folder = self.downloader.download_dir
        self.placement_file_folder = self.downloader.placement_dir

        # TODO: get sepp path
        self.sepp_execute_command = "run_sepp.py"

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


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Run AutoLineage')
    parser.add_argument("--autolineage_folder", dest="autolineage_folder", help="Path to AutoLineage folder", type=str,
                        required=True)
    parser.add_argument("--sepp_output_directory", dest="sepp_output_directory", help="Path to SEPP output directory",
                        type=str, required=True)
    parser.add_argument("--threads", dest="threads", help="Number of threads to use", type=int, default=1)
    parser.add_argument("--library_path", dest="library_path", help="Path to the library folder", default=None,
                        required=False)
    args = parser.parse_args()

    auto_lineage = AutoLineager(args)
    lineage = auto_lineage.Run()
    print(lineage)
