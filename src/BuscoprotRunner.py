import os.path
import argparse
import shutil

# 1. run miniprot on lineage ["archaea_odb10", "bacteria_odb10", "eukaryota_odb10"]
# 2. run AnalisisMiniprot on the output of miniprot
# 3. get the highest complete gene percentage
# 4. run repp on the lineage with the highest complete gene percentage
# 5. get the most likely sub-lineage
# 6. run miniprot on the most likely sub-lineage
# 7. run AnalisisMiniprot on the output of miniprot
# 8. get final output

from .RunMiniprot import MiniprotRunner
from .AnalysisMiniprot import MiniprotAlignmentParser
from .DownloadLineage import Downloader
from .AutoLineage import AutoLineager
from .utils import MinibuscoLogger

logger = MinibuscoLogger().getlog(__name__)

class BuscoprotRunner:
    def __init__(self, config):
        autolineage = config.autolineage
        library_path = config.library_path
        output_folder = config.output_dir
        assembly_path = config.assembly_path
        if config.lineage is not None:
            self.lineage = config.lineage + "_odb10"
        else:
            self.lineage = config.lineage
        # TODO: get miniprot path
        miniprot_execute_command = "miniprot"

        self.autolineage = autolineage
        self.output_folder = output_folder
        self.assembly_path = assembly_path
        self.config = config

        self.miniprot_runner = MiniprotRunner(miniprot_execute_command, config)
        self.downloader = Downloader(library_path)

        sepp_output_path = os.path.join(output_folder, "sepp_output")
        sepp_tmp_path = os.path.join(output_folder, "sepp_tmp")
        self.lineage_searcher = AutoLineager(sepp_output_path, sepp_tmp_path, config)

    def Run(self):
        if self.autolineage:
            lineage = "eukaryota_odb10"
        else:
            lineage = self.lineage

        self.downloader.download_lineage(lineage)
        logger.info("lineage: {}".format(lineage))
        lineage_filepath = os.path.join(self.downloader.lineage_description[lineage][3], "refseq_db.faa.gz")
        output_dir = os.path.join(self.output_folder, lineage)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        miniprot_output_path = self.miniprot_runner.run_miniprot(self.assembly_path, lineage_filepath, output_dir)
        miniprot_alignment_parser = MiniprotAlignmentParser(output_dir, miniprot_output_path, lineage_filepath,
                                                            self.config)
        if os.path.exists(miniprot_alignment_parser.completeness_output_file):
            os.remove(miniprot_alignment_parser.completeness_output_file)
        miniprot_alignment_parser.Run()
        if self.autolineage:
            marker_genes_filapath = miniprot_alignment_parser.marker_gene_path
            best_match_lineage = self.lineage_searcher.Run(marker_genes_filapath)
            logger.info("best_match_lineage: {}".format(best_match_lineage))
            if best_match_lineage == lineage:
                return
            self.downloader.download_lineage(best_match_lineage)
            lineage = best_match_lineage
            lineage_filepath = os.path.join(self.downloader.lineage_description[lineage][3], "refseq_db.faa.gz")
            output_dir = os.path.join(self.output_folder, lineage)
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            miniprot_output_path = self.miniprot_runner.run_miniprot(self.assembly_path, lineage_filepath, output_dir)
            miniprot_alignment_parser = MiniprotAlignmentParser(output_dir, miniprot_output_path, lineage_filepath,
                                                                self.config)
            miniprot_alignment_parser.Run()
        if os.path.exists("logs"):
            shutil.move("logs", os.path.join(self.output_folder))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output_dir", type=str, help="Run miniprot folder", required=True)
    parser.add_argument("-a", "--assembly_path", type=str, help="Assembly file path", required=True)
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to use")
    parser.add_argument("-l", "--lineage", type=str, help="Lineage file path", default=None)
    parser.add_argument("--library_path", type=str, help="Library path", default="downloads")
    parser.add_argument("--autolineage", action="store_true", help="Auto lineage")

    parser.add_argument("--min_diff",
                        help="The thresholds for the best matching and second best matching. (1st-2nd)/2nd >= d, [0, 1]",
                        type=float, default=0.2)
    parser.add_argument("--min_identity", help="The identity threshold for valid mapping results. [0, 1]", type=float,
                        default=0.4)
    parser.add_argument("--min_length_percent",
                        help="The protein sequence length threshold for valid mapping results. (mapped_gene_length/full_gene_length)>=l, [0, 1]",
                        type=float, default=0.6)
    parser.add_argument("--min_complete",
                        help="The length threshold for complete gene. (mapped_gene_length/full_gene_length)>=c, [0, 1]",
                        type=float, default=0.9)
    parser.add_argument("--min_rise",
                        help="Minimum length threshold to make dupicate take precedence over single or fragmented over single/duplicate. l1>=l2*(1+s), [0, 1]",
                        type=float, default=0.5)

    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    if args.autolineage is False:
        assert args.lineage is not None, "lineage name is required when auto is False! e.g. -l eukaryota"

    buscoprot_runner = BuscoprotRunner(args)
    buscoprot_runner.Run()


if __name__ == "__main__":
    main()
