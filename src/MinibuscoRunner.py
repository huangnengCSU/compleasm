import os.path
import argparse
import shutil
import time

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


# from .utils import MinibuscoLogger

# logger = MinibuscoLogger(__name__).getlog()

class MinibuscoRunner:
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
        output_dir = os.path.join(self.output_folder, lineage)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        run_miniprot_start_time = time.time()
        miniprot_output_path = self.miniprot_runner.run_miniprot(self.assembly_path, lineage_filepath, output_dir)
        run_miniprot_end_time = time.time()
        analysis_miniprot_start_time = time.time()
        miniprot_alignment_parser = MiniprotAlignmentParser(output_dir, miniprot_output_path, lineage_filepath,
                                                            self.config)
        if os.path.exists(miniprot_alignment_parser.completeness_output_file):
            os.remove(miniprot_alignment_parser.completeness_output_file)
        miniprot_alignment_parser.Run()
        analysis_miniprot_end_time = time.time()
        if self.autolineage:
            autolineage_start_time = time.time()
            marker_genes_filapath = miniprot_alignment_parser.marker_gene_path
            best_match_lineage = self.lineage_searcher.Run(marker_genes_filapath)
            print("best_match_lineage: {}".format(best_match_lineage))
            autolineage_end_time = time.time()
            if best_match_lineage == lineage:
                return
            self.downloader.download_lineage(best_match_lineage)
            lineage = best_match_lineage
            lineage_filepath = os.path.join(self.downloader.lineage_description[lineage][3], "refseq_db.faa.gz")
            output_dir = os.path.join(self.output_folder, lineage)
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            second_run_miniprot_start_time = time.time()
            miniprot_output_path = self.miniprot_runner.run_miniprot(self.assembly_path, lineage_filepath, output_dir)
            second_run_miniprot_end_time = time.time()
            second_analysis_miniprot_start_time = time.time()
            miniprot_alignment_parser = MiniprotAlignmentParser(output_dir, miniprot_output_path, lineage_filepath,
                                                                self.config)
            miniprot_alignment_parser.Run()
            second_analysis_miniprot_end_time = time.time()
        if os.path.exists("logs"):
            shutil.move("logs", os.path.join(self.output_folder))
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


def run(args):
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    if args.autolineage is False:
        assert args.lineage is not None, "lineage name is required when auto is False! e.g. -l eukaryota"

    minibusco_runner = MinibuscoRunner(args)
    minibusco_runner.Run()


def download(args):
    dest = args.destination
    lineage = args.lineage + "_odb10"

    if not os.path.exists(dest):
        os.mkdir(dest)

    d = Downloader(dest)
    d.download_lineage(lineage)


def analysis(args):
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    full_table_output_file = os.path.join(args.output_dir, "full_table.tsv")
    completeness_output_file = os.path.join(args.output_dir, "gene_completeness.tsv")

    if args.path_to_proteins is not None:
        MiniprotAlignmentParser.Local_Run(gff_file=args.gff,
                                          full_table_output_file=full_table_output_file,
                                          completeness_output_file=completeness_output_file,
                                          min_diff=args.min_diff,
                                          min_identity=args.min_identity,
                                          min_length_percent=args.min_length_percent,
                                          min_complete=args.min_complete,
                                          min_rise=args.min_rise,
                                          lineage_file=args.path_to_proteins,
                                          specified_contigs=args.specified_contigs)
    else:
        MiniprotAlignmentParser.Local_Run(gff_file=args.gff,
                                          full_table_output_file=full_table_output_file,
                                          completeness_output_file=completeness_output_file,
                                          min_diff=args.min_diff,
                                          min_identity=args.min_identity,
                                          min_length_percent=args.min_length_percent,
                                          min_complete=args.min_complete,
                                          min_rise=args.min_rise,
                                          specified_contigs=args.specified_contigs)


def main():
    parser = argparse.ArgumentParser()
    subparser = parser.add_subparsers(dest="command", help="Minibusco modules help", required=True)

    run_parser = subparser.add_parser("run",
                                      help="Run minibusco including miniprot alignment and completeness evaluation")
    run_parser.add_argument("-a", "--assembly_path", type=str, help="Input genome file in FASTA format.", required=True)
    run_parser.add_argument("-o", "--output_dir", type=str, help="The output folder.", required=True)
    run_parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to use")
    run_parser.add_argument("-l", "--lineage", type=str, help="Specify the name of the BUSCO lineage to be used. (e.g. eukaryota, primates, saccharomycetes etc.)", default=None)
    run_parser.add_argument("--library_path", type=str,
                            help="Folder path to download lineages or already downloaded lineages. "
                                 "If not specified, a folder named \"downloads\" will be created on the current running path by default to store the downloaded lineage files.",
                            default="downloads")
    run_parser.add_argument("--autolineage", action="store_true", help="Automatically search for the best matching lineage without specifying lineage file.")
    run_parser.add_argument("--specified_contigs", help="Specify the contigs to be evaluated, e.g. chr1 chr2 chr3. If not specified, all contigs will be evaluated.", type=str, nargs='+', default=None)
    run_parser.add_argument("--min_diff",
                            help="The thresholds for the best matching and second best matching.",
                            type=float, default=0.2)
    run_parser.add_argument("--min_identity", help="The identity threshold for valid mapping results.",
                            type=float,
                            default=0.4)
    run_parser.add_argument("--min_length_percent",
                            help="The fraction of protein for valid mapping results.",
                            type=float, default=0.6)
    run_parser.add_argument("--min_complete",
                            help="The length threshold for complete gene.",
                            type=float, default=0.9)
    run_parser.add_argument("--min_rise",
                            help="Minimum length threshold to make dupicate take precedence over single or fragmented over single/duplicate.",
                            type=float, default=0.5)
    run_parser.set_defaults(func=run)

    download_parser = subparser.add_parser("download", help="Download specified BUSCO lineage")
    download_parser.add_argument("-l", "--lineage", type=str,
                                 help="Specify the name of the BUSCO lineage to be downloaded. (e.g. eukaryota, primates, saccharomycetes etc.)", required=True)
    download_parser.add_argument("-d", "--destination", type=str, help="Folder path to download folder", required=True)
    download_parser.set_defaults(func=download)

    analysis_parser = subparser.add_parser("analysis", help="Analysis with miniprot result")
    analysis_parser.add_argument("-g", "--gff", type=str, help="Miniprot output gff file", required=True)
    analysis_parser.add_argument("-o", "--output_dir", type=str, help="Output analysis folder", required=True)
    analysis_parser.add_argument("-p", "--path_to_proteins", type=str, help="Path to protein sequence file in lineage for count the number of total genes", default=None)
    analysis_parser.add_argument("--specified_contigs", help="Specify the contigs to be evaluated, e.g. chr1 chr2 chr3. If not specified, all contigs will be evaluated.",
                            type=str, nargs='+', default=None)
    analysis_parser.add_argument("--min_diff",
                                 help="The thresholds for the best matching and second best matching.",
                                 type=float, default=0.2)
    analysis_parser.add_argument("--min_identity", help="The identity threshold for valid mapping results. [0, 1]",
                                 type=float,
                                 default=0.4)
    analysis_parser.add_argument("--min_length_percent",
                                 help="The fraction of protein for valid mapping results.",
                                 type=float, default=0.6)
    analysis_parser.add_argument("--min_complete",
                                 help="The length threshold for complete gene.",
                                 type=float, default=0.9)
    analysis_parser.add_argument("--min_rise",
                                 help="Minimum length threshold to make dupicate take precedence over single or fragmented over single/duplicate.",
                                 type=float, default=0.5)
    analysis_parser.set_defaults(func=analysis)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
