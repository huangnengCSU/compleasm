import os.path
import subprocess
import shlex


class MiniprotRunner:
    def __init__(self, miniprot_execute_command, run_folder, config):
        self.miniprot_execute_command = miniprot_execute_command
        self.run_folder = run_folder
        self.threads = config.threads

    def run_miniprot(self, assembly_filepath, lineage_filepath, output_filename):
        output_filepath = os.path.join(self.run_folder, output_filename)
        fout = open(output_filepath, "w")
        miniprot_process = subprocess.Popen(shlex.split(
            "{} -I --outs=0.95 -t {} --gff {} {}".format(self.miniprot_execute_command, self.threads,
                            assembly_filepath, lineage_filepath, output_filepath)), stdout=fout, bufsize=8388608)
        miniprot_process.wait()
        return output_filepath


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-e", "--execute", type=str, help="Miniprot execute command", required=True)
    parser.add_argument("-r", "--run_folder", type=str, help="Run folder", required=True)
    parser.add_argument("-a", "--assembly", type=str, help="Assembly file path", required=True)
    parser.add_argument("-l", "--lineage", type=str, help="Lineage file path", required=True)
    parser.add_argument("-t", "--threads", type=int, default=16, help="Number of threads to use")

    args = parser.parse_args()

    miniprot_execute_command = args.execute
    run_folder = args.run_folder
    miniprot_runner = MiniprotRunner(miniprot_execute_command, run_folder, args)
    miniprot_runner.run_miniprot(args.assembly, args.lineage, "test.gff")
