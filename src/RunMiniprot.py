import os.path
import subprocess
import shlex



class MiniprotRunner:
    def __init__(self, miniprot_execute_command, config):
        self.miniprot_execute_command = miniprot_execute_command
        self.threads = config.threads
        self.autolineage = config.autolineage

    def run_miniprot(self, assembly_filepath, lineage_filepath, output_dir):
        output_filepath = os.path.join(output_dir, "miniprot_output.gff")
        fout = open(output_filepath, "w")
        if self.autolineage:
            miniprot_process = subprocess.Popen(shlex.split(
                "{} -I --outs=0.95 -t {} --aln --gff {} {}".format(self.miniprot_execute_command, self.threads,
                                                                   assembly_filepath, lineage_filepath,
                                                                   output_filepath)), stdout=fout, bufsize=8388608)
        else:
            miniprot_process = subprocess.Popen(shlex.split(
                "{} -I --outs=0.95 -t {} --gff {} {}".format(self.miniprot_execute_command, self.threads,
                                                             assembly_filepath, lineage_filepath, output_filepath)),
                stdout=fout, bufsize=8388608)
        miniprot_process.wait()
        fout.close()
        return output_filepath


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-e", "--execute", type=str, help="Miniprot execute command", required=True)
    parser.add_argument("-r", "--output_dir", type=str, help="Run miniprot folder", required=True)
    parser.add_argument("-a", "--assembly", type=str, help="Assembly file path", required=True)
    parser.add_argument("-l", "--lineage", type=str, help="Lineage file path", required=True)
    parser.add_argument("-t", "--threads", type=int, default=16, help="Number of threads to use")
    parser.add_argument("--autolineage", action="store_true", help="Use autolineage")

    args = parser.parse_args()

    miniprot_execute_command = args.execute
    miniprot_runner = MiniprotRunner(miniprot_execute_command, args)
    miniprot_runner.run_miniprot(args.assembly, args.lineage, args.output_dir)
