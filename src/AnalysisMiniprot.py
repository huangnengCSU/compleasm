import argparse
import os
import gzip
import pandas as pd
from enum import Enum
from Bio import SeqIO
from collections import defaultdict

AminoAcid = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
             "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "*"]


class GeneLabel(Enum):
    Single = 1
    Duplicate = 2
    Fragmented = 3
    Missing = 4


class MiniprotGffItems:
    def __init__(self):
        self.atn_seq = ""
        self.att_seq = ""
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

    def show(self):
        return [self.atn_seq,
                self.att_seq,
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
                self.codons]

    def print(self):
        print(self.show())


def get_region_clusters(regions):
    # TODO: use contig_id, check whether the algorithm is correct
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


def load_protein_seqs(fasta_file):
    protein_seqs = {}
    if fasta_file.endswith(".gz"):
        with gzip.open(fasta_file, "rt") as f:
            for record in SeqIO.parse(f, "fasta"):
                protein_seqs[record.id] = str(record.seq)
    else:
        with open(fasta_file, "r") as f:
            for record in SeqIO.parse(f, "fasta.gz"):
                protein_seqs[record.id] = str(record.seq)
    return protein_seqs


class OutputFormat:
    def __init__(self):
        self.gene_label = None
        self.single_complete_gene_id = None


class MiniprotAlignmentParser:
    def __init__(self, run_folder, gff_file, lineage_file, config):
        self.run_folder = run_folder
        self.gff_file = gff_file
        self.lineage_file = lineage_file
        self.min_il = config.min_il
        self.min_length_percent = config.min_length_percent
        self.min_diff = config.min_diff
        self.min_identity = config.min_identity
        self.min_complete = config.min_complete
        self.min_rise = config.min_rise
        self.marker_gene_path = os.path.join(self.run_folder, "gene_marker.fasta")

    @staticmethod
    def parse_miniprot_records(gff_file):
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
                    items.contig_stop = int(fields[8])
                    items.score = int(fields[9])
                    # Score = fields[12].strip().split(":")[2]

                    """ # This is for parse the ATN/ATT/AAS/AQA sequences in miniprot --aln
                    atn_line = gff.readline()
                    ata_line = gff.readline()
                    aas_line = gff.readline()
                    aqa_line = gff.readline()
                    items.atn_seq = atn_line.strip().split("\t")[1].replace("-", "")
                    items.ata_seq = ata_line.strip().split("\t")[1]
                    new_ata = []
                    for i in range(len(ata_seq)):
                        if ata_seq[i].upper() not in AminoAcid:
                            continue
                        else:
                            new_ata.append(ata_seq[i])
                    att_seq = "".join(new_ata)
                    items.att_seq = att_seq
                    """
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
                        items.codons.append([codon_start, codon_end, codon_strand])

    @staticmethod
    def record_1st_gene_label(dataframe, min_identity, min_complete):
        # check records with same tid of the best record
        output = OutputFormat()
        gene_id = dataframe.iloc[0]["Target_id"]
        dataframe = dataframe[dataframe["Identity"] >= min_identity]
        if dataframe.shape[0] == 0:
            output.gene_label = GeneLabel.Missing
            output.single_complete_gene_id = None
            return output
        elif dataframe.shape[0] == 1:
            if dataframe.iloc[0]["Protein_mapped_rate"] >= min_complete:
                output.gene_label = GeneLabel.Single
                output.single_complete_gene_id = gene_id
                return output
            else:
                output.gene_label = GeneLabel.Fragmented
                output.single_complete_gene_id = None
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
                output.single_complete_gene_id = None
                return output
            elif len(complete_regions) == 1:
                output.gene_label = GeneLabel.Single
                output.single_complete_gene_id = gene_id
                return output
            else:
                ctgs = [x[0] for x in complete_regions]
                if len(set(ctgs)) > 1:
                    output.gene_label = GeneLabel.Duplicate
                    output.single_complete_gene_id = None
                    return output
                regions = [(x[1], x[2]) for x in complete_regions]
                clusters = get_region_clusters(regions)
                if len(clusters) == 1:
                    output.gene_label = GeneLabel.Single
                    output.single_complete_gene_id = gene_id
                    return output
                else:
                    output.gene_label = GeneLabel.Duplicate
                    output.single_complete_gene_id = None
                    return output

    @staticmethod
    def record_1st_2nd_gene_label(dataframe_1st, dataframe_2nd, min_identity, min_complete, min_rise):
        # check top 1st and 2nd records whether they are the same gene
        output = OutputFormat()
        dataframe_1st = dataframe_1st[dataframe_1st["Identity"] >= min_identity]
        dataframe_2nd = dataframe_2nd[dataframe_2nd["Identity"] >= min_identity]
        if dataframe_1st.shape[0] >= 1 and dataframe_2nd.shape[0] == 0:
            out = MiniprotAlignmentParser.record_1st_gene_label(dataframe_1st, min_identity, min_complete)
            return out
        if dataframe_1st.shape[0] == 0 and dataframe_2nd.shape[0] >= 1:
            out = MiniprotAlignmentParser.record_1st_gene_label(dataframe_2nd, min_identity, min_complete)
            return out
        if dataframe_1st.shape[0] == 0 and dataframe_2nd.shape[0] == 0:
            output.gene_label = GeneLabel.Missing
            output.single_complete_gene_id = None
            return output
        else:
            gene_id1 = dataframe_1st.iloc[0]["Target_id"]
            gene_id2 = dataframe_2nd.iloc[0]["Target_id"]
            protein_length1 = dataframe_1st.iloc[0]["Protein_length"]
            protein_length2 = dataframe_2nd.iloc[0]["Protein_length"]

            label_length = defaultdict(list)
            out1 = MiniprotAlignmentParser.record_1st_gene_label(dataframe_1st, min_identity, min_complete)
            label_length[out1.gene_label].append([protein_length1, out1.single_complete_gene_id])
            out2 = MiniprotAlignmentParser.record_1st_gene_label(dataframe_2nd, min_identity, min_complete)
            label_length[out2.gene_label].append([protein_length2, out2.single_complete_gene_id])
            if label_length.keys() == {GeneLabel.Single}:
                output.gene_label = GeneLabel.Single
                output.single_complete_gene_id = gene_id1
                return output
            elif label_length.keys() == {GeneLabel.Fragmented}:
                output.gene_label = GeneLabel.Fragmented
                output.single_complete_gene_id = None
                return output
            elif label_length.keys() == {GeneLabel.Duplicate}:
                output.gene_label = GeneLabel.Duplicate
                output.single_complete_gene_id = None
                return output
            elif label_length.keys() == {GeneLabel.Single, GeneLabel.Fragmented}:
                if label_length[GeneLabel.Fragmented][0][0] > label_length[GeneLabel.Single][0][0] * (1 + min_rise):
                    output.gene_label = GeneLabel.Fragmented
                    output.single_complete_gene_id = None
                    return output
                else:
                    output.gene_label = GeneLabel.Single
                    output.single_complete_gene_id = label_length[GeneLabel.Single][0][1]
                    return output
            elif label_length.keys() == {GeneLabel.Single, GeneLabel.Duplicate}:
                if label_length[GeneLabel.Duplicate][0][0] > label_length[GeneLabel.Single][0][0] * (1 + min_rise):
                    output.gene_label = GeneLabel.Duplicate
                    output.single_complete_gene_id = None
                    return output
                else:
                    output.gene_label = GeneLabel.Single
                    output.single_complete_gene_id = label_length[GeneLabel.Single][0][1]
                    return output
            elif label_length.keys() == {GeneLabel.Fragmented, GeneLabel.Duplicate}:
                if label_length[GeneLabel.Fragmented][0][0] > label_length[GeneLabel.Duplicate][0][0] * (1 + min_rise):
                    output.gene_label = GeneLabel.Fragmented
                    output.single_complete_gene_id = None
                    return output
                else:
                    output.gene_label = GeneLabel.Duplicate
                    output.single_complete_gene_id = None
                    return output
            else:
                print("Error")
                print(label_length)
                raise ValueError

    @staticmethod
    def Ost_eval(dataframe, difficial_rate, min_identity, min_complete, min_rise):
        if dataframe.shape[0] == 0:
            output = OutputFormat()
            output.gene_label = GeneLabel.Missing
            output.single_complete_gene_id = None
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

    def Run(self):
        single_genes = []
        duplicate_genes = []
        fragmented_genes = []
        missing_genes = []
        records = []
        single_complete_proteins = []
        gff_file = self.gff_file
        protein_file = self.lineage_file
        protein_seqs = load_protein_seqs(protein_file)
        try:
            reader = iter(self.parse_miniprot_records(gff_file))
            for items in reader:
                (Atn_seq, Att_seq, Target_id, Contig_id, Protein_length, Protein_Start, Protein_End, Start, Stop, Strand, Score, Rank, Identity,
                 Positive, Codons) = items.show()
                Target_species = Target_id.split("_")[0]
                # items.print()
                records.append([Target_species, Target_id, Contig_id, Protein_length, Protein_Start, Protein_End,
                                Protein_End - Protein_Start, (Protein_End - Protein_Start) / Protein_length, Start,
                                Stop, Stop - Start, Strand, Rank, Identity, Positive,
                                (Protein_End - Protein_Start) / Protein_length + Identity])
        except StopIteration:
            pass
        records_df = pd.DataFrame(records, columns=["Target_species", "Target_id", "Contig_id", "Protein_length",
                                                    "Protein_Start", "Protein_End", "Protein_mapped_length",
                                                    "Protein_mapped_rate", "Start", "Stop", "Genome_mapped_length",
                                                    "Strand", "Rank", "Identity", "Positive", "I+L"])
        all_species = records_df["Target_species"].unique()
        grouped_df = records_df.groupby(["Target_species"])
        for gene_id in all_species:
            mapped_records = grouped_df.get_group(gene_id)
            mapped_records = mapped_records.sort_values(by=["I+L"], ascending=False)

            if mapped_records.iloc[0]["I+L"] >= self.min_il or \
                    mapped_records.iloc[0]["Protein_mapped_rate"] >= self.min_length_percent:
                output = self.Ost_eval(mapped_records, self.min_diff, self.min_identity, self.min_complete,
                                       self.min_rise)
                if output.gene_label == GeneLabel.Single:
                    single_complete_proteins.append(output.single_complete_gene_id)
            else:
                output = OutputFormat()
                output.gene_label = GeneLabel.Missing
                output.single_complete_gene_id = None

            if output.gene_label == GeneLabel.Single:
                single_genes.append(gene_id)
            elif output.gene_label == GeneLabel.Duplicate:
                duplicate_genes.append(gene_id)
            elif output.gene_label == GeneLabel.Fragmented:
                fragmented_genes.append(gene_id)
            elif output.gene_label == GeneLabel.Missing:
                missing_genes.append(gene_id)
            else:
                print("Error")
                raise ValueError

        d = len(all_species) - len(single_genes) - len(duplicate_genes) - len(fragmented_genes) - len(missing_genes)
        print("S:{:.2f}%, {}".format(len(single_genes) / len(all_species) * 100, len(single_genes)))
        print("D:{:.2f}%, {}".format(len(duplicate_genes) / len(all_species) * 100, len(duplicate_genes)))
        print("F:{:.2f}%, {}".format(len(fragmented_genes) / len(all_species) * 100, len(fragmented_genes)))
        print("M:{:.2f}%, {}".format((len(missing_genes) + d) / len(all_species) * 100, len(missing_genes) + d))
        print("N:{}".format(len(all_species)))
        # print("Duplicate genes:")
        # print(duplicate_genes)
        # print("Fragmented genes:")
        # print(fragmented_genes)
        # print("Missing genes:")
        # print(missing_genes)

        with open(self.marker_gene_path, "w") as fout:
            for protein_id in single_complete_proteins:
                fout.write(">{}\n{}\n".format(protein_id, protein_seqs[protein_id]))


if __name__ == "__main__":
    parser_a = argparse.ArgumentParser(description="Ost_eval")
    parser_a.add_argument("--run_folder", dest="run_folder", help="Path to running folder", type=str, required=True)
    parser_a.add_argument("--lineage_file", "--lineage_file", help="Gene library file", required=True)
    parser_a.add_argument("-g", "--gff", help="GFF file", required=True)
    parser_a.add_argument("-d", "--min_diff",
                          help="The thresholds for the best matching and second best matching. (1st-2nd)/2nd >= d, [0, 1]",
                          type=float, default=0.2)
    parser_a.add_argument("-i", "--min_identity", help="The identity threshold for valid mapping results. [0, 1]",
                          type=float, default=0.6)
    parser_a.add_argument("-l", "--min_length_percent",
                          help="The protein sequence length threshold for valid mapping results. (mapped_gene_length/full_gene_length)>=l, [0, 1]",
                          type=float, default=0.64)
    parser_a.add_argument("-e", "--min_il",
                          help="The thresholds for sum of identity and mapped length for valid mapping results. identity+mapped_rate >= e, [0, 2]",
                          type=float, default=1.4)
    parser_a.add_argument("-c", "--min_complete",
                          help="The length threshold for complete gene. (mapped_gene_length/full_gene_length)>=c, [0, 1]",
                          type=float, default=0.9)
    parser_a.add_argument("-s", "--min_rise",
                          help="Minimum length threshold to make dupicate take precedence over single or fragmented over single/duplicate. l1>=l2*(1+s), [0, 1]",
                          type=float, default=1.5)
    args = parser_a.parse_args()

    miniprot_alignment_parser = MiniprotAlignmentParser(args.run_folder, args.gff, args.lineage_file, args)
    miniprot_alignment_parser.Run()
