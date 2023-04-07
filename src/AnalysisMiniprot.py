import argparse
import os
import gzip
import re

import numpy as np
import pandas as pd
from enum import Enum
from Bio import SeqIO
from collections import defaultdict
# from .utils import MinibuscoLogger

# logger = MinibuscoLogger(__name__).getlog()

AminoAcid = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
             "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]




class GeneLabel(Enum):
    Single = 1
    Duplicate = 2
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
            for record in SeqIO.parse(f, "fasta"):
                protein_seqs[record.id] = str(record.seq)
    return protein_seqs


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


class MiniprotAlignmentParser:
    def __init__(self, run_folder, gff_file, lineage_file, config):
        self.autolineage = config.autolineage
        self.run_folder = run_folder  # output_dir/lineage/
        self.completeness_output_file = os.path.join(self.run_folder, "..", "gene_completeness.tsv")
        self.full_table_output_file = os.path.join(self.run_folder, "full_table.tsv")
        self.gff_file = gff_file
        self.lineage_file = lineage_file
        self.min_length_percent = config.min_length_percent
        self.min_diff = config.min_diff
        self.min_identity = config.min_identity
        self.min_complete = config.min_complete
        self.min_rise = config.min_rise
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
                        # items.codons.append([codon_start, codon_end, codon_strand])

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
                    output.gene_label = GeneLabel.Duplicate
                    output.data_record = dataframe.iloc[0]
                    return output
                regions = [(x[1], x[2]) for x in complete_regions]
                clusters = get_region_clusters(regions)
                if len(clusters) == 1:
                    output.gene_label = GeneLabel.Single
                    output.data_record = dataframe.iloc[0]
                    return output
                else:
                    output.gene_label = GeneLabel.Duplicate
                    output.data_record = dataframe.iloc[0]
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
            elif label_length.keys() == {GeneLabel.Duplicate}:
                output.gene_label = GeneLabel.Duplicate
                output.data_record = dataframe_1st.iloc[0]
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
            elif label_length.keys() == {GeneLabel.Single, GeneLabel.Duplicate}:
                if label_length[GeneLabel.Duplicate][0] > label_length[GeneLabel.Single][0] * (1 + min_rise):
                    output.gene_label = GeneLabel.Duplicate
                    if out1.gene_label == GeneLabel.Duplicate:
                        output.data_record = dataframe_1st.iloc[0]
                    elif out2.gene_label == GeneLabel.Duplicate:
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
            elif label_length.keys() == {GeneLabel.Fragmented, GeneLabel.Duplicate}:
                if label_length[GeneLabel.Fragmented][0] > label_length[GeneLabel.Duplicate][0] * (1 + min_rise):
                    output.gene_label = GeneLabel.Fragmented
                    if out1.gene_label == GeneLabel.Fragmented:
                        output.data_record = dataframe_1st.iloc[0]
                    elif out2.gene_label == GeneLabel.Fragmented:
                        output.data_record = dataframe_2nd.iloc[0]
                    else:
                        raise ValueError
                    return output
                else:
                    output.gene_label = GeneLabel.Duplicate
                    if out1.gene_label == GeneLabel.Duplicate:
                        output.data_record = dataframe_1st.iloc[0]
                    elif out2.gene_label == GeneLabel.Duplicate:
                        output.data_record = dataframe_2nd.iloc[0]
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
        protein_file = self.lineage_file
        protein_seqs = load_protein_seqs(protein_file)
        protein_names = defaultdict(int)
        for pid in protein_seqs.keys():
            pid = pid.split("_")[0]
            protein_names[pid] += 1
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
        grouped_df = records_df.groupby(["Target_species"])
        full_table_writer = open(self.full_table_output_file, "w")
        full_table_writer.write(
            "Gene\tStatus\tSequence\tGene Start\tGene End\tStrand\tScore\tLength\tIdentity\tFraction\tFrameshift events\tBest gene\tCodons\n")
        for gene_id in all_species:
            mapped_records = grouped_df.get_group(gene_id)
            mapped_records = mapped_records.sort_values(by=["I+L"], ascending=False)
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
                lower_bound = mean_v - 1.5*std_v
                upper_bound = mean_v + 1.5*std_v
                pass_tids = length_df[(length_df["Protein_length"]>lower_bound) &
                                      (length_df["Protein_length"]<upper_bound)].index.tolist()

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
                full_table_writer.write("{}\t{}\n".format(gene_id, output.gene_label))
            else:
                assert gene_id == output.data_record["Target_species"]
                full_table_writer.write(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(output.data_record["Target_species"],
                                                                  output.gene_label,
                                                                  output.data_record["Contig_id"],
                                                                  output.data_record["Start"],
                                                                  output.data_record["Stop"],
                                                                  output.data_record["Strand"],
                                                                  output.data_record["Score"],
                                                                  output.data_record["Protein_length"],
                                                                  output.data_record["Identity"],
                                                                  output.data_record["Protein_mapped_rate"],
                                                                  output.data_record["Frameshift_events"],
                                                                  output.data_record["Target_id"],
                                                                  output.data_record["Codons"]))

            if output.gene_label == GeneLabel.Single:
                single_genes.append(gene_id)
            elif output.gene_label == GeneLabel.Duplicate:
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

        total_genes = len(protein_names.keys())
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
            fout.write("## lineage: {}\n".format(os.path.dirname(self.lineage_file).split("/")[-1]))
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

    @staticmethod
    def Local_Run(gff_file, full_table_output_file, completeness_output_file, min_diff, min_identity,
                  min_complete, min_rise, min_length_percent, lineage_file=None, autolineage=False):
        single_genes = []
        duplicate_genes = []
        fragmented_genes = []
        interspaced_genes = []
        missing_genes = []
        records = []
        try:
            reader = iter(MiniprotAlignmentParser.parse_miniprot_records(gff_file, autolineage))
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
        grouped_df = records_df.groupby(["Target_species"])
        full_table_writer = open(full_table_output_file, "w")
        full_table_writer.write(
            "Gene\tStatus\tSequence\tGene Start\tGene End\tStrand\tScore\tLength\tIdentity\tFraction\tFrameshift events\tCodons\n")
        for gene_id in all_species:
            mapped_records = grouped_df.get_group(gene_id)
            mapped_records = mapped_records.sort_values(by=["I+L"], ascending=False)
            pass_tids = mapped_records[(mapped_records["Protein_mapped_rate"] >= min_length_percent) | (
                    mapped_records["Identity"] >= min_identity)]["Target_id"].unique()

            # mean std filter
            if len(pass_tids) >= 5:
                lengths_dict = {}
                for tid in pass_tids:
                    lengths_dict[tid] = mapped_records[mapped_records["Target_id"] == tid].iloc[0]["Protein_length"]
                length_df = pd.DataFrame.from_dict(lengths_dict, orient="index", columns=["Protein_length"])
                mean_v = length_df["Protein_length"].mean()
                std_v = length_df["Protein_length"].std()
                lower_bound = mean_v - 1.5*std_v
                upper_bound = mean_v + 1.5*std_v
                pass_tids = length_df[(length_df["Protein_length"]>lower_bound) &
                                      (length_df["Protein_length"]<upper_bound)].index.tolist()

            if len(pass_tids)>0:
                mapped_records = mapped_records[mapped_records["Target_id"].isin(pass_tids)]
                output = MiniprotAlignmentParser.Ost_eval(mapped_records, min_diff, min_identity, min_complete,
                                                          min_rise)
                if output.gene_label == GeneLabel.Fragmented:
                    output = MiniprotAlignmentParser.refine_fragmented(mapped_records)
            else:
                output = OutputFormat()
                output.gene_label = GeneLabel.Missing

            if output.gene_label == GeneLabel.Missing:
                full_table_writer.write("{}\t{}\n".format(gene_id, output.gene_label))
            else:
                assert gene_id == output.data_record["Target_species"]
                full_table_writer.write(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(output.data_record["Target_species"],
                                                                  output.gene_label,
                                                                  output.data_record["Contig_id"],
                                                                  output.data_record["Start"],
                                                                  output.data_record["Stop"],
                                                                  output.data_record["Strand"],
                                                                  output.data_record["Score"],
                                                                  output.data_record["Protein_length"],
                                                                  output.data_record["Identity"],
                                                                  output.data_record["Protein_mapped_rate"],
                                                                  output.data_record["Frameshift_events"],
                                                                  output.data_record["Target_id"],
                                                                  output.data_record["Codons"]))

            if output.gene_label == GeneLabel.Single:
                single_genes.append(gene_id)
            elif output.gene_label == GeneLabel.Duplicate:
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

        if lineage_file is not None:
            protein_seqs = load_protein_seqs(lineage_file)
            protein_names = defaultdict(int)
            for pid in protein_seqs.keys():
                pid = pid.split("_")[0]
                protein_names[pid] += 1
            total_genes = len(protein_names.keys())
        else:
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
        with open(completeness_output_file, 'a') as fout:
            fout.write("## lineage: xx_xx\n")
            fout.write("S:{:.2f}%, {}\n".format(len(single_genes) / total_genes * 100, len(single_genes)))
            fout.write("D:{:.2f}%, {}\n".format(len(duplicate_genes) / total_genes * 100, len(duplicate_genes)))
            fout.write("F:{:.2f}%, {}\n".format(len(fragmented_genes) / total_genes * 100, len(fragmented_genes)))
            fout.write("I:{:.2f}%, {}\n".format(len(interspaced_genes) / total_genes * 100, len(interspaced_genes)))
            fout.write(
                "M:{:.2f}%, {}\n".format((len(missing_genes) + d) / total_genes * 100, len(missing_genes) + d))
            fout.write("N:{}\n".format(total_genes))


if __name__ == "__main__":
    parser_a = argparse.ArgumentParser()
    parser_a.add_argument("--run_folder", dest="run_folder", help="Path to running folder", type=str)
    parser_a.add_argument("--lineage_file", "--lineage_file", help="Gene library file", type=str)
    parser_a.add_argument("-g", "--gff", help="GFF file", required=True)
    parser_a.add_argument("--complete_file", help="Complete file", type=str)
    parser_a.add_argument("--full_table_file", help="Full table file", type=str)
    parser_a.add_argument("-d", "--min_diff",
                          help="The thresholds for the best matching and second best matching. (1st-2nd)/2nd >= d, [0, 1]",
                          type=float, default=0.2)
    parser_a.add_argument("-i", "--min_identity", help="The identity threshold for valid mapping results. [0, 1]",
                          type=float, default=0.4)
    parser_a.add_argument("-l", "--min_length_percent",
                          help="The protein sequence length threshold for valid mapping results. (mapped_gene_length/full_gene_length)>=l, [0, 1]",
                          type=float, default=0.6)
    parser_a.add_argument("-c", "--min_complete",
                          help="The length threshold for complete gene. (mapped_gene_length/full_gene_length)>=c, [0, 1]",
                          type=float, default=0.9)
    parser_a.add_argument("-s", "--min_rise",
                          help="Minimum length threshold to make dupicate take precedence over single or fragmented over single/duplicate. l1>=l2*(1+s), [0, 1]",
                          type=float, default=0.5)
    args = parser_a.parse_args()

    # miniprot_alignment_parser = MiniprotAlignmentParser(args.run_folder, args.gff, args.lineage_file, args)
    # miniprot_alignment_parser.Run()
    MiniprotAlignmentParser.Local_Run(gff_file=args.gff,
                                      full_table_output_file=args.full_table_file,
                                      completeness_output_file=args.complete_file,
                                      min_diff=args.min_diff,
                                      min_identity=args.min_identity,
                                      min_length_percent=args.min_length_percent,
                                      min_complete=args.min_complete,
                                      min_rise=args.min_rise)
