from argparse import ArgumentParser
from collections import defaultdict
from csv import reader, writer
from itertools import product
from pathlib import Path
from typing import Any, Iterator, Union

class GtfManager:
    def __init__(self, gtf_path: Path) -> None:
        self.gtf_path = gtf_path

    def run(self) -> list[dict[Any]]:
        cds_info = list()
        with self.gtf_path.open() as inhandle:
            reader_iterator = reader(inhandle, delimiter="\t")
            for line in reader_iterator:
                if not self.check_if_cds(line):
                    continue
                if (cds_info_singlet := self.extract_cds_info(line)) is None:
                    continue
                cds_info.append(cds_info_singlet)
        return cds_info

    @staticmethod
    def check_if_cds(line: list[str]) -> bool:
        try:
            feature = line[2]
        except IndexError:
            return False

        if feature == "CDS":
            return True
        else:
            return False

    @staticmethod
    def extract_cds_info(line: list[str]) -> Union[None, dict[Any]]:
        chromosome = line[0]
        start = int(line[3]) # 1-based
        end = int(line[4]) # 1-based
        strand = line[6]
        frame = int(line[7])

        attributes = [info.strip() for info in line[-1].split(";")][:-1]
        attribute_pairs = {}
        for attribute in attributes:
            pair = attribute.split('"')
            try:
                attribute_pairs[pair[0].strip()] = pair[1]
            except IndexError:
                return None

        gene_id = attribute_pairs["gene_id"]
        exon_number = int(attribute_pairs["exon_number"])

        keys = ["chromosome", "start", "end", "strand", "frame", "gene_id", "exon_number"]
        cds_info = dict(zip(keys, map(eval, keys)))
        return cds_info

class CdsManager:
    def __init__(self, cds_path: Path, outdir_path: Path) -> None:
        self.cds_path = cds_path
        self.outdir_path = outdir_path
        self.codon_usage_path = self.outdir_path / f"{cds_path.stem}_codon_usage.csv"
        self.codon_log_path = self.outdir_path / f"{cds_path.stem}_codon_usage_stats.csv"

        self.codon_dict = self.set_codon_dict()

        self.total_cds = 0
        self.non_redundant_cds = 0
        self.not_divisible_cds = []
        self.unconforming_codons = defaultdict(int)

    @classmethod
    def set_codon_dict(cls) -> dict[int]:
        codons = cls.set_codons()
        return {codon: 0 for codon in codons}

    @staticmethod
    def set_codons() -> list[str]:
        nucleotides = ["A", "T", "C", "G"]
        return [''.join(p) for p in product(nucleotides, repeat=3)]

    def run(self) -> None:
        print(f"Starting on: {self.cds_path.stem}")
        unique_ids = self.extract_ids_for_non_redundant_sequences(self.cds_path)

        for fasta_feature in self.fasta_chunker(self.cds_path):
            self.total_cds += 1
            fasta_name = fasta_feature[0]
            if fasta_name not in unique_ids:
                continue
            self.non_redundant_cds += 1
            fasta_seq = "".join(fasta_feature[1:])

            if (remainder := len(fasta_seq) % 3) != 0:
                self.not_divisible_cds.append(fasta_name)
                continue

            self.process_cds_codons(fasta_seq, self.codon_dict, self.unconforming_codons)

        codon_usage = self.convert_codon_counts_to_proportion(self.codon_dict)
        self.write_codon_usage(self.codon_usage_path, codon_usage)

        log_info = {"Total_CDS": self.total_cds,
                    "Non-redundant_CDS": self.non_redundant_cds,
                    "Not_divisible_CDS": len(self.not_divisible_cds)}
        self.write_codon_stats(self.codon_log_path, log_info, self.unconforming_codons)

    @classmethod
    def extract_ids_for_non_redundant_sequences(cls, fasta_path: Path) -> set[str]:
        sequence_info = defaultdict(lambda: defaultdict(int))
        for fasta_feature in cls.fasta_chunker(fasta_path):
            fasta_name = fasta_feature[0]
            gene = fasta_name.split("[")[1].strip().replace("gene=", "")[:-1] # NOTE: Should clean up
            fasta_seq_len = len("".join(fasta_feature[1:]))

            lead_seq_length = sequence_info[gene]["length"]
            if fasta_seq_len > lead_seq_length:
                sequence_info[gene]["length"] = fasta_seq_len
                sequence_info[gene]["id"] = fasta_name

        ids_for_non_redundant_sequences = set()
        for gene, info in sequence_info.items():
            unique_id = info["id"]
            ids_for_non_redundant_sequences.add(unique_id)

        return ids_for_non_redundant_sequences

    @staticmethod
    def fasta_chunker(fasta_path: Path) -> Iterator[list[str]]:
        fasta_seq = []
        first_chunk = True
        with fasta_path.open() as inhandle:
            reader_iterator = reader(inhandle)
            for line in reader_iterator:
                try:
                    line = line[0]
                except IndexError:
                    if len(line) == 0:
                        continue
                    else:
                        raise Exception
                if not line.startswith(">"):
                    fasta_seq.append(line)
                else:
                    if first_chunk:
                        fasta_seq.append(line)
                        first_chunk = False
                        continue
                    yield fasta_seq
                    fasta_seq = [line]
            if fasta_seq:
                yield fasta_seq

    @classmethod
    def process_cds_codons(cls, fasta_seq: str, codons: dict[int], unconforming_codons: defaultdict[int]) -> None:
        for codon in cls.codon_chunker(fasta_seq):
            try:
                codons[codon] += 1
            except KeyError:
                unconforming_codons[codon] += 1

    @staticmethod
    def codon_chunker(dna_seq: str) -> Iterator[str]:
        codon_size = 3
        for i in range(0, len(dna_seq), codon_size):
            yield dna_seq[i:i+codon_size]

    @staticmethod
    def convert_codon_counts_to_proportion(codons: dict[int]) -> dict[float]:
        codon_total = sum(codons.values())
        codon_proportions = {codon: round(count/codon_total, 4) for codon, count in codons.items()}
        return codon_proportions

    @staticmethod
    def write_codon_usage(codon_usage_path: Path, codon_usage: dict[float]) -> None:
        with codon_usage_path.open("w") as outhandle:
            csv_writer = writer(outhandle)
            for codon, usage in codon_usage.items():
                csv_writer.writerow([codon, usage])

    @staticmethod
    def write_codon_stats(codon_log_path: Path, log_info: dict[Any], unconforming_codons: defaultdict[int]) -> None:
        with codon_log_path.open("w") as outhandle:
            csv_writer = writer(outhandle)
            for key, value in log_info.items():
                csv_writer.writerow([key, value])
            for codon, count in unconforming_codons.items():
                csv_writer.writerow([codon, count])

class CdsArrayManager:
    def __init__(self) -> None:
        pass

    @classmethod
    def get_cds_file(cls, input_dir: Path, file_index: int) -> Path:
        cds_files = sorted([cds for cds in input_dir.iterdir()])
        return cds_files[file_index]

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-c", "--cds", type=str, required=True)
    parser.add_argument("-o", "--outdir", type=str, required=True)
    parser.add_argument("-j", "--file_index", type=int, required=False)
    args = parser.parse_args()

    if (file_index := args.file_index) is not None:
        cds_file_path = CdsArrayManager().get_cds_file(Path(args.cds), file_index)
    else:
        cds_file_path = Path(args.cds)

    cm = CdsManager(cds_file_path, Path(args.outdir))
    cm.run()