from argparse import ArgumentParser
from collections import defaultdict
from csv import reader
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
    def __init__(self, cds_path: Path) -> None:
        self.cds_path = cds_path
        self.codon_dict = self.set_codon_dict()

    @classmethod
    def set_codon_dict(cls) -> dict[int]:
        codons = cls.set_codons()
        return {codon: 0 for codon in codons}

    @staticmethod
    def set_codons() -> list[str]:
        nucleotides = ["A", "T", "C", "G"]
        return [''.join(p) for p in product(nucleotides, repeat=3)]

    def run(self) -> None:
        total = 0
        not_divisible = []
        bad_codons = defaultdict(int)
        for fasta_feature in self.fasta_chunker(self.cds_path):
            total += 1
            fasta_name = fasta_feature[0]
            fasta_seq = "".join(fasta_feature[1:])

            if (remainder := len(fasta_seq) % 3) != 0:
                not_divisible.append(fasta_name)
                continue

            for codon in self.codon_chunker(fasta_seq):
                try:
                    self.codon_dict[codon] += 1
                except KeyError:
                    bad_codons[codon] += 1
        codon_usage = self.convert_codon_counts_to_proportion(self.codon_dict)

        print(total)
        print(len(not_divisible))
        print(self.codon_dict)
        print(bad_codons)
        print(codon_usage)

    @staticmethod
    def fasta_chunker(fasta_path: Path) -> Iterator[list[str]]:
        fasta_seq = []
        first_chunk = True
        with fasta_path.open() as inhandle:
            reader_iterator = reader(inhandle)
            for line in reader_iterator:
                line = line[0]
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

if __name__ == "__main__":
    parser = ArgumentParser()
    # parser.add_argument("-g", "--gtf", type=str, required=True)
    parser.add_argument("-c", "--cds", type=str, required=True)
    args = parser.parse_args()

    cm = CdsManager(Path(args.cds))
    cm.run()