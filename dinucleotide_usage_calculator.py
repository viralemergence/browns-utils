from argparse import ArgumentParser
from collections import Counter, defaultdict
from itertools import product
from pathlib import Path
from csv import reader
from typing import Iterator

class SingleSequenceDinucleotideBiasCalculator:
    def __init__(self) -> None:#rna add true and 
        pass

    def run(self, nucleotide_sequence: str) -> None:
        nucleotides = self.set_nucleotides()
        dinucleotides = self.set_dinucleotides(nucleotides)

        nucleotide_counts = self.count_nucleotides(nucleotide_sequence)
        total_nucleotides = sum(nucleotide_counts.values())
        dinucleotide_counts = self.count_observed_dinucleotides(nucleotide_sequence)

        observed_frequencies = self.calculate_observed_frequencies(dinucleotide_counts, total_nucleotides - 1)
        expected_frequencies = self.calculate_expected_frequencies(nucleotide_counts, total_nucleotides, dinucleotides)
        oe_ratios = self.calculate_observed_expected_ratio(observed_frequencies, expected_frequencies, dinucleotides)

        observed_frequencies = self.round_dictionary(observed_frequencies, 4)
        expected_frequencies = self.round_dictionary(expected_frequencies, 4)
        oe_ratios = self.round_dictionary(oe_ratios, 4)

        print(observed_frequencies)
        print()
        print(expected_frequencies)
        print()
        print(oe_ratios)

    @staticmethod
    def set_nucleotides() -> list[str]:
        return ['A', 'T', 'C', 'G']
    
    @staticmethod
    def set_dinucleotides(nucleotides: list[str]) -> list[str]:
        return [''.join(p) for p in product(nucleotides, repeat=2)]

    @staticmethod
    def count_nucleotides(sequence: str) -> defaultdict[int]:
        return defaultdict(int, Counter(sequence))

    @staticmethod
    def count_observed_dinucleotides(sequence: str) -> defaultdict[int]:
        return defaultdict(int, Counter(sequence[i:i+2] for i in range(len(sequence) - 1)))

    @staticmethod
    def calculate_observed_frequencies(dinucleotide_counts: dict[str, int], total_dinucleotides: int) -> defaultdict[float]:
        return defaultdict(int, {dinuc: count / total_dinucleotides for dinuc, count in dinucleotide_counts.items()})
    
    @staticmethod
    def calculate_expected_frequencies(nucleotide_counts: dict[int], total_nucleotides: int, dinucleotides: list[str]) -> dict[float]:
        nucleotide_frequencies = defaultdict(int, {n: nucleotide_counts[n] / total_nucleotides for n in nucleotide_counts})
        return {dinuc: nucleotide_frequencies[dinuc[0]] * nucleotide_frequencies[dinuc[1]] for dinuc in dinucleotides}

    @classmethod
    def calculate_observed_expected_ratio(cls, observed_frequencies: defaultdict[float], expected_frequencies: dict[float], dinucleotides: list[str]) -> dict[float]:
        return {dinuc: cls.divide_with_zero_denominator(observed_frequencies[dinuc], expected_frequencies[dinuc]) for dinuc in dinucleotides}

    @staticmethod
    def divide_with_zero_denominator(numerator: float, denominator: float) -> float:
        try:
            return numerator / denominator
        except ZeroDivisionError:
            return 0

    @staticmethod
    def round_dictionary(dictionary: dict[float], decimals: int) -> dict[float]:
        return {key: round(value, decimals) for key, value in dictionary.items()}
    
    @staticmethod
    def remove_n_nucleotide_keys(dictionary: dict[float], default_dictionary: bool) -> dict:
        cleaned_dictionary = {}
        for key, value in dictionary.items():
            if "N" in key:
                continue
            cleaned_dictionary[key] = value

        if default_dictionary:
            return defaultdict(int, cleaned_dictionary)
        else:
            return cleaned_dictionary

class CdsSequenceDinucleotideBiasCalculator(SingleSequenceDinucleotideBiasCalculator):
    def __init__(self, cds_path: Path, header_type: int) -> None:
        self.cds_path = cds_path
        self.header_type = header_type

    def run(self) -> dict[float]:
        print(f"Starting on: {self.cds_path.stem}")
        unique_ids = self.extract_ids_for_non_redundant_sequences(self.cds_path, self.header_type)

        nucleotides = self.set_nucleotides()
        dinucleotides = self.set_dinucleotides(nucleotides)

        all_nucleotide_counts = defaultdict(int)
        all_dinucleotide_counts = defaultdict(int)

        for fasta_feature in self.fasta_chunker(self.cds_path):
            fasta_name = fasta_feature[0]
            if fasta_name not in unique_ids:
                continue
            fasta_seq = "".join(fasta_feature[1:])

            nucleotide_counts = self.count_nucleotides(fasta_seq)
            dinucleotide_counts = self.count_observed_dinucleotides(fasta_seq)

            all_nucleotide_counts = self.combine_dictionary_values(all_nucleotide_counts, nucleotide_counts, True)
            all_dinucleotide_counts = self.combine_dictionary_values(all_dinucleotide_counts, dinucleotide_counts, True)

        all_nucleotide_counts = self.remove_n_nucleotide_keys(all_nucleotide_counts, True)
        all_dinucleotide_counts = self.remove_n_nucleotide_keys(all_dinucleotide_counts, True)

        total_nucleotides = sum(all_nucleotide_counts.values())

        observed_frequencies = self.calculate_observed_frequencies(all_dinucleotide_counts, total_nucleotides - 1)
        expected_frequencies = self.calculate_expected_frequencies(all_nucleotide_counts, total_nucleotides, dinucleotides)
        oe_ratios = self.calculate_observed_expected_ratio(observed_frequencies, expected_frequencies, dinucleotides)

        observed_frequencies = self.round_dictionary(observed_frequencies, 4)
        expected_frequencies = self.round_dictionary(expected_frequencies, 4)
        oe_ratios = self.round_dictionary(oe_ratios, 4)

        print(observed_frequencies)
        print()
        print(expected_frequencies)
        print()
        print(oe_ratios)

    @classmethod
    def extract_ids_for_non_redundant_sequences(cls, fasta_path: Path, header_type: int) -> set[str]:
        sequence_info = defaultdict(lambda: defaultdict(int))
        for fasta_feature in cls.fasta_chunker(fasta_path):
            fasta_name = fasta_feature[0]
            if header_type == 0:
                # >lcl|NC_014373.1_cds_YP_003815432.1_1 [gene=NP] [locus_tag=BDBVp1] [db_xref=GeneID:9487269] [protein=nucleoprotein] [protein_id=YP_003815432.1] [location=458..2677] [gbkey=CDS]
                gene = fasta_name.split("[")[1].strip().replace("gene=", "")[:-1] # NOTE: Should clean up
            if header_type == 1:
                # >NC_055939.1:85..3213 |glycoprotein [Scaldis River bee virus]
                gene = fasta_name.split("|")[1].split("[")[0].strip()
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

    @staticmethod
    def combine_dictionary_values(dictionary_1: dict[float], dictionary_2: dict[float], default_dictionary: bool) -> dict:
        combined_counter = Counter(dictionary_1) + Counter(dictionary_2)
        if default_dictionary:
            return defaultdict(int, combined_counter)
        else:
            return dict(combined_counter)

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-c", "--cds", type=str, required=True)
    parser.add_argument("-header", type=int, required=False)
    args = parser.parse_args()

    #dna_sequence = "ATCGATCGTAGCTAGCTAGCTA"
    # dna_sequence = "ATATATATA"
    # ssdbc = SingleSequenceDinucleotideBiasCalculator()
    # ssdbc.run(dna_sequence)

    csdbc = CdsSequenceDinucleotideBiasCalculator(Path(args.cds), args.header)
    csdbc.run()