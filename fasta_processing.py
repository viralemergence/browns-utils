from argparse import ArgumentParser
from Bio import SeqIO
from csv import reader
from pathlib import Path

class FastaSampler:
    def __init__(self, fasta_path: str, sampling_number: int) -> None:
        self.fasta_path = Path(fasta_path)
        self.sampling_number = sampling_number

    def run(self) -> None:
        fasta_outpath = self.set_fasta_outpath(self.fasta_path, self.sampling_number)
        self.head_fasta(self.fasta_path, fasta_outpath, self.sampling_number)

    @classmethod
    def head_fasta(cls, fasta_path: Path, fasta_outpath: Path, sampling_number: int) -> None:
        with fasta_path.open() as inhandle, fasta_outpath.open("w") as outhandle:
            fasta_reader = SeqIO.parse(inhandle, "fasta")
            for i, record in enumerate(fasta_reader, 1):
                if i > sampling_number:
                    break
                SeqIO.write(record, outhandle, "fasta")

    @staticmethod
    def set_fasta_outpath(fasta_path: Path, sampling_number: int) -> Path:
        return fasta_path.parent / f"{fasta_path.stem}_{sampling_number}_sampled.fasta"

class FastaRenamer:
    def __init__(self, fasta_path: str) -> None:
        self.fasta_path = Path(fasta_path)
        self.fasta_outpath = self.set_fasta_outpath(self.fasta_path, "_renamed.fasta")

    def run(self) -> None:
        with self.fasta_path.open() as inhandle, self.fasta_outpath.open("w") as outhandle:
            reader_iterator = reader(inhandle)
            i = 0
            for line in reader_iterator:
                if line[0].startswith(">"):
                    i += 1
                    new_name = f">seq_{i}"
                    outhandle.write(new_name + "\n")
                else:
                    outhandle.write(line[0] + "\n")

    @staticmethod
    def set_fasta_outpath(fasta_path: Path, suffix: str) -> Path:
        return fasta_path.parent / f"{fasta_path.stem}{suffix}"

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-f", "--fasta", type=str, required=True)
    args = parser.parse_args()
    
    fr = FastaRenamer(args.fasta)
    fr.run()