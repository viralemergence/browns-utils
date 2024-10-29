from argparse import ArgumentParser
from csv import reader
from pathlib import Path
from typing import Iterator

class FastaSampler:
    def __init__(self, fasta_path: str, sampling_number: int) -> None:
        self.fasta_path = Path(fasta_path)
        self.sampling_number = sampling_number

    def run(self) -> None:
        fasta_outpath = self.set_fasta_outpath(self.fasta_path, self.sampling_number)
        self.head_fasta(self.fasta_path, fasta_outpath, self.sampling_number)

    @staticmethod
    def set_fasta_outpath(fasta_path: Path, sampling_number: int) -> Path:
        return fasta_path.parent / f"{fasta_path.stem}_{sampling_number}_sampled.fasta"

    @classmethod
    def head_fasta(cls, fasta_path: Path, fasta_outpath: Path, sampling_number: int) -> None:
        with fasta_outpath.open("w") as outhandle:
            for i, fasta_seq in enumerate(cls.fasta_chunker(fasta_path), 1):
                if i > sampling_number:
                    break
                for line in fasta_seq:
                    outhandle.write(line + "\n")

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
    
    #fr = FastaRenamer(args.fasta)
    fr = FastaSampler(args.fasta, 10_000)
    fr.run()