from Bio import SeqIO
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

if __name__ == "__main__":
    fasta_path = ""
    sampling_number = 10_000
    
    fs = FastaSampler(fasta_path, sampling_number)
    fs.run()