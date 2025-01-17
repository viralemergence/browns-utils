from argparse import ArgumentParser
from csv import reader
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


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-g", "--gtf", type=str, required=True)
    args = parser.parse_args()
    
    gm = GtfManager(Path(args.gtf))
    cds_info = gm.run()