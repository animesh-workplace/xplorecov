import string
import random
import pandas
import argparse
from Bio import SeqIO
from pathlib import Path
from mpire import WorkerPool
from snakemake.shell import shell
from more_itertools import chunked
from tempfile import TemporaryDirectory
from mpire.utils import make_single_arguments

parser = argparse.ArgumentParser()
parser.add_argument(
    "--sequences",
    required=True,
    help="FASTA file of sequences to partition into smaller chunks",
)
parser.add_argument(
    "--sequences-per-group",
    required=True,
    type=int,
    help="number of sequences to include in each group",
)
parser.add_argument(
    "--threads",
    required=True,
    type=int,
    help="threads",
)
parser.add_argument(
    "--output",
    required=True,
    help="Output filename",
)

args = parser.parse_args()
sequences = SeqIO.to_dict(SeqIO.parse(args.sequences, "fasta"))
tempdir = TemporaryDirectory()


def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return "".join(random.choice(chars) for _ in range(size))


def run_pango_usher(seq_items):
    seq = list(map(sequences.get, seq_items))
    base_name = Path(f"{tempdir.name}/sequence_{id_generator()}")
    fasta_file = base_name / f"{base_name}.fasta"
    lineage_report = base_name / f"{base_name}.csv"
    base_name.mkdir(parents=True, exist_ok=True)
    SeqIO.write(seq, fasta_file, "fasta")
    print(f"Written successfully {fasta_file}")
    # threads = int(args.threads / 8)
    print(f"Total threads allowed {args.threads}")
    threads = 4
    shell(
        """
            pangolin {fasta_file} --outfile {lineage_report} -t {threads} --analysis-mode fast        
        """
    )
    pango_usher_output = pandas.read_csv(lineage_report, sep=",")
    print(f"Pangolin Usher completed successfully {base_name}")
    return pango_usher_output


with WorkerPool(n_jobs=1) as pool:
    output = pool.map(
        run_pango_usher,
        make_single_arguments(
            list(chunked(sequences.keys(), args.sequences_per_group)), generator=False
        ),
    )

combined_metadata = pandas.DataFrame()
for metadata in output:
    combined_metadata = pandas.concat([combined_metadata, metadata])

combined_metadata.to_csv(args.output, sep="\t", index=False)
