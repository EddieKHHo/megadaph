import csv
from uuid import uuid4

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from click.testing import CliRunner
from fmbiopy.bio import rand_dna_seq
from fmbiopy.system import capture_stdout
from plumbum import local
from plumbum.cmd import list_csomes
from pytest import fixture

from filter_mitochondrial_scaffolds import main


@fixture(
    scope="session",
    params=[
        dict(MT=15000, NOT_MT1=150),
        dict(NOT_MT1=900, NOT_MT2=150),
        dict(MT=13000, NUMT=10000),
    ],
    ids=["single_hit", "no_hit", "numt_hit"],
)
def blast_tsv(request, tempdir):
    tmpfile = local.path(str(tempdir)) / (uuid4().hex + ".tsv")
    with tmpfile.open("w") as tsv:
        writer = csv.writer(tsv, delimiter="\t")
        for name, length in request.param.items():
            writer.writerow(
                [
                    "FI-XINB3_mtDNA",
                    name,
                    "99.699",
                    length,
                    "45",
                    "0",
                    "1",
                    "14969",
                    "59042",
                    "74010",
                    "0.0",
                    "27394",
                ]
            )
    return dict(file=tmpfile, params=request.param)


@fixture(scope="module")
def assembly(tempdir):
    tmpfile = local.path(str(tempdir)) / (uuid4().hex + ".fa")

    scaffolds = [
        SeqRecord(Seq(rand_dna_seq(16000), IUPAC.ambiguous_dna), id="MT"),
        SeqRecord(Seq(rand_dna_seq(20000), IUPAC.ambiguous_dna), id="NOT_MT1"),
        SeqRecord(Seq(rand_dna_seq(20000), IUPAC.ambiguous_dna), id="NOT_MT2"),
        SeqRecord(Seq(rand_dna_seq(50000), IUPAC.ambiguous_dna), id="NUMT"),
    ]

    with tmpfile.open("w") as f:
        for scaffold in scaffolds:
            SeqIO.write(scaffold, f, "fasta")

    return tmpfile


def test_filter_mitochondrial_scaffolds(blast_tsv, assembly, tempdir):
    outfasta = tempdir / (uuid4().hex + ".fa")
    outtxt = tempdir / (uuid4().hex + ".txt")
    runner = CliRunner()
    runner.invoke(
        main,
        [
            "--blast",
            str(blast_tsv["file"]),
            "--assembly",
            str(assembly),
            "--outfasta",
            str(outfasta),
            "--outtxt",
            str(outtxt),
        ],
    )
    output_scaffolds = capture_stdout(list_csomes[outfasta])

    with outtxt.open("r") as f:
        filtered_scaffolds = f.readlines()
        filtered_scaffolds = [x.rstrip() for x in filtered_scaffolds]

    for seqid in ["NOT_MT1", "NOT_MT2", "NUMT"]:
        assert seqid in output_scaffolds
        assert seqid not in filtered_scaffolds
    if "MT" in blast_tsv["params"].keys():
        assert "MT" not in output_scaffolds
        assert "MT" in filtered_scaffolds
    else:
        assert "MT" in output_scaffolds
        assert "MT" not in filtered_scaffolds
