#!/usr/bin/env python
"""Megadaph Decontamination Pipeline

Built using Ruffus pipelining infrastructure (See http://www.ruffus.org.uk/)

Pipeline steps are as follows:
    Index Fasta Files -> Align Fastq files ->

Usage:
    megadaph_decontam_pipeline.py [--flowchart=<file.svg>] [--verbose]
        [--just-print] (--threads=<n>) (--assemblies=<directory>)
        (--fwd_reads=<directory>) (--rev_reads=<directory>)
    megadaph_decontam_pipeline.py (-h | --help)


Options:
    -h --help  Show this screen
    -j --just-print  Show actions that pipeline would take if run
    -c --flowchart=<file.svg>  Create a flowchart image describing the pipeline
    -v --verbose  Increased logging verbosity
    -p --threads=<nthreads>  Number of threads to use
    -a --assemblies=ASSEMBLY_DIR  Location of directory holding .fasta assembly
files.
    -x --fwd_reads=FORWARD_READ_DIR  Location of directory holding forward
reads
    -y --rev_reads=REVERSE_READ_DIR  Location of directory holding reverse
reads. Can be the same as x - direction will be distinguished using 1./2. in
filenames.


Required Input Files:
    1. A set of paired Fastq files
    2. A corresponding set of assemblies produced from the Fastq files in
    Fasta format

Notes:
    Matched Fastq and Fasta files must share the same prefix and be named using
    dot filename convention (E.g FA_SC.1.trim.fastq). Forward reads must have
    '.1' immediately after shared prefix, reverse reads must have '.2'.
"""

from pathlib import Path
from typing import Any
from typing import Dict
from typing import List
from typing import Tuple

from docopt import docopt
import ruffus

import fmbiopy.biofile as biofile
import fmbiopy.fmpaths as fmpaths
import fmbiopy.fmruffus as fmruffus


def _get_gzip_status(files: List[Path]) -> bool:
    """Get the gzip status of a group of files

    Files are expected to have the same Gzip status.

    Returns
    -------
    True if all `files` are gzipped, False if all `files` are unzipped.

    Raises
    ------
    ValueError
        If only some files are gzipped
    """
    gzipped = [f.suffix == '.gz' for f in files]
    unzipped = [not x for x in gzipped]
    if all(gzipped):
        return True
    elif all(unzipped):
        return False
    raise ValueError("Only some of the input files are gzipped")


def _validate_args(
        assemblies: List[Path],
        fwd_reads: List[Path],
        rev_reads: List[Path],
        )-> None:
    """Check that the given input files are valid

    Paired fastq files must share the same prefix and use .1/.2 to distinguish
    read direction.

    Parameters
    ----------
    fwd_reads, rev_reads
        The user inputted fastq file paths
    assemblies
        The user inputted assembly files
    """

    inputs = [fwd_reads, rev_reads, assemblies]
    gzipped = [_get_gzip_status(input) for input in inputs]
    biofile.MatchedPrefixGroup([
        biofile.BiofileGroup(fwd_reads, filetype='fastq', gzipped=gzipped[0]),
        biofile.BiofileGroup(rev_reads, filetype='fastq', gzipped=gzipped[1]),
        biofile.BiofileGroup(assemblies, filetype='fasta', gzipped=gzipped[2])
        ]).validate()


def _extract_input_files(
        opt: Dict[str, Any]
        )-> List[Tuple[Path, Path, Path]]:
    """Extract the input files from the given directories

    Input file prefixes are checked to ensure they can be matched correctly and
    have the correct file extensions.

    Parameters
    ----------
    opt
        A dictionary of parsed input parameters formatted as in `docopt`

    Returns
    -------
    A dictionary of the given input files. Returns None if no input files were
    given.
    """

    # Gather the inputs
    inputs = [
            opt['--assemblies'], opt['--fwd_reads'], opt['--rev_reads']]

    # Extract input files
    input_dirs = [Path(d).resolve() for d in inputs]

    # The expected filetypes
    extensions = [
            biofile.Fasta.accepted_extensions,
            biofile.Fastq.accepted_extensions,
            biofile.Fastq.accepted_extensions]
    # The expected substrings
    substrings = [None, '1.', '2.']

    # Get the paths to all input files
    file_paths: List[List[Path]] = []
    for path, typ, substring in zip(input_dirs, extensions, substrings):
        file_paths.append(fmpaths.find(path, typ, substring))

    assemblies = file_paths[0]
    fwd_reads = file_paths[1]
    rev_reads = file_paths[2]
    _validate_args(assemblies, fwd_reads, rev_reads)

    return list(zip(assemblies, fwd_reads, rev_reads))


def _define_pipeline(
        input_files: List[Tuple[Path, Path, Path]],
        name: str,
        )-> ruffus.Pipeline:
    """Define the pipeline

    Parameters
    ----------
    input_files
        Input files produced by `_extract_input_files`
    name
        The name of the pipeline

    Returns
    -------
    The pipeline definition
    """
    path_strings = [fmpaths.as_str(paths) for paths in input_files]
    # Define the pipeline
    pipe = ruffus.Pipeline(name=name)

    # Add symlink step
    pipe.transform(
            task_func=fmruffus.PairedBowtie2Align,
            name='alignment',
            filter=ruffus.suffix('fasta'),
            input=path_strings,
            output='bam',
            output_dir='PairedBowtie2align')
    return pipe


def _main(opt: Dict[str, Any], name: str = 'decontam') -> None:
    """Main function

    Parse the command line arguments and decide how to proceed. Either run the
    pipeline, print the pipeline, and/or produce a flowchart.

    Parameters
    ----------
    opt
        A dictionary of parsed input parameters formatted as in `docopt`
    name: Optional
        The name of the pipeline
    """
    # Unpack the arguments
    if opt['--threads']:
        nthreads = int(opt['--threads'])

    input_files = _extract_input_files(opt)

    pipe = _define_pipeline(input_files, name=name)

    if opt['--flowchart']:
        pipe.printout_graph(opt['--flowchart'])

    if opt['--just-print']:
        pipe.printout()
    elif input_files:
        pipe.run(
            multithread=nthreads,
            logger=fmruffus.ROOT_LOGGER)

if __name__ == '__main__':
    opt = docopt(__doc__)
    _main(opt)
