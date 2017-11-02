#!/usr/bin/env python
"""Megadaph Decontamination Pipeline

Built using Ruffus pipelining infrastructure (See http://www.ruffus.org.uk/)

Pipeline steps are as follows:
    Index Fasta Files -> Align Fastq files ->

Usage:
    megadaph_decontam_pipeline.py [--flowchart=<file.svg>] [--verbose]
        [--just-print] (--threads=<n>)
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

Notes:
    Matched Fastq files must share the same prefix and be named using dot
    filename convention (E.g FA_SC.1.trim.fastq). Forward reads must have '.1'
    immediately after shared prefix, reverse reads must have '.2'.
"""

from pathlib import Path
from typing import (
        Any,
        Dict,
        Iterable,
        List,
        )

from docopt import docopt
from ruffus import (
        formatter,
        Pipeline,
        ruffus_utility,
        suffix,
        )
from fmbiopy.biofile import (
        BiofileGroup,
        Fastq,
        MatchedPrefixGroup,
        )
from fmbiopy.fmpaths import (
        as_strs,
        find,
        )
from fmbiopy.fmruffus import (
        apply_,
        Centrifuge,
        format_input,
        output_format,
        PairedBowtie2Align,
        SymlinkInputs,
        )

"""The extra parameters for each task in the pipeline"""
PARAM = {
        'PairedBowtie2Align' : ['-X', '1000'],
        }

"""The working directory"""
WD = Path.cwd()

"""The output directory"""
OUTPUT_DIR = WD / 'pipe'

def _get_gzip_status(files: Iterable[Path]) -> bool:
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
    """

    inputs = [fwd_reads, rev_reads]
    gzipped = [_get_gzip_status(inp) for inp in inputs]
    MatchedPrefixGroup([
        BiofileGroup(fwd_reads, filetype='fastq', gzipped=gzipped[0]),
        BiofileGroup(rev_reads, filetype='fastq', gzipped=gzipped[1])
        ]).validate()


def _extract_input_files(
        opt: Dict[str, Any]
        )-> Dict[str, List[Path]]:
    """Extract the input files from the given directories

    Input file prefixes are checked to ensure they can be matched correctly and
    have the correct file extensions.

    Parameters
    ----------
    opt
        A dictionary of parsed input parameters formatted as in `docopt`

    Returns
    -------
    A dictionary of the input files
    """

    # Gather the inputs
    inputs = [opt['--fwd_reads'], opt['--rev_reads']]

    # Extract input files
    input_dirs = [Path(d).resolve() for d in inputs]

    # The expected filetypes
    extensions = [Fastq.accepted_extensions] * 2
    # The expected substrings
    substrings = ['1.', '2.']

    # Get the paths to all input files
    file_paths: List[List[Path]] = []
    for path, typ, substring in zip(input_dirs, extensions, substrings):
        file_paths.append(find(path, typ, substring))

    fwd_reads = file_paths[0]
    rev_reads = file_paths[1]
    _validate_args(fwd_reads, rev_reads)

    output_dict = {
            'fwd_reads' : fwd_reads,
            'rev_reads' : rev_reads}

    return output_dict


def _define_pipeline(
        input_files: Dict[str, List[Path]],
        name: str,
        )-> Pipeline:
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
    path_strings = [as_strs(paths) for paths in input_files.values()]
    fwd_reads = path_strings[0]
    rev_reads = path_strings[1]
    gzipped = _get_gzip_status(input_files['fwd_reads'])

    link_suffixes = (['.R1.fastq', '.R2.fastq'])
    if gzipped:
        link_suffixes = [''.join([suff, '.gz']) for suff in link_suffixes]

    # Define the pipeline
    pipe = Pipeline(name=name)

    # Symlink Inputs
    fwd_read_links = pipe.transform(
            task_func=SymlinkInputs().run,
            filter=format_input(Fastq.accepted_extensions),
            input=fwd_reads,
            output=output_format([link_suffixes[0]], dirname='symlink_inputs'))
    rev_read_links = pipe.transform(
            task_func=SymlinkInputs().run,
            filter=format_input(Fastq.accepted_extensions),
            input=rev_reads,
            output=output_format([link_suffixes[1]], dirname='symlink_inputs'))
    # Run Fastqc
    # Trim Reads
    # Rerun Fastqc
    # Summarize Trimming
    # Assemble Reads
    # Build Centrifuge Database
    # Simulate Minimum hit length
    # Run Centrifuge
    # cent = pipe.transform(
    #         task_func=Centrifuge().run,
    #         filter=suffix('fasta'),
    #         output='tsv',
    #         output_dir='centrifuge')
    # align = pipe.transform(
    #         task_func=PairedBowtie2Align(PARAM['PairedBowtie2Align']).run,
    #         name='Alignment',
    #         filter=suffix('fasta'),
    #         input=[assembly_links, fwd_read_links, rev_read_links],
    #         output='bam',
    #         output_dir='align_to_assemblies')
    return pipe


def _main(opt: Dict[str, Any], name: str = 'decontam') -> None:
    """Main function

    Parse the command line arguments and decide how to proceed. Either run the
    pipeline, print the pipeline, and/or produce a flowchart.

    Parameters
    ----------
    opt
        A dictionary of parsed input parameters formatted as in `docopt`
    name: optional
        The name of the pipeline
    """
    # Unpack arguments
    if opt['--threads']:
        nthreads = int(opt['--threads'])
    input_files = _extract_input_files(opt)

    # Define the pipeline
    pipe = _define_pipeline(input_files, name=name)

    # Run or print pipeline
    if opt['--flowchart']:
        pipe.printout_graph(opt['--flowchart'])
    if opt['--just-print']:
        pipe.printout()
    elif input_files:
        pipe.run(multithread=nthreads, exceptions_terminate_immediately=True)


if __name__ == '__main__':
    opt = docopt(__doc__)
    _main(opt)
