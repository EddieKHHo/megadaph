"""Test the Megadaph Decontamination Pipeline using py.test"""
import uuid

from docopt import docopt
import pytest
from ruffus import formatter

import fmbiopy.fmlist as fmlist
from fmbiopy.fmtest import *
import megadaph_decontam_pipeline as pipe


@pytest.fixture(scope="module")
def doc():
    return pipe.__doc__


def test_docopt_parsing(doc):
    expected_keys = [
            '--threads', '--fwd_reads', '--rev_reads']
    values = '10', 'bar', 'foobar'
    argv = fmlist.interleave(expected_keys, values)
    parsed = docopt(doc, argv)
    for opt, val in zip(expected_keys, values):
        assert parsed[opt] == val


@pytest.fixture()
def opts(small):
    subdirs = ['fwd_reads', 'rev_reads']
    paths = [str((small / d).resolve()) for d in subdirs]
    return {
            '--threads' : 1,
            '--fwd_reads' : paths[0],
            '--rev_reads' : paths[1],
            '--flowchart' : False,
            '--just-print' : False,
            '--help' : False}


@pytest.fixture()
def print_only_opts(opts):
    opts['--just-print'] = True
    return opts


@pytest.fixture()
def flowchart_only_opts(print_only_opts, sandbox):
    print_only_opts['--flowchart'] = str(sandbox / 'flow.svg')
    return print_only_opts


@pytest.fixture()
def input_file_dict(dat):
    return {
            'fwd_reads' : sorted(dat['small']['fwd_reads']),
            'rev_reads' : sorted(dat['small']['rev_reads'])
    }


def test_check_input_files(opts, input_file_dict):
    assert pipe._extract_input_files(opts) == input_file_dict


def test_flowchart_only(flowchart_only_opts, randstr):
    pipe._main(flowchart_only_opts, name=randstr())
    assert Path(flowchart_only_opts['--flowchart']).exists()


def test_print_only(print_only_opts, randstr):
    pipe._main(print_only_opts, name=randstr())


def test_logging():
    pass

def test_run_pipeline(opts):
    pipe._main(opts)
