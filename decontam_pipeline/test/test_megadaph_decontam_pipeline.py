"""Test the Megadaph Decontamination Pipeline using py.test"""
from uuid import uuid4

from docopt import docopt
import pytest

import fmbiopy.fmlist as fmlist
from fmbiopy.fmtest import *
import megadaph_decontam_pipeline as pipe


@pytest.fixture(scope="module")
def doc():
    return pipe.__doc__


def test_docopt_parsing(doc):
    expected_keys = [
            '--threads', '--assemblies', '--fwd_reads', '--rev_reads']
    values = '10', 'foo', 'bar', 'foobar'
    argv = fmlist.interleave(expected_keys, values)
    parsed = docopt(doc, argv)
    for opt, val in zip(expected_keys, values):
        assert parsed[opt] == val


@pytest.fixture()
def opts():
    # TODO: refactor to use dat keys
    location = Path('sandbox/small')
    subdirs = ['assemblies', 'fwd_reads', 'rev_reads']
    paths = [str((location / d).resolve()) for d in subdirs]
    return {
            '--threads' : 1,
            '--assemblies' : paths[0],
            '--fwd_reads' : paths[1],
            '--rev_reads' : paths[2],
            '--flowchart' : False,
            '--just-print' : False,
            '--help' : False}


@pytest.fixture()
def print_only_opts(opts):
    opts['--just-print'] = True
    return opts


@pytest.fixture()
def flowchart_only_opts(print_only_opts):
    print_only_opts['--flowchart'] = 'sandbox/flow.svg'
    return print_only_opts


@pytest.fixture()
def zipped_dat_inputs(dat):
    return list(zip(
        sorted(dat['small']['assemblies']),
        sorted(dat['small']['fwd_reads']),
        sorted(dat['small']['rev_reads'])))


@pytest.fixture()
def uid():
    return uuid4().hex


def test_check_input_files(opts, zipped_dat_inputs):
    assert pipe._extract_input_files(opts) == zipped_dat_inputs


def test_flowchart_only(flowchart_only_opts, uid):
    pipe._main(flowchart_only_opts, name=uid)
    assert Path(flowchart_only_opts['--flowchart']).exists()


def test_print_only(print_only_opts, uid):
    pytest.set_trace()
    pipe._main(print_only_opts, name=uid)


def test_logging():
    pass

# def test_run_pipeline(opts):
#     _main(opts)
