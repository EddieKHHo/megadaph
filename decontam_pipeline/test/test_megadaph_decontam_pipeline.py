"""Test the Megadaph Decontamination Pipeline using py.test"""
import os

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
            '<nthreads>' : 16,
            '--assemblies' : paths[0],
            '--fwd_reads' : paths[1],
            '--rev_reads' : paths[2]}

@pytest.fixture()
def zipped_dat_inputs(dat):
    return list(zip(
        dat['small']['assemblies'],
        dat['small']['fwd_reads'],
        dat['small']['rev_reads']))

def test_check_input_files(opts, zipped_dat_inputs):
    assert pipe._extract_input_files(opts) == zipped_dat_inputs
