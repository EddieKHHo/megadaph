"""Test the Megadaph Decontamination Pipeline using py.test"""
from docopt import docopt
import fmbiopy.fmstring as fmstring
import megadaph_decontam_pipeline as pipe
import pytest

@pytest.fixture(scope="module")
def doc():
    return pipe.__doc__

def test_docopt_parsing(doc):
    options = ['-p', '-a', '-r1', '-r2']
    expected_keys = ['<nthreads>', '--assemblies', 'FORWARD_READ_DIR',
            'REVERSE_READ_DIR']
    values = '10', 'foo', 'bar', 'foobar'
    argv = fmstring.interleave(options, values)
    parsed = docopt(doc, argv)
    for opt, val in zip(expected_keys, values):
        assert parsed[opt] == val
