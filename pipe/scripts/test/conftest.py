from uuid import uuid4

from plumbum import local
from pytest import fixture


@fixture(scope='session')
def tempdir(tmpdir_factory):
    return local.path(str(tmpdir_factory.mktemp(uuid4().hex)))
