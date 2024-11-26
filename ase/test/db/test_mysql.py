import os

import pytest

from ase import Atoms
from ase.build import molecule
from ase.calculators.emt import EMT
from ase.db import connect


@pytest.fixture(scope='module')
def url(mysql_port):
    pytest.importorskip('pymysql')

    on_ci_server = 'CI_PROJECT_DIR' in os.environ

    if on_ci_server:
        # CI server configured to use non-standard port 3307
        # instead of the default 3306 port. This is to test
        # for proper passing of the port for creating the
        # mysql connection
        db_url = f'mysql://root:ase@mysql:{mysql_port}/testase_mysql'
        # HOST = 'mysql'
        # USER = 'root'
        # PASSWD = 'ase'
        # DB_NAME = 'testase_mysql'
    else:
        db_url = os.environ.get('MYSQL_DB_URL')
        # HOST = os.environ.get('MYSQL_HOST', None)
        # USER = os.environ.get('MYSQL_USER', None)
        # PASSWD = os.environ.get('MYSQL_PASSWD', None)
        # DB_NAME = os.environ.get('MYSQL_DB_NAME', None)

    if db_url is None:
        msg = ('Not on GitLab CI server. To run this test '
               'host, username, password and database name '
               'must be in the environment variables '
               'MYSQL_HOST, MYSQL_USER, MYSQL_PASSWD and '
               'MYSQL_DB_NAME, respectively.')
        pytest.skip(msg)
    return db_url


@pytest.fixture()
def db(url):
    return connect(url)


@pytest.fixture()
def h2o():
    return molecule('H2O')


def test_connect(db):
    db.delete([row.id for row in db.select()])


def test_write_read(db):
    co = Atoms('CO', positions=[(0, 0, 0), (0, 0, 1.1)])
    uid = db.write(co, tag=1, type='molecule')

    co_db = db.get(id=uid)
    atoms_db = co_db.toatoms()

    assert len(atoms_db) == 2
    assert atoms_db[0].symbol == co[0].symbol
    assert atoms_db[1].symbol == co[1].symbol
    assert co_db.tag == 1
    assert co_db.type == 'molecule'


def test_write_read_with_calculator(db, h2o):
    calc = EMT(dummy_param=2.4)
    h2o.calc = calc

    uid = db.write(h2o)

    h2o_db = db.get(id=uid).toatoms()

    # Back in the days we allowed reconstructing calculators.
    # For security we don't anymore.
    assert h2o_db.calc is None

    # Check that get_atoms function works
    db.get_atoms(H=2)
    # XXX We should assert something should we not?


def test_update(db, h2o):
    h2o = molecule('H2O')

    uid = db.write(h2o, type='molecule')
    db.update(id=uid, type='oxide')

    atoms_type = db.get(id=uid).type

    assert atoms_type == 'oxide'


def test_delete(db, h2o):
    h2o = molecule('H2O')
    uid = db.write(h2o, type='molecule')

    # Make sure that we can get the value
    db.get(id=uid)
    db.delete([uid])

    with pytest.raises(KeyError):
        db.get(id=uid)


def test_read_write_bool_key_value_pair(db, h2o):
    # Make sure we can read and write boolean key value pairs
    uid = db.write(h2o, is_water=True, is_solid=False)
    row = db.get(id=uid)
    assert row.is_water
    assert not row.is_solid
