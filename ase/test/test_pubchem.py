from ase.data.pubchem import (
    pubchem_atoms_conformer_search,
    pubchem_atoms_search,
    pubchem_conformer_search,
    pubchem_search,
)


def test_pubchem():

    # check class functionality
    data = pubchem_search('ammonia', mock_test=True)
    data.get_atoms()
    data.get_pubchem_data()
    # XXX maybe verify some of this data?

    # check the various entry styles and the functions that return atoms
    pubchem_search(cid=241, mock_test=True).get_atoms()
    pubchem_atoms_search(smiles='CCOH', mock_test=True)
    pubchem_atoms_conformer_search('octane', mock_test=True)
    # (maybe test something about some of the returned atoms)

    # check conformer searching
    confs = pubchem_conformer_search('octane', mock_test=True)
    for _ in confs:
        pass
    try:  # check that you can't pass in two args
        pubchem_search(name='octane', cid=222, mock_test=True)
        raise Exception('Test Failed')
    except ValueError:
        pass

    try:  # check that you must pass at least one arg
        pubchem_search(mock_test=True)
        raise Exception('Test Failed')
    except ValueError:
        pass
