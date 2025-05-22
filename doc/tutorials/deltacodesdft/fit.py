import json
from pathlib import Path
from typing import Tuple

from ase.eos import EquationOfState as EOS
from ase.io import read


def fit(symbol: str) -> Tuple[float, float, float, float]:
    V = []
    E = []
    for atoms in read(f'{symbol}.traj@:'):
        V.append(atoms.get_volume() / len(atoms))
        E.append(atoms.get_potential_energy() / len(atoms))
    eos = EOS(V, E, 'birchmurnaghan')
    eos.fit(warn=False)
    e0, B, Bp, v0 = eos.eos_parameters
    return e0, v0, B, Bp


data = {}  # Dict[str, Dict[str, float]]
for path in Path().glob('*.traj'):
    symbol = path.stem
    e0, v0, B, Bp = fit(symbol)
    data[symbol] = {
        'emt_energy': e0,
        'emt_volume': v0,
        'emt_B': B,
        'emt_Bp': Bp,
    }

Path('fit.json').write_text(json.dumps(data))
