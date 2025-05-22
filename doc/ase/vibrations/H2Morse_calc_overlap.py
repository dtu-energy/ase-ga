from ase.calculators.h2morse import H2Morse, H2MorseExcitedStatesCalculator
from ase.vibrations.resonant_raman import ResonantRamanCalculator

atoms = H2Morse()
rmc = ResonantRamanCalculator(
    atoms, H2MorseExcitedStatesCalculator, overlap=lambda x, y: x.overlap(y)
)
rmc.run()
