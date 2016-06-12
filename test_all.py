#! python
'''
pytest
'''
from positronium import constants 
from positronium import Bohr
from positronium import Ferrell
from positronium import Ps


def test_answer():
    assert Bohr.energy(1, 2, unit='au') == 0.1875
    assert round((Ferrell.energy_level(1, 0, 0, 0, unit='J') -\
                  Ferrell.energy_level(1, 0, 1, 1, unit='J'))/constants.h) \
                  == -204386630477.0
    x1 = Ps(n=1, l=0, S=1, J=1)
    assert round(1e9*x1.annihilation(unit='s'), 3) == 142.037