#! python
""" 
    (approx) transition intervals between two states.
""" 
from .constants import En_h, h, c

def interval(state_1, state_2):
    "The energy interval between state_1 and state_2 in J"
    return abs(state_1.e0 - state_2.e0) * En_h

def frequency(state_1, state_2):
    """ The frequency interval between state_1 and state_2 in GHz.
    """
    return 1e-9 * interval(state_1, state_2) / h

def wavelength(state_1, state_2):
    """ The photon wavelength corresponding to the energy interval between 
        state_1 and state_2 in nm.
    """
    return 1e9 * h * c / interval(state_1, state_2)

def wavenumbers(state_1, state_2):
    """ The energy interval between state_1 and state_2 in wavemeters (cm^-1).
    """
    return 1e-2 * interval(state_1, state_2) / (h * c)