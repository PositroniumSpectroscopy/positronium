#! python
'''
The Bohr model of positronium
'''
from __future__ import print_function, division
import numpy as np
from .constants import au_energy, au_distance

def energy(n1=1, n2=float('inf'), **kwargs):
    '''
        The energy interval between states n1 and n2.

        Parameters
        ----------
        n1 : int
            principal quantum number of state 1
        n2 : int
            principal quantum number of state 2
        **kwargs:
            unit="eV"
                "J", "eV", "meV", "ueV", "au", "Hartree",     [energy]
                "Hz", "kHz", "MHz", "GHz", "THz",             [frequency]
                "m", "cm", "mm", "um", "nm", "A", "pm", "fm", [vacuum wavelength]
                "m^-1", "cm^-1".                              [wavenumber]

        Returns
        -------
        float64
    '''
    unit = kwargs.get('unit', 'eV')
    interval = 0.25 * (1.0 / n1**2.0 - 1.0 / n2**2.0)
    if unit in ['au', 'Hartree']:
        return interval
    elif unit in au_energy:
        return au_energy[unit](interval)
    else:
        raise KeyError(unit + ' is not recognised as a suitable unit. See' + \
                              ' docstring for unit list.')

def radius(n=1, **kwargs):
    '''
        The Bohr radius.

        Parameters
        ----------
        n : int
            principal quantum number of the state
        **kwargs
            unit="m":
                "au", "Bohr"
                "m", "cm", "mm", "um", "nm", "A", "pm", "fm",

        Returns
        -------
        float64
    '''
    unit = kwargs.get('unit', 'm')
    rad = 2.0 * n**2.0
    if (unit == 'au') or (unit == 'Bohr'):
        return rad
    elif unit in au_distance:
        return au_distance[unit](rad)
    else:
        raise KeyError(unit + ' is not recognised as a suitable unit. See' + \
                              ' docstring for unit list.')
