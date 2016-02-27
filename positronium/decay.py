from __future__ import print_function, division
#! python
"""
Created on Sat Feb 27 11:35:54 2016

@author: Adam
"""
import numpy as np
import positronium.constants as constants
t_0 = 3*constants.hbar / (2.0 * constants.alpha**5 * constants.reduced_mass_Ps * constants.c**2)

def radiative(n, l=0, unit='s'):
    '''
    A universal formula for the radiative mean lifetime of hydrogenlike states (n,l).
    The formula is accurate to at least 6% for the lowest states and to a much higher
    degree of accuracy for highly excited states.
    
        Semiclassical estimation of the radiative mean lifetimes of hydrogenlike states
        Hermann Marxer and Larry Spruch
        Phys. Rev. A 43, 1268
        https://dx.doi.org/10.1103/PhysRevA.43.1268
    
    kwargs:
        unit:
            s, ms, us, ns, ps,                      [lifetime]
            Hz, kHz, MHz, GHz, THz, PHz, EHz,       [rate]

    defaults:
        n = 1
        l = 0
        unit='s'       
    '''
    rescale = {'s': (lambda x: x),
               'ms': (lambda x: x*1e3),
               'us': (lambda x: x*1e6),
               'ns': (lambda x: x*1e9),
               'ps': (lambda x: x*1e12),
               'Hz': (lambda x: 1.0/ x),
               'kHz': (lambda x: 1.0/ x * 1e-3),
               'MHz': (lambda x: 1.0/ x * 1e-6),
               'GHz': (lambda x: 1.0/ x * 1e-9),
               'THz': (lambda x: 1.0/ x * 1e-12),
               'PHz': (lambda x: 1.0/ x * 1e-15),
               'EHz': (lambda x: 1.0/ x * 1e-18),
               }
    if unit not in rescale:
        raise KeyError('"' + unit + '" is not recognised as a suitable unit. See' +
                           ' docstring for unit list.')
    # make it possible to check values from arrays
    n = np.array([n]).flatten()
    l = np.array([l]).flatten()
    if not (l >= 0).any():
        raise ValueError("'l' must be greater than or equal to zero")
    if not np.greater(n, l).all():
        raise ValueError("'n' must be greater than 'l'")
    else:
        lifetime = t_0*np.power(n,3)*np.multiply(l, np.add(l, 1))
        try:
            result = rescale[unit](lifetime)
        except ZeroDivisionError:
            result = float('inf')
        except:
            raise
        finally:
            if len(result) == 1:
                return result[0]
            else:
                return result