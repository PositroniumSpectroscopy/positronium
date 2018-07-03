#! python
"""
    Wavefunction for positronium in state |n, l, m>  
"""
from __future__ import print_function, division
from math import factorial
from scipy.special import sph_harm, hyp1f1
import numpy as np

def angular_wf(l, m):
    """ 
        The angular part of the wavefunction.
        
        Parameters
        ----------
        l : int
            orbital angular momentum quantum number
        m : int
            magnetic quantum number
        
        Returns
        -------
        Y_lm(theta, phi)
    """
    return lambda theta, phi: sph_harm(m, l, phi, theta) * (-1.0)**m

def radial_wf(n, l):
    """ 
        The radial part of the wavefunction, see
 
        > Quantum Mechanics of one and two electron atoms,
        > H. A. Bethe and E. E. Salpeter 1957
        > (page 15)
 
        Parameters
        ----------
        n : int
            principal quantum number
        l : int
            orbital angular momentum quantum number
        
        Returns
        -------
        R_nl(r)

        Notes
        -----
        r in units of the Bohr radius, a_Ps.
    """
    # hydrogen wf
    Z = 1.0
    epsilon = 2.0 * Z / n
    c1 = np.sqrt(factorial(n + l) / (2.0 * n * factorial(n - l - 1))) * \
         1.0 / (factorial(2.0 * l + 1.0))  
    return lambda r: c1 * epsilon**(3.0 / 2.0) * np.exp(-0.5 * epsilon * r) * (epsilon * r)**l *\
                     hyp1f1(-(n - l - 1), 2 * l + 2, (epsilon * r))

def wf(n, l, m):
    """ 
        Solution to the Schrodinger equation for a 1/r potential.

        > Quantum Mechanics of one and two electron atoms,
        > H. A. Bethe and E. E. Salpeter 1957

        Parameters
        ----------
        n : int
            principal quantum number
        l : int
            orbital angular momentum quantum number
        m : int
            magnetic quantum number
            
        Returns
        -------
        u_nlm(r, theta, phi) = R_nl(r) * Y_lm(theta, phi)

        Notes
        -----
        r in units of the Bohr radius, a_Ps.
    """
    rad = radial_wf(n, l)
    ang = angular_wf(l, m)
    return lambda r, theta, phi: rad(r) * ang(theta, phi)