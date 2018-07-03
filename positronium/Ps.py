#! python
'''
Atomic state of positronium
'''
from __future__ import print_function, division
import numpy as np
from . import constants, Ferrell

TAU_0 = 3.0*constants.hbar / (2.0 * constants.alpha**5.0 * \
        constants.reduced_mass_Ps * constants.c**2.0)

class Ps(object):
    """ Ps atomic state described by the quantum numbers n, l, m, S and J.
    """
    def __init__(self, n, l=0, m=0, S=0, J=None):
        self.n = n
        self.l = l
        self.m = m
        self.S = S
        if J is None:
            self.J = self.l + self.S
        else:
            self.J = J
        self.check_quantum_numbers()

    def check_quantum_numbers(self):
        """" Validate chosen quantum numbers. """
        if not isinstance(self.n, (int)) or self.n <= 0:
            raise ValueError("The principal quantum number 'n' must be a " + \
                             "positive integer")
        if not isinstance(self.l, (int)) or self.l < 0 or self.l >= self.n:
            raise ValueError("The orbital angular momentum quantum number 'l' " + \
                             "must be a positive or zero integer, and must be " + \
                             "less than the principal quantum number 'n'.")
        if not isinstance(self.S, (int)) or not (self.S == 0 or self.S == 1):
            raise ValueError("The total spin quantum number 'S' " + \
                             "must be an integer equal to 1 or 0.")
        if not isinstance(self.J, (int)) or not (self.l - self.S) <= self.J <= (self.l + self.S):
            raise ValueError("The total angular momentum quantum number 'J' " + \
                             "must be in the range l - 1 < J < l + 1")
        if self.l == 0 and self.J != self.S:
            raise ValueError("The total angular momentum quantum number 'J' " + \
                             "must be equal to S if l = 0.")
        return "ok"

    # ------
    # energy
    # ------

    def energy(self, unit='J'):
        """ Energy levels of positronium.  Includes fine structure. Does not
            include the Lamb shift/ radiative corrections.

            Richard A. Ferrell (1951) Phys. Rev. 84, 858
            http://dx.doi.org/10.1103/PhysRev.84.858

            quantum numbers:
                n, l, S, J
        """
        return Ferrell.energy_level(self.n, self.l, self.S, self.J, unit=unit)

    # ----------
    # ionisation
    # ----------

    def ionization_efield(self, unit='V m^-1'):
        """ Ionization in an electric field, for a given n.

            For each value of n the ionization field for the outermost Stark state with
            a negative Stark shift is approximately equal to the classical ionization
            field.

            T. F. Gallagher, Rydberg Atoms (Cambridge University Press, Cambridge,
            England, 1994).

            quantum numbers:
                n

            kwargs:
                unit:
                    V m^-1, V cm^-1      [electric field]

            defaults:
                n = 1
                unit='V cm^-1'
        """
        rescale = {'V m^-1': (lambda x: x),
                   'V cm^-1': (lambda x: 0.01 * x)}
        if unit not in rescale:
            raise KeyError('"' + unit + '" is not recognised as a suitable unit. See' +
                           ' docstring for unit list.')
        else:
            efield = np.power(self.n, -4.0) * 2.0 * constants.Ryd_Ps * \
                               constants.h * constants.c / (constants.e * constants.a_Ps * 9.0)
            return rescale[unit](efield)

    # ---------------
    # radiative decay
    # ---------------

    def radiative_decay(self, unit='Hz'):
        """ A universal formula for the radiative mean rate/ lifetime of hydrogenlike
            states. The formula is accurate to at least 6% for the lowest
            states and to a much higher degree of accuracy for highly excited states.

            Semiclassical estimation of the radiative mean lifetimes of hydrogenlike states
            Hermann Marxer and Larry Spruch
            Phys. Rev. A 43, 1268
            https://dx.doi.org/10.1103/PhysRevA.43.1268

            quantum numbers:
                n, l>0

            kwargs:
                unit:
                    s, ms, us, ns, ps,                      [lifetime]
                    Hz, kHz, MHz, GHz, THz, PHz, EHz,       [rate]

            defaults:
                unit='Hz'
        """
        if self.l == 0:
            raise ValueError('This semiclassical estimate does not work for l=0 states.')
        else:
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
                       'EHz': (lambda x: 1.0/ x * 1e-18)}
            if unit not in rescale:
                raise KeyError('"' + unit + '" is not recognised as a suitable unit. See' +
                               ' docstring for unit list.')
            else:
                lifetime = TAU_0 * np.multiply(np.power(self.n, 3.0),
                                               np.multiply(self.l, np.add(self.l, 1)))
                return rescale[unit](lifetime)

    # -----------------
    # self annihilation
    # -----------------

    def annihilation(self, unit='Hz'):
        """ Calculate the annihilation rate/ lifetime for S states (l=0).

            quantum numbers:
                n, l=0
        """
        if self.l != 0:
            raise ValueError('This calculation only works for l=0 states.')
        else:
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
                       'EHz': (lambda x: 1.0/ x * 1e-18)}
            if unit not in rescale:
                raise KeyError('"' + unit + '" is not recognised as a suitable unit.')
            else:
                if self.S == 0:
                    return constants.lifetime_pPs * self.n**3.0
                elif self.S == 1:
                    return constants.lifetime_oPs * self.n**3.0
                else:
                    raise ValueError("The total spin quantum number 'S' must be 0 or 1.")

    # -----
    # LaTeX
    # -----

    def tex(self):
        """ Convert the quantum numbers into a tex string of the form n^{2S + 1}L_{J}

            quantum numbers:
                n, l, S, J

            In IPython use Latex to render, e.g.

            >>> from IPython.display import Latex
            >>> from positronium import Ps
            >>> x = Ps(n=1, l=0, S=1, J=1)
            >>> Latex(x.tex())
        """
        L = 'SPDFGHIKLMNOQRTUVWXYZ'[int(self.l%22)]
        return r'$%d^{%d}'%(self.n, 2*self.S + 1) + L + r'_{%d}$'%(self.J)
