#! python
'''
Atomic state of positronium

Energy levels are calculated using equations from

> Richard A. Ferrell (1951) Phys. Rev. 84, 858
> http://dx.doi.org/10.1103/PhysRev.84.858

(does not include the Lamb shift/ radiative corrections).

'''
from __future__ import print_function, division
import numpy as np
import attr
from .constants import au_energy, alpha, mu_me

@attr.s()
class Ps(object):
    """ attrs class to represent the quantum state |n l S J [MJ]>.
    """
    n = attr.ib(convert=int)
    @n.validator
    def check_n(self, attribute, value):
        if not value > 0:
            raise ValueError("n must be a positive integer.")
    l = attr.ib(convert=int)
    @l.validator
    def check_l(self, attribute, value):
        if not value < self.n:
            raise ValueError("l must be an integer less than n.")
    S = attr.ib(convert=int)
    @S.validator
    def check_S(self, attribute, value):
        if value not in [0, 1]:
            raise ValueError("S must be 0 or 1.")
    J = attr.ib(convert=int)
    @J.validator
    def check_J(self, attribute, value):
        if self.l == 0 and not value == self.S:
            raise ValueError("If l == 0, J must be equal to S.")
        elif self.S == 0 and not value == self.l:
            raise ValueError("If S == 0, J must be equal to l.")
        elif (not self.l - self.S <= value <= self.l + self.S):
            raise ValueError("J must be in range l - S <= J <= l + S.")
    # MJ is optional (levels are degenerate in zero field)
    MJ = attr.ib(default=None)
    @MJ.validator
    def check_MJ(self, attribute, value):
        if value is not None and (not -self.J <= value <= self.J):
            raise ValueError("MJ must be in the range of -J <= MJ <= J.")

    def __attrs_post_init__(self):
        """ check if MJ is used and calculate the zero-field energy of the state.
        """
        if self.MJ is not None:
            self.MJ = int(self.MJ)
        self.e0 = energy(self.n, self.l, self.S, self.J)
        ## no quantum defect
        self.n_eff = float(self.n) #np.sqrt(- mu_me * 0.5 / self.e0)

    def __str__(self):
        """ print quantum numbers like |n l S J [MJ] >
        """
        if self.MJ is None:
            return u"\u007C {} {} {} {} \u27E9".format(self.n, self.l, self.S, self.J)
        else:
            return u"\u007C {} {} {} {} {} \u27E9".format(self.n, self.l, self.S, self.J, self.MJ)

    def asdict(self):
        """ quantum numbers as a dictionary.
        """
        return attr.asdict(self)

    def tex(self, show_MJ=True):
        """ Tex string of the form n^{2S + 1}L_{J} (M_J = {MJ})
        """
        if self.MJ is None:
            show_MJ = False
        L = 'SPDFGHIKLMNOQRTUVWXYZ'[int(self.l%22)]
        tex_str = r'$%d^{%d}'%(self.n, 2*self.S + 1) + L + r'_{%d}'%(self.J)
        if show_MJ:
            tex_str = tex_str + r'\,' + r'(M_J = %d)$'%self.MJ
        else:
            tex_str = tex_str + r'$'
        return tex_str
    
    def energy(self, unit='eV'):
        """
        """
        if unit in ['au', 'Hartree']:
            return self.e0
        elif unit in au_energy:
            return au_energy[unit](self.e0)
        else:
            raise KeyError(unit + ' is not recognised as a suitable unit.')

def epsilon(l, S, J):
    """ scaling of the fine structure shift.
    """
    if S == 0:
        # singlet
        epsilon = 0.0
    elif S == 1:
        # triplet
        delta = int(l == 0)
        if J == l + 1:
            omega = (3*l + 4)/ ((l + 1) * (2*l + 3))
        elif J == l:
            omega =  -1.0 / (l*(l + 1))
        elif J == l - 1:
            omega = - (3*l - 1.0)/ (l*(2*l - 1))
        else:
            raise ValueError("The total angular momentum quantum number 'J' must " + \
                             "be in the range l - 1 < J < l + 1")
        epsilon = 7.0 / 6.0 * delta + (1 - delta) / (2.0 * (2 * l + 1)) * omega
    else:
        raise ValueError("The total spin quantum number 'S' must be 0 or 1.")
    return epsilon

def energy_fs(n, l, S, J):
    """ first-order fine structure shift for state |n l S J >
    
        > H. A. Bethe and E. E. Salpeter (1957)
        > Quantum Mechanics of One- and Two-Electron Systems
    """
    # TODO - special case for l=0 and l=1.
    en = (11.0/ (32 * n**4) + (epsilon(l, S, J) - 1.0 / (2*l + 1)) * 1.0/ (n**3.0)) * alpha**2.0
    return mu_me * en

def energy(n, l, S, J):
    """ energy levels, including fine structure.
    """
    en = - mu_me * 0.5 / (n**2.0) + energy_fs(n, l, S, J)
    return en
