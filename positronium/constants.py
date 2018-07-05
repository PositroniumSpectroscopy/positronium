#! python
'''
Physical constants

'''
from __future__ import print_function, division
from math import floor, log10
import webbrowser

metric_prefix = dict([(18, 'E'), (15, 'P'), (12, 'T'), (9, 'G'),(6, 'M'), (3, 'k'), (0, ''),
                      (-3, 'm'), (-6, 'u'), (-9, 'n'), (-12, 'p'), (-15, 'f'), (-18, 'a')])

#CODATA 2014, DOI: 10.1103/RevModPhys.88.035009
c = 299792458.0 ## speed of light in vacuum
h = 6.626070040e-34
h_uncertainty = 8.1e-42
hbar = 1.054571800e-34
Ry = 10973731.568508
e = 1.6021766208e-19
m_e = 9.10938356e-31
alpha = 7.2973525664e-3
m_u = 1.660539040e-27
En_h = alpha**2.0 * m_e * c**2.0
a_0 = hbar/ (m_e * c * alpha)
mu_B = e * hbar / (2.0 * m_e)

# Ps
mass_Ps = m_Ps  = 2.0 * m_e
reduced_mass_Ps = mu_Ps = m_e/ 2.0
mu_me = 0.5
Rydberg_Ps = Ry_Ps = Ry / 2.0
a_Ps = 2.0 * a_0

class MeasuredValue(float):
    """ subclass of built-in type float
    
        Parameters
        ----------
        value : float
        uncertainty : float
        unit : str

        Attributes
        ----------
        source : str
            source reference
        url : str
            url to source

        Methods
        -------
        article()
            opens source url in a web browser
    """
    def __new__(cls, value, uncertainty, unit):
        measured = float.__new__(cls, value)
        measured.value = value
        measured.uncertainty = uncertainty
        measured.unit = unit
        measured.source = None
        measured.url = None
        return measured

    def __str__(self):
        ''' format for print() by rounding the uncertainty to 2 significant figures and 
            rescaling value by nearest metric prefix.
        '''
        exp = floor(log10(abs(self.value)))
        si = 3 * floor(exp / 3.0)
        if si not in metric_prefix:
            si = 0
        scale = 10**(-si)
        dp = -floor(log10(abs(self.uncertainty * scale))) + 1
        uncertainty_str = format(self.uncertainty * scale, '.{}f'.format(dp))
        value_str = format(self.value * scale, '.%df'%dp)
        return u"{} \u00B1 {} {}{}".format(value_str, uncertainty_str, metric_prefix[si], self.unit)

    def article(self):
        ''' open url to constant reference'''
        if self.url is not None:
            webbrowser.open(self.url)
        else:
            raise LookupError('The source url has not been specified.')

# ground-state decay rate / lifetime
## para-positronium
decay_pPs = MeasuredValue(7990.9e6, 1.7e6, 'Hz')
lifetime_pPs = MeasuredValue(1.0 / decay_pPs, decay_pPs.uncertainty/ decay_pPs**2.0, 's')
for con in [decay_pPs, lifetime_pPs]:
    setattr(con, 'source', 'Al-Ramadhan, A. H. and Gidley, D. W. (1994) ' + \
                           'Phys. Rev. Lett. 72, 1632 ')
    setattr(con, 'url', 'http://dx.doi.org/10.1103/PhysRevLett.72.1632')

## ortho-positronium
decay_oPs = MeasuredValue(7.0404e6, 0.0018e6, 'Hz')
lifetime_oPs = MeasuredValue(1.0 / decay_oPs, decay_oPs.uncertainty/ decay_oPs**2.0, 's')
for con in [decay_oPs, lifetime_oPs]:
    setattr(con, 'source', 'R. S. Vallery, P. W. Zitzewitz, and D. W. Gidley (2003) ' + \
                           'Phys. Rev. Lett. 90, 203402')
    setattr(con, 'url', 'http://dx.doi.org/10.1103/PhysRevLett.90.203402')

# ground-state hyperfine splitting
frequency_hfs = MeasuredValue(203.38910e9, 000.00074e9, 'Hz')
energy_hfs = MeasuredValue(h * frequency_hfs,
                           ((h * frequency_hfs.uncertainty)**2.0 + (h_uncertainty * frequency_hfs)**2.0)**0.5,
                           'J')
for con in [frequency_hfs, energy_hfs]:
    setattr(con, 'source', 'M. W. Ritter, P. O. Egan, V. W. Hughes, and ' +
                           'K. A. Woodle (1984) Phys. Rev. A 30, 1331')
    setattr(con, 'url', 'http://dx.doi.org/10.1103/PhysRevA.30.1331')

# 1^3S-2^3S interval
frequency_1s2s = MeasuredValue(1233607216.4e6, 3.2e6, 'Hz')
energy_1s2s = MeasuredValue(h * frequency_1s2s,
                            ((h * frequency_1s2s.uncertainty)**2.0 + (h_uncertainty * frequency_1s2s)**2.0)**0.5,
                            'J')
for con in [frequency_1s2s, energy_1s2s]:
    setattr(con, 'source', 'Fee, M.S. et al. (1993) Phys. Rev. Lett. 70, 1397')
    setattr(con, 'url', 'http://dx.doi.org/10.1103/PhysRevLett.70.1397')

# rescale atomic units
## energy / energy interval / wavelength/ wavenumbers
au_energy = dict([('J', lambda x: x * En_h),
                  ('eV', lambda x: x * En_h / e),
                  ('meV', lambda x: 1e3 * En_h / e),
                  ('ueV', lambda x: 1e6 * En_h / e),
                  ('Hz', lambda x: x * 2.0 * Ry * c),
                  ('kHz', lambda x: 1e-3 * x * 2.0 * Ry * c),
                  ('MHz', lambda x: 1e-6 * x * 2.0 * Ry * c),
                  ('GHz', lambda x: 1e-9 * x * 2.0 * Ry * c),
                  ('THz', lambda x: 1e-12 * x * 2.0 * Ry * c),
                  ('m', lambda x: 1.0 / (x * 2.0 * Ry)),
                  ('cm', lambda x: 1e2 / (x * 2.0 * Ry)),
                  ('mm', lambda x: 1e3 / (x * 2.0 * Ry)),
                  ('um', lambda x: 1e6 / (x * 2.0 * Ry)),
                  ('nm', lambda x: 1e9 / (x * 2.0 * Ry)),
                  ('A', lambda x: 1e10 / (x * 2.0 * Ry)),
                  ('pm', lambda x: 1e12 / (x * 2.0 * Ry)),
                  ('fm', lambda x: 1e15 / (x * 2.0 * Ry)),
                  ('m^-1', lambda x: x * 2.0 * Ry),
                  ('cm^-1', lambda x: 1e-2 * x * 2.0 * Ry)])
## distance
au_distance = dict([('m', lambda x: x * a_0),
                    ('cm', lambda x: x * a_0 * 1e2),
                    ('mm', lambda x: x * a_0 * 1e3),
                    ('um', lambda x: x * a_0 * 1e6),
                    ('nm', lambda x: x * a_0 * 1e9),
                    ('A', lambda x: x * a_0 * 1e10),
                    ('pm', lambda x: x * a_0 * 1e12),
                    ('fm', lambda x: x * a_0 * 1e15)])
