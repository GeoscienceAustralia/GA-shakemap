# -*- coding: utf-8 -*-
# vim: tabstop=4 shiftwidth=4 softtabstop=4
#
# Copyright (C) 2015-2019 GEM Foundation
#
# OpenQuake is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# OpenQuake is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with OpenQuake. If not, see <http://www.gnu.org/licenses/>.

"""
Module exports :class:'AtkinsonEtAl2014Cal',
                      'AtkinsonEtAl2014CEUS'
"""
import numpy as np

from openquake.hazardlib.gsim.base import IPE, CoeffsTable
from openquake.hazardlib import const
from openquake.hazardlib.imt import MMI

class AtkinsonEtAl2014Cal(IPE):
    """
    Implements the Intensity Prediction Equation of Atkinson, Worden and 
    Wald (2014). Intensity prediction equations for North America, Bull. 
    Seismol. Soc. Am. 104, doi: 10.1785/0120140178.


    This class implements the version for California neglecting
    site amplification
    """
    #: The GMPE is derived from induced earthquakes
    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.ACTIVE_SHALLOW_CRUST

    #: Supported intensity measure types are peak ground acceleration
    #: and peak ground velocity
    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([
        MMI,
    ])

    #: Supported intensity measure component is not considered for IPEs, so
    #: we assume equivalent to 'average horizontal'
    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = const.IMC.AVERAGE_HORIZONTAL

    #: Supported standard deviation types is total.
    DEFINED_FOR_STANDARD_DEVIATION_TYPES = set([
        const.StdDev.TOTAL
    ])

    #: No required site parameters (in the present version)
    REQUIRES_SITES_PARAMETERS = set()

    #: Required rupture parameters are magnitude (ML is used)
    REQUIRES_RUPTURE_PARAMETERS = set(('mag', ))

    #: Required distance measure is rupture distance
    REQUIRES_DISTANCES = set(('rhypo',))

    def get_mean_and_stddevs(self, sites, rup, dists, imt, stddev_types):
        """
        See :meth:`superclass method
        <.base.GroundShakingIntensityModel.get_mean_and_stddevs>`
        for spec of input and result values.
        """
        C = self.COEFFS[imt]
        mean = (self._compute_magnitude_term(C, rup.mag) +
                self._compute_distance_term(C, dists.rhypo, rup.mag))
        stddevs = self._get_stddevs(stddev_types)
        return mean, stddevs

    def _compute_magnitude_term(self, C, mag):
        """
        Returns the magnitude scaling term
        """
        return C["c1"] + (C["c2"] * mag)

    def _compute_distance_term(self, C, rhypo, mag):
        """
        Returns the distance scaling term
        """
        R = np.sqrt(rhypo**2 + 14.**2)
        
        B = np.zeros_like(R)
        R50 = np.log10(R / 50.)
        idx = R50 > 0
        B[idx] = R50[idx]
        
        return C["c3"]*np.log10(R) + C["c4"]*R + C["c5"]*B + C["c6"]*mag*np.log10(R)

    def _get_stddevs(self, stddev_types):
        """
        Returns the total standard deviation
        """
        stddevs = []
        for stddev_type in stddev_types:
            assert stddev_type in self.DEFINED_FOR_STANDARD_DEVIATION_TYPES
            if stddev_type == const.StdDev.TOTAL:
                sigma = 0.15
                stddevs.append(sigma)
        return stddevs

    COEFFS = CoeffsTable(sa_damping=5, table="""
    IMT     c1     c2      c3     c4    c5     c6
    mmi  0.309  1.864  -1.672  -0.00219  1.77  -0.383
    """)

class AtkinsonEtAl2014CEUS(AtkinsonEtAl2014Cal):
    """
    This class implements the version for the CEUS neglecting
    site amplification
    """
    #: Required distance measure is hypocentral distance
    REQUIRES_DISTANCES = set(('rhypo',))

    def get_mean_and_stddevs(self, sites, rup, dists, imt, stddev_types):
        """
        See :meth:`superclass method
        <.base.GroundShakingIntensityModel.get_mean_and_stddevs>`
        for spec of input and result values.
        """
        C = self.COEFFS[imt]
        mean = (self._compute_magnitude_term(C, rup.mag) +
                self._compute_distance_term(C, dists.rhypo, rup.mag) +
                self._compute_ceus_term(C, dists.rhypo, rup.hypo_depth))
        stddevs = self._get_stddevs(stddev_types)
        return mean, stddevs

    def _compute_ceus_term(self, C, rhypo, hypo_depth):
        """
        Returns the term to adjust California MMI to CEUS MMI
        """  
        repi = np.sqrt(rhypo**2 - hypo_depth**2)
        
        R1 =  np.ones_like(repi) * 150.
        idx = repi < R1
        R1[idx] = repi[idx]
        R_tmp = 0.8 * np.log10(R1/50.)
        
        R2 = np.zeros_like(repi)
        idx = R_tmp > R2
        R2[idx] = R_tmp[idx]
                
        return 0.7 + 0.001 * repi + R2
