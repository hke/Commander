#================================================================================
#
# Copyright (C) 2020 Institute of Theoretical Astrophysics, University of Oslo.
#
# This file is part of Commander3.
#
# Commander3 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Commander3 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Commander3. If not, see <https://www.gnu.org/licenses/>.
#
#================================================================================

import h5py
import commander_tools.tod_tools.commander_instrument as inst
import os
import numpy as np
import healpy as hp
import random
import math

class hfi(object):   
 
    freqs = [100, 143, 217, 353, 545, 857]
    dets = {100:['1a', '1b', '2a', '2b', '3a', '3b', '4a', '4b'], 143:['1a', '1b', '2a', '2b', '3a', '3b', '4a', '4b', '5', '6', '7'], 217:['1', '2', '3', '4', '5a', '5b', '6a', '6b', '7a', '7b', '8a', '8b'], 353:['1', '2', '3a', '3b', '4a', '4b', '5a', '5b', '6a', '6b', '7', '8'], 545:['1', '2', '4'], 857:['1', '2', '3', '4']}
    npsi = 4096
    ntodsigma = 100
    nsides = {100:1024, 143:1024, 217:1024, 353:1024, 545:2048, 857:2048}
    #compression arrays 
    huffman = ['huffman', {'dictNum':1}]
    huffTod = ['huffman', {'dictNum':2}]
    rice = ['rice', {'k':0}]
    psiDigitize = ['digitize', {'min':0, 'max':2*np.pi,'nbins':npsi}]
    todDtype = ['dtype', {'dtype':'int32'}]

#    ['100-1a': , '100-1b':, '100-2a':, '100-2b':, '100-3a':, '100-3b':, '100-4a':, '100-4b':, 
#                    '143-1a':, '143-1b':, '143-2a':, '143-2b':, '143-3a':, '143-3b':, '143-4a':, '143-4b':, '143-5':, '143-6':, '143-7':, '143-8':,
#                    '217-1':, '217-2':, '217-3':, '217-4':, '217-5a':, '217-5b':, '217-6a':, '217-6b':, '217-7a':, '217-7b':, '217-8a':, '217-8b':, 
#                    '353-1':, '353-2':, '353-3a':, '353-3b':, '353-4a':, '353-4b':, '353-5a':, '353-5b':, '353-6a':, '353-6b':, '353-7':, '353-8':,
#                    '545-1':, '545-2':, '545-3':, '545-4':,
#                    '857-1':, '857-2':, '857-3':, '857-4':   ]

    # nu_eff from https://wiki.cosmos.esa.int/planck-legacy-archive/index.php/Spectral_response
    cent_freqs = {'100-1a': 100.28, '100-1b':100.87, '100-2a':101.34, '100-2b':101.19, '100-3a':101.64, '100-3b':100.63, '100-4a':101.77, '100-4b':101.91, 
                    '143-1a':141.71, '143-1b':142.29, '143-2a':141.79, '143-2b':142.5, '143-3a':140.51, '143-3b':141.63, '143-4a':142.71, '143-4b':142.19, '143-5':144.24, '143-6':143.0, '143-7':144.46, '143-8':143.55, 
                    '217-1':222.817, '217-2':223.231, '217-3':223.116, '217-4':222.717, '217-5a':220.421, '217-5b':220.655, '217-6a':220.619, '217-6b':220.619, '217-7a':220.766, '217-7b':220.332, '217-8a':220.51, '217-8b':220.712, 
                    '353-1':360.289, '353-2':360.866, '353-3a':359.59, '353-3b':359.65, '353-4a':362.224, '353-4b':362.212, '353-5a':358.73, '353-5b':358.84, '353-6a':359.91, '353-6b':356.06, '353-7':363.35, '353-8':365.1,
                    '545-1':559.83, '545-2':556.05, '545-3':557.4, '545-4':556.85,
                    '857-1':866.05, '857-2':860.55, '857-3':864.92, '857-4':854.75}

    #psi_uv from https://www.aanda.org/articles/aa/full_html/2016/10/aa25818-15/T5.html
    #mbangs = 

    
    def __init__():
        return


    def init_extra_flags(flagfile):
        flags = {}

        flagdata = open(flagfile)

        for line in flagdata.readlines():
            if line[0] == '#':
                continue
            if(len(line) > 1):
                what, start, stop = line.split(' ')
                if what not in flags.keys():
                    flags[what] = []
                flags[what].append([float(start)*1e9, float(stop)*1e9])

        flagdata.close()

        for entry in flags.keys():
            flags[entry] = np.array(flags[entry])

        return flags 

    def get_extra_flags(times, det, flags):

        outFlags = np.zeros(len(times), dtype=np.uint8)

        for key in flags.keys():
            if det not in key:
                continue

            array = flags[key]

            reduced = array[np.logical_and(array[:,0] < times[-1], array[:,1] > times[0])]
            for entry in reduced:
                outFlags[np.logical_and(times > entry[0],times < entry[1])] = 128
        return outFlags

    @staticmethod
    def instrument_filename(version):
        return 'HFI_instrument_v' + str(version) + '.h5'

    @staticmethod
    def ring_outer_product(theta, phi):
        outAng = [0, 0]
        nsamps = min(100, len(theta))
        pair1 = random.sample(range(len(theta)), nsamps)
        pair2 = random.sample(range(len(theta)), nsamps)
        vecs = hp.pixelfunc.ang2vec(theta, phi)
        for a1, a2 in zip(pair1, pair2):
            crossP = np.cross(vecs[a1], vecs[a2])
            if(crossP[0] < 0):
                crossP *= -1
            theta1, phi1 = hp.vec2ang(crossP)
            if(not math.isnan(theta1) and not math.isnan(phi1)):
                outAng[0] += theta1/nsamps
                outAng[1] += phi1/nsamps
            else:
                outAng[0] += outAng[0]/nsamps
                outAng[1] += outAng[1]/nsamps
        return outAng

    def freq2dets(self, freq):
        return self.dets[freq]

    @staticmethod
    def verify_instrument_file(outDir, version):
        f = inst.commander_instrument(outDir, hfi.instrument_filename(version), version, 'r')

        if version < 1:
            raise ValueError("Version number must be positive")
        
        if version == 1:
            print("Should check the version here")

        if version > 1:
            raise ValueError("Version " + str(version) + " of HFI instrument file has not yet been defined.")
