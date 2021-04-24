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

class lfi(object):   
 
    freqs = [30, 44, 70]
    horns = {30:[27, 28], 44:[24, 25, 26], 70:[18, 19, 20, 21, 22, 23]}
    hornTypes = ['M', 'S']
    diodeTypes = {'M':['00', '01'], 'S':['10', '11']}
    npsi = 4096
    ntodsigma = 1
    nsides = {30:512, 44:512, 70:1024}
    #compression arrays 
    huffman = ['huffman', {'dictNum':1}]
    huffTod = ['huffman', {'dictNum':2}]
    psiDigitize = ['digitize', {'min':0, 'max':2*np.pi,'nbins':npsi}]
    todDtype = ['dtype', {'dtype':'f4'}]
    todSigma = ['sigma', {'sigma0':None, 'nsigma':ntodsigma, 'offset':None}] 
    #fwhm, elipticity and psi_ell from https://www.aanda.org/articles/aa/full_html/2016/10/aa25809-15/T6.html
    fwhms = {'18M':13.44, '18S':13.5, '19M':13.14, '19S':13.07, '20M':12.84, '20S':12.84, '21M':12.77, '21S':12.87, '22M':12.92, '22S':12.97, '23M':13.35, '23S':13.36, '24M':23.18, '24S':23.04, '25M':30.23, '25S':30.94, '26M':30.29, '26S':30.64, '27M':32.02, '27S':33.11, '28M':33.1, '28S':33.09}

    elips = {'18M':1.23, '18S':1.27, '19M':1.25, '19S':1.28, '20M':1.27, '20S':1.29, '21M':1.28, '21S':1.29, '22M':1.27, '22S':1.28, '23M':1.23, '23S':1.28, '24M':1.39, '24S':1.34, '25M':1.19, '25S':1.19, '26M':1.19, '26S':1.19, '27M':1.37, '27S':1.38, '28M':1.37, '28S':1.37}

    psis = {'18M':85, '18S':86, '19M':78, '19S':79, '20M':71, '20S':72, '21M':107, '21S':106, '22M':101, '22S':101, '23M':92, '23S':92, '24M':89, '24S':89, '25M':114, '25S':117, '26M':62, '26S':61, '27M':101, '27S':101, '28M':78, '28S':78}

    b_effs = {'18M':0.9921, '18S':0.9887, '19M':0.9883, '19S':0.9898, '20M':0.9885, '20S':0.9881, '21M':0.9894, '21S':0.9882, '22M':0.9916, '22S':0.9915, '23M':0.9926, '23S':0.9919, '24M':0.9972, '24S':0.9973, '25M':0.9975, '25S':0.9976, '26M':0.9974, '26S':0.9977, '27M':0.9904, '27S':0.9889, '28M':0.9907, '28S':0.9879}

    cent_freqs = {'18M':70.4, '18S':70.4, '19M':70.4, '19S':70.4, '20M':70.4, '20S':70.4, '21M':70.4, '21S':70.4, '22M':70.4, '22S':70.4, '23M':70.4, '23S':70.4, '24M':44.1, '24S':44.1, '25M':44.1, '25S':44.1, '26M':44.1, '26S':44.1, '27M':28.4, '27S':28.4, '28M':28.4, '28S':28.4}

    #psi_uv from https://www.aanda.org/articles/aa/full_html/2016/10/aa25818-15/T5.html
    mbangs = {27:-22.46, 28:22.45, 24:0.01, 25:-113.23, 26:113.23, 18:22.15, 19:22.4, 20:22.38, 21:-22.38, 22:-22.34, 23:-22.08}

    
    def __init__():
        return

    @staticmethod
    def instrument_filename(version):
        return 'LFI_instrument_v' + str(version) + '.h5'

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

    def horn2freq(self, horn):
         return

    def freq2horns(self, freq):
        return self.horns[freq]

    @staticmethod
    def verify_instrument_file(outDir, version):
        f = inst.commander_instrument(outDir, lfi.instrument_filename(version), version, 'r')

        if version <= 3:
            raise ValueError("Versions of LFI instrument files <= 3 are lost to time")
        if version == 4:
            print("Should check the version here")

        if version > 4:#new beam information 
            assert(f.h5file['/28M/beammmax'][()] == 14)
            assert(f.h5file['/28M/beamlmax'][()] == 3000)

        if version > 5:
            raise ValueError("Version " + str(version) + " of LFI instrument file has not yet been defined.")

    @staticmethod
    def mbang(horn):
        mbangs = {27:-22.46, 28:22.45, 24:0.01, 25:-113.23, 26:113.23, 18:22.15, 19:22.4, 20:22.38, 21:-22.38, 22:-22.34, 23:-22.08}
        if len(str(horn)) == 3:
            horn = int(str(horn[:-1]))
        return mbangs[int(horn)]

    @staticmethod
    def complex2realAlms(data, mmax):
        lmax = lfi.getLmax(len(data), mmax)
        outData = np.zeros((lmax+1)**2)

        for l in range(0, lmax):
            for m in range(0, mmax):
                if(m > l):
                    continue
                #TODO: figure this out
                scaling = np.sqrt(2)
                if(m == 0):
                    scaling = 1
                healpixI = hp.sphtfunc.Alm.getidx(lmax, l, m)
                outI = lfi.getOutidx(l, m)
                outJ = lfi.getOutidx(l, -1*m)
                outData[outI] = np.real(data[healpixI]) * scaling
                if(m != 0):
                    outData[outJ] = np.imag(data[healpixI]) * scaling

        return outData

    @staticmethod
    def getLmax(N, mmax):
        return int((2.0*(N - 1.0)/mmax + mmax -1)/(2.0*(1.0 + 1.0/mmax)))

    @staticmethod
    def getOutidx(l, m):
        return l**2 + l + m

