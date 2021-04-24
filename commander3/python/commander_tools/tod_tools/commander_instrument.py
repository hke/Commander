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
import os
import numpy as np

class commander_instrument:

    def __init__(self, path, fileName, version, mode):
        self.version = version
        self.outPath = path
        self.h5file = h5py.File(os.path.join(path, fileName), mode)
        
    #Raw field writing
    def add_field(self, fieldName, data):
        self.h5file.create_dataset(fieldName, data=data)

    #Write a matrix field with info about the column names
    def add_matrix(self, fieldName, data, columnInfo):
        dims = np.shape(data)

        if(len(dims) != 2):
            raise TypeError('Call to add matrix with an object of shape ' + str(shape))

        if(dims[0] != len(columnInfo)):
            if(dims[1] == len(columnInfo)):
                data = np.transpose(data)
            else:
                raise ValueError('Data is shape ' + str(shape) + ' but column headers have length ' + str(len(columnInfo)))

        self.h5file.create_dataset(fieldName, data=data)
        self.h5file[fieldName].attrs['index'] = columnInfo

    
    #add bandpass data for detector det
    def add_bandpass(self, det, freqs, response):
        if len(freqs) != len(response):
            raise ValueError('Bandpass frequency x and f(x) must be the same length, len(freqs): ' + str(len(freqs)) + ' len(f(x)): ' + str(len(response)))
        det = self.parse_det(det)

        self.add_field(det + '/bandpassx', freqs)
        self.add_field(det + '/bandpass', response)
       
    def add_alms(self, det, almType, lmax, mmax, T, E, B):
       
        if(mmax > lmax):
            mmax = lmax
 
        nalm = self.nalm(lmax, mmax)
        if (len(T) != nalm) or (E is not None and len(E) != nalm) or (B is not None and len(B) != nalm):
            errstr = 'Length of alms is wrong for lmax=' + str(lmax) + ', mmax=' +str(mmax) + '. nalm = ' + str(nalm) + ' Lengths: ' + str(len(T))
            if(E is not None and B is not None):
                errstr += ', ' + str(len(E)) + ', ' + str(len(B))

            raise ValueError(errstr)
               
        det = self.parse_det(det) 

        self.add_field(det + '/' + almType + '/T', T)
        #polarized experiment
        if(E is not None and B is not None):
            self.add_field(det + '/' + almType + '/E', E)
            self.add_field(det + '/' + almType + '/B', B)

        self.add_field(det + '/'+almType + 'mmax', [mmax])
        self.add_field(det + '/'+almType + 'lmax', [lmax])

    #currently does nothing but maybe it will in the future, who knowssss
    def finalize(self):
        self.add_field('common/version', self.version)

    def parse_det(self, det):
        if type(det) is not str:
            det = str(det).zfill(3)
        return det

    def nalm(self, lmax, mmax):
        return lmax**2 + 2*lmax +1
