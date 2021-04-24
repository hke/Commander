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

import numpy as np


# k = 0 means autoselect k
def encode(data, k=0):

    assert(np.issubdtype(data.dtype, np.integer))
    #bitstream = bitarray()
    bitstream = ''

    offset = int(np.mean(data))
    data = data - offset
    data = np.insert(data, 0, offset)

    if(k < 1): #autoselect k
        k = 1
        mean = np.mean(np.abs(data[1:]))
        if mean > 0:
            k = int(np.log2(mean))

    #print(offset, data, len(data), np.mean(np.abs(data)), pow(2, k), k)

    m = pow(2,k)        

    for entry in data:

        negBit = '1'
        if entry > 0:
            negBit = '0'

        entry = abs(entry)

        q = int(entry/m)
        for i in range(q):    
            #bitstream.append(1)
            bitstream += '1'
        #bitstream.append(0)
        bitstream += '0' + negBit
        #bitstream.append(negBit)   
 
        binary = bin(entry)[2:].zfill(k)
        bits = binary[-k:]
        #for bit in bits:
        #    bitstream.append(int(bit))
        bitstream += bits

        #print(entry, bitstream, q, binary, bits, negBit)

    padding = 8 - len(bitstream) % 8
    bitstream += padding*'0'
    

    b = bytearray()
    for i in range(0, len(bitstream), 8):
        byte = bitstream[i:i+8]
        b.append(int(byte, 2))

    #return bitstream, k
    return b, k

def decode(data, k):
    data = bytearray(data)
    binary_txt = ''.join(bin(i)[2:].rjust(8,'0') for i in data)

    k = int(k[:])
    m = pow(2, k)
    outData = []
    i = 0

    while i<len(binary_txt):

        if(len(binary_txt) - i < 8 and '1' not in binary_txt[i:]):
            break

        num = 0
        negBit = 1

        while(binary_txt[i] == '1'):
            num += m
            i += 1
        i += 1 #skip the placeholder 0

        if(binary_txt[i] == '1'): #read the negative bit
            negBit = -1

        i += 1
  
        num += int(binary_txt[i:i+k], 2)
        i += k

        outData.append(num*negBit)

    return outData
