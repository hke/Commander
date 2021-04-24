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
import commander_tools.tod_tools.huffman as huffman
import commander_tools.tod_tools.rice as rice
import healpy as hp
import numpy as np
import multiprocessing as mp
import os
import sys

class commander_tod:

    def __init__(self, outPath, name, version=None, dicts=None, overwrite=False):
        self.outPath = outPath
        self.filelists = dicts
        self.version = version
        self.name = name
        #TODO: something with the version number
        self.overwrite = overwrite

    #initilizes a file for a single od
    def init_file(self, freq, od, mode='r'):
        self.huffDict = {}
        self.raggedDict = {}
        self.attrDict = {}
        self.encodings = {}
        self.pids = {}

        self.od = od
        self.freq = freq
        self.outName = os.path.join(self.outPath, self.name+ '_' + str(freq).zfill(3) + '_' + str(od).zfill(6) + '.h5')

        self.exists = False
        if os.path.exists(self.outName):
            self.exists = True
        if mode == 'w':
            if self.exists and self.overwrite:
                os.remove(self.outName)
            try:
                self.outFile = h5py.File(self.outName, 'a')

                if self.exists and not self.overwrite:
                    for pid in self.load_field('/common/pids'):
                        loadBalance = self.load_field('/' + str(pid).zfill(6) + '/common/load')
                        self.pids[pid] = str(float(loadBalance[0])) + ' ' + str(float(loadBalance[1]))
            except (KeyError, OSError):
                if(hasattr(self, 'outFile')):
                    self.outFile.close()
                os.remove(self.outName)
                self.exists = False
                self.outFile = h5py.File(self.outName, 'a')
        
        if mode == 'r':
            if not self.exists:
                raise OSError('Cannot find file ' + self.outName)
            self.outFile = h5py.File(self.outName, 'r')


    # Write a possibly compressed, possibly ragged matrix
    def add_matrix(self, fieldName, data, columnInfo, compression=None):

        writeField = True
        dims = np.shape(data)

        if(len(dims) != 2):
            raise TypeError('Call to add matrix with an object of shape ' + str(shape))

        if(dims[0] != len(columnInfo)):
            if(dims[1] == len(columnInfo)):
                data = np.transpose(data)
            else:
                raise ValueError('Data is shape ' + str(shape) + ' but column headers have length ' + str(len(columnInfo)))  

        if compression == None or compression == []:
            print(type(data))
            if type(data) == np.ndarray:
                self.outFile.create_dataset(fieldName, data=data)
            else: #ragged array
                dt = h5py.vlen_dtype(np.dtype(data[0]))
                dset = f.create_dataset(fieldName, dims[0], dtype=dt)
                for i,vec in zip(range(len(data)), data):
                    dset[i] = data
        else: #compressed array
            for i, vec in zip(range(len(data)), data):
                self.add_field(fieldName+':' + str(i), data[i], compression)

        self.add_attribute(fieldName, 'index', columnInfo)
        self.add_attribute(fieldName, 'matrix', True)

    #Single field write
    def add_field(self, fieldName, data, compression=None):
        writeField = True
        if(compression is not None and compression is not []):
            compInfo = ''
            if(len(compression) == 2 and type(compression[0]) == str):
                #catch case with only one compression argument not in array
                compression = [compression]
            for compArr in compression:
                compInfo += compArr[0] + ' '
                if compArr[0] == 'dtype':
                    data=np.array(data, dtype=compArr[1]['dtype'])

                elif compArr[0] == 'sigma':
                    data = np.int32(compArr[1]['nsigma'] * (data-compArr[1]['offset'])/(compArr[1]['sigma0']))
                    metaName = '/common/n' + fieldName.split('/')[-1] + 'sigma'
                    self.encodings[metaName] = compArr[1]['nsigma']
 
                    self.add_attribute(fieldName, 'nsigma', compArr[1]['nsigma'])
                    self.add_attribute(fieldName, 'sigma0', compArr[1]['sigma0'])
                    self.add_attribute(fieldName, 'offset', compArr[1]['offset'])

                elif compArr[0] == 'digitize':
                    bins = np.linspace(compArr[1]['min'], compArr[1]['max'], num = compArr[1]['nbins'])
                    data = np.digitize(data, bins)
                    metaName = '/common/n' + fieldName.split('/')[-1]
                    self.add_encoding(metaName, compArr[1]['nbins'])
                    self.add_attribute(fieldName, 'min', compArr[1]['min'])
                    self.add_attribute(fieldName, 'max', compArr[1]['max'])
                    self.add_attribute(fieldName, 'nbins', compArr[1]['nbins'])
                

                elif compArr[0] == 'huffman': #differenced huffman
                    dictNum = compArr[1]['dictNum']
                    if dictNum not in self.huffDict.keys():
                        self.huffDict[dictNum] = {}
                    delta = np.diff(data)
                    delta = np.insert(delta, 0, data[0])
                    self.huffDict[dictNum][fieldName] = delta
                    #print("adding " + fieldName + " to dict, contents ", delta[delta != 0], data[data != 0])
                    self.add_attribute(fieldName, 'huffmanDictNumber', dictNum)
                    writeField = False 

                elif compArr[0] == 'rice': #rice encoding
                    k = compArr[1]['k']
                    data, k = rice.encode(data, k)
                    data = np.void(bytes(data))
                    self.add_encoding(fieldName.rsplit('/',1)[0] + '/riceK', k)
                    self.add_attribute(fieldName, 'riceK', k)

                else:
                    raise ValueError('Compression type ' + compArr[0] + ' is not a recognized compression')
            self.add_attribute(fieldName, 'compression', compInfo)
            #print("adding " + compInfo + ' to ' + fieldName)

        if ':' in fieldName:
            writeField = False
            dictName, index = fieldName.rsplit(':', 1)
            if dictName not in self.raggedDict.keys():
                self.raggedDict[dictName] = {}

            self.raggedDict[dictName][index] = data

        if writeField:
            try:
                self.outFile.create_dataset(fieldName, data=data)
            except OSError as e:
                if self.overwrite:
                    del self.outFile[fieldName]
                    self.outFile.create_dataset(fieldName, data=data)
                else:
                    raise OSError(e)
            for attr in self.attrDict.copy().keys():
                fName, attrName = attr.split('@')
                if fName == fieldName:
                    self.add_attribute(fieldName, attrName, self.attrDict.pop(attr))

    def add_attribute(self, fieldName, attrName, data):
        try:
            self.outFile[fieldName].attrs[attrName] = data
        except KeyError as k:
            self.attrDict[fieldName + '@' + attrName] = data   
 
    def add_encoding(self, encoding, value):
        if encoding not in self.encodings.keys():
            self.encodings[encoding] = value
        else:
            if self.encodings[encoding] != value:
               print('Warning: Inconsistant encoding value ' + encoding + ' is set to ' + str(self.encodings[encoding]) + ' but wants to be ' + str(value))

    def finalize_file(self):

        if(not self.exists or self.overwrite):
            for encoding in self.encodings.keys():
                self.add_field(encoding, [self.encodings[encoding]])
                #print('adding ' + encoding + ' to file ' + self.outName)

            self.add_field('/common/version', self.version)
            self.add_field('/common/pids', list(self.pids.keys()))

        if self.filelists is not None:
            for pid in self.pids.keys():
                self.filelists[self.freq]['id' + str(pid)] = str(pid) + ' "' + os.path.abspath(self.outName) + '" ' + '1 ' + self.pids[pid] + '\n'       
 
        return

    def finalize_chunk(self, pid, loadBalance=[0,0]):
        if len(loadBalance) != 2:
            raise ValueError('Load Balancing numbers must be length 2')
        for key in self.huffDict.keys():
            h = huffman.Huffman()
            h.GenerateCode(list(self.huffDict[key].values()))
            huffArray = np.append(np.append(np.array(h.node_max), h.left_nodes), h.right_nodes)
            numStr = str(key)
            if(key == 1):
                numStr = ''

            self.add_field('/' + str(pid).zfill(6) + '/common/hufftree' + numStr, huffArray)
            self.add_field('/' + str(pid).zfill(6) + '/common/huffsymb' + numStr, h.symbols)
            #with np.printoptions(threshold=np.inf):
            #    print(huffArray, len(huffArray), len(h.symbols))
            for field in self.huffDict[key].keys():
                if ':' not in field:
                    self.add_field(field, np.void(bytes(h.byteCode(self.huffDict[key][field]))))
                    #print(self.huffDict[key][field], len(h.byteCode(self.huffDict[key][field])))
                else:
                    fieldName, index = field.rsplit(':', 1)
                    self.raggedDict[fieldName][index] = np.void(bytes(h.byteCode(self.huffDict[key][field])))


        # handle the compressed ragged matrices
        for field in self.raggedDict.keys():
            ragged = False 
            length = len(list(self.raggedDict[field].values())[0])
            for entry in self.raggedDict[field].keys():
                if len(self.raggedDict[field][entry]) != length or length < 1:
                    # np.void objects seem to have length 0
                    ragged = True

            if ragged == False:
                dataType = np.dtype(list(self.raggedDict[field].values())[0][0])
                array = np.ndarray((len(self.raggedDict[field].keys()), length), dtype=dataType)

                for entry in self.raggedDict[field].keys():
                    array[int(entry)] = self.raggedDict[field][entry] 

                self.outFile.create_dataset(field, data=array)

            else:
                dt = h5py.special_dtype(vlen=np.dtype('uint8'))
                dset = self.outFile.create_dataset(field, len(list(self.raggedDict[field].keys())), dtype=dt)
                for entry in self.raggedDict[field].keys():
                    dset[int(entry)] = np.frombuffer(self.raggedDict[field][entry], dtype='uint8')
 
        for attr in self.attrDict.copy():
            field, attrName = attr.split('@', 1)
            if ':' in field:
                fieldName, index = field.split(':', 1)
            else:
                fieldName = field

            self.add_attribute(fieldName, attrName, self.attrDict.pop(attr))
                
 
        if len(self.attrDict) > 0:
            print("Attributes left unassigned at chunk end", self.attrDict)

 
        self.raggedDict = {}
        self.huffDict = {}
        self.add_field('/' + str(pid).zfill(6) + '/common/load', loadBalance)
        self.pids[pid] = str(float(loadBalance[0])) + ' ' + str(float(loadBalance[1]))

    def compute_version(self):
        return

    def make_filelists(self):
        for freq in self.filelists.keys():
            outfile = open(os.path.join(self.outPath, 'filelist_' + str(freq) + '.txt'), 'w')
            outfile.write(str(len(self.filelists[freq])) + '\n')
            for buf in self.filelists[freq].values():
                #print(buf, len(buf))
                outfile.write(buf)

            outfile.close()

        return

    #File Reading Functions
    def load_field(self, fieldName):
        try:
            compStr = self.outFile[fieldName].attrs['compression']
        except KeyError:
            return self.outFile[fieldName]
        
        return self.decompress(fieldName, compression=compStr)       

    def load_all_fields(self):
        return

    def read_across_files(self, fieldName):
        return

    def decompress(self, field, compression=''):
        comps = compression.split(' ')
        data = self.outFile[field][:]
        print(data, len(data), len(data[0]), len(data[1]), comps)
        ndim = 1
        try:
            matrix = self.outFile[field].attrs['matrix']
        except KeyError:
            ndim = 1

        if(matrix):
            ndim = len(data) 
        
        for i in range(ndim):
    
            if ndim > 1:
                dataBuf = data[i]
            else:
                dataBuf = data

            for comp in comps[::-1]: # apply the filters in the reverse order
                if comp == '':
                    #residual from str.split()
                    pass
                elif comp == 'dtype':
                    pass
                elif comp == 'sigma':
                    sigma0 = self.outFile[field].attrs['sigma0']
                    nsigma = self.outFile[field].attrs['nsigma']
                    offset = self.outFile[field].attrs['offset']
                    dataBuf = dataBuf * sigma0/nsigma + offset                

                elif comp == 'digitize':
                    nbins = self.outFile[field].attrs['nbins']
                    nmin = self.outFile[field].attrs['min']
                    nmax = self.outFile[field].attrs['max']

                    bins = np.linspace(nmin, nmax, num = nbins)
                    dataBuf = bins[dataBuf]

                elif comp == 'huffman':
                    pid = field.split('/')[1]
                    try:
                        huffNum = str(self.outFile[field].attrs['huffmanDictNumber'])
                    except KeyError:
                        huffNum = ""
                    if huffNum == '1':
                        huffNum = ''
                    huffTree = self.load_field('/' + pid + '/common/hufftree' + huffNum)
                    huffSymb = self.load_field('/' + pid + '/common/huffsymb' + huffNum)
                    h = huffman.Huffman(tree=huffTree, symb=huffSymb)
                    dataBuf = h.Decoder(np.array(dataBuf))

                elif comp == 'rice':
                    k = self.load_field(field.rsplit('/', 1)[0] + '/riceK')
                    dataBuf = np.array(rice.decode(np.array(dataBuf), k))
                    offset = dataBuf[0]
                    dataBuf += offset
                    dataBuf = dataBuf[1:]                
                    
                else:
                    raise ValueError('Decompression type ' + comp + ' is not a recognized operation')
       
            if ndim > 1:
                data[i] = dataBuf
            else:
                data = dataBuf
 
        return data


