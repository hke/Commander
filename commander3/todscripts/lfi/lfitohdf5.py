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

from commander_tools.tod_tools.lfi import lfi
from commander_tools.tod_tools import commander_tod as tod
import argparse
import multiprocessing as mp
import os
import numpy as np
import math
from astropy.io import fits
import healpy as hp        
import sys
import random
import h5py

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('planck_dir', type=str, action='store', help='path to the legacy planck data in hdf format')
    #/mn/stornext/d16/cmbco/bp/data

    parser.add_argument('--gains-dir', type=str, action='store', help='path to a directory with the initial gain estimates', default='/mn/stornext/d16/cmbco/bp/data/npipe_gains')

    parser.add_argument('--velocity-file', type=str, action='store', help='path to a file with the satelite velocities', default='/mn/stornext/d16/cmbco/bp/data/auxiliary_data/satellite_velocity.fits')

    parser.add_argument('--position-file', type=str, action='store', help='path to the on disk satellite position file', default='/mn/stornext/d16/cmbco/bp/data/auxiliary_data/planck_xyz.txt')

    parser.add_argument('--rimo', type=str, action='store', help='path to on disk rimo file', default='/mn/stornext/d14/bp/data/auxiliary_data/LFI_RIMO_R3.31.fits')

    parser.add_argument('--out-dir', type=str, action='store', default=os.getcwd(), help='path to output data structure you want to generate')

    parser.add_argument('--num-procs', type=int, action='store', default=1, help='number of processes to use')

    parser.add_argument('--freqs', type=int, nargs='+', default=lfi.freqs, help='which lfi frequencies to operate on')

    parser.add_argument('--ods', type=int, nargs=2, default=[91, 1604], help='the operational days to operate on')

    parser.add_argument('--no-compress', action='store_true', default=False, help='Produce uncompressed data output')

    parser.add_argument('--no-compress-tod', action='store_true', default=False, help='should we compress the tod field')

    parser.add_argument('--restart', action='store_true', default=False, help="restart from a previous run that didn't finish")

    parser.add_argument('--produce-filelist', action='store_true', default=False, help='force the production of a filelist even if only some files are present')

    parser.add_argument('--differenced-data', action='store_true', default=False, help='store the differenced data produced by the DPC instead of the L1 data')

    in_args = parser.parse_args()

    in_args.version = 5
    if(in_args.no_compress):
        in_args.version = 4
    if(not in_args.differenced_data):
        in_args.version += 2
    
    random.seed()

    os.environ['OMP_NUM_THREADS'] = '1'

    pool = mp.Pool(processes=in_args.num_procs)
    manager = mp.Manager()

    ods = range(in_args.ods[0], in_args.ods[1], 1)

    manager = mp.Manager()
    dicts = {30:manager.dict(), 44:manager.dict(), 70:manager.dict()}

    comm_tod = tod.commander_tod(in_args.out_dir, 'LFI', in_args.version, dicts, not in_args.restart)

    x = [[pool.apply_async(make_od, args=[comm_tod, freq, od, in_args]) for freq in in_args.freqs] for od in ods]

    for res1 in np.array(x):
        for res in res1:
            res.get()

    pool.close()
    pool.join()

    if ((in_args.ods[0] == 91 and in_args.ods[1] == 1604) or in_args.produce_filelist) :
        comm_tod.make_filelists()
        #write file lists 

def make_od(comm_tod, freq, od, args):

    print(freq, od)

    nside = lfi.nsides[freq]

    comm_tod.init_file(freq, od, mode='w')

    if(args.restart and comm_tod.exists):
        comm_tod.finalize_file()
        print('Skipping existing file ' + comm_tod.outName)
        return

    rimo = fits.open(args.rimo)

    if args.velocity_file is not None:
        velFile = fits.open(args.velocity_file)

    if args.position_file is not None:
        posArray = np.loadtxt(args.position_file, comments='*').transpose()
        #Julian Date to Modified Julian Date
        posArray[0] -= 2400000.5

    #make common group for things we only read once
    #polang, mbeamang, nside, fsamp
    prefix = '/common'
    rimo_i = np.where(rimo[1].data.field('detector').flatten() == 'LFI' + str(lfi.horns[freq][0]) + 'M')

    #sampling frequency
    fsamp = rimo[1].data.field('f_samp')[rimo_i]
    comm_tod.add_field(prefix + '/fsamp', fsamp)

    #nside
    comm_tod.add_field(prefix + '/nside', [nside])

    detNames = ''
    polangs = []
    mainbeamangs = []
    for horn in lfi.horns[freq]:
        for hornType in lfi.hornTypes:
            rimo_i = np.where(rimo[1].data.field('detector').flatten() == 'LFI' + str(horn) + hornType)

            detNames += str(horn) + hornType + ', '
            polangs.append(math.radians(rimo[1].data.field('psi_pol')[rimo_i]))
            mainbeamangs.append(math.radians(lfi.mbangs[horn]))

    compArr=[lfi.huffman]
    if args.no_compress:
        compArr = None

    #make detector names lookup
    comm_tod.add_field(prefix + '/det', np.string_(detNames[0:-2]))

    diodeNames = 'M:sky00,ref00,sky01,ref01.S:sky10,ref10,sky11,ref11'
    comm_tod.add_field(prefix + '/diodes', np.string_(diodeNames))

    #make polarization angle
    comm_tod.add_field(prefix + '/polang', polangs)
    comm_tod.add_attribute(prefix + '/polang', 'index', detNames[0:-2])

    #make main beam angle
    comm_tod.add_field(prefix + '/mbang', mainbeamangs)
    comm_tod.add_attribute(prefix + '/mbang', 'index', detNames[0:-2])

    try:
        exFile = h5py.File(os.path.join(args.planck_dir, 'L2Data', 'LFI_0' + str(freq) + '_' + str(lfi.horns[freq][0]) + '_L2_002_OD' + str(od).zfill(4) +'.h5'), 'r')
    except (OSError):
        print("Failed to open file " + os.path.join(args.planck_dir, 'L2Data', 'LFI_0' + str(freq) + '_' + str(lfi.horns[freq][0]) + '_L2_002_OD' + str(od).zfill(4) +'.h5'))
        return

    #per pid
    for pid, index in zip(exFile['AHF_info/PID'], range(len(exFile['AHF_info/PID']))):
        startIndex = np.where(exFile['Time/OBT'] > exFile['AHF_info/PID_start'][index])
        endIndex = np.where(exFile['Time/OBT'] > exFile['AHF_info/PID_end'][index])
        #print(pid)
        if len(startIndex[0]) > 0:
            pid_start = startIndex[0][0]
        else:#catch days with no pids
            continue
        if len(endIndex[0]) != 0:
            pid_end = endIndex[0][0]
        else:#catch final pid per od
            pid_end = len(exFile['Time/OBT'])
        if pid_start == pid_end:#catch chunks with no data like od 1007
            print('skipping pid ' + str(pid) + ' because it is empty')
            continue

        obt = exFile['Time/OBT'][pid_start]

        #common fields per pid
        prefix = str(pid).zfill(6) + '/common'

        #time field
        comm_tod.add_field(prefix + '/time', [exFile['Time/MJD'][pid_start], exFile['Time/OBT'][pid_start], exFile['Time/SCET'][pid_start]])
        comm_tod.add_attribute(prefix + '/time','index','MJD, OBT, SCET')

        #length of the tod
        comm_tod.add_field(prefix + '/ntod', [pid_end - pid_start])

        #velocity field
        velIndex = np.where(velFile[1].data.scet > exFile['Time/SCET'][pid_start])[0][0]
        #rotate from ecliptic to galactic 
        r = hp.Rotator(coord=['E', 'G'])
        comm_tod.add_field(prefix + '/vsun', r([velFile[1].data.xvel[velIndex], velFile[1].data.yvel[velIndex], velFile[1].data.zvel[velIndex]])) 
      
        #add some metadata so someone might be able to figure out what is going on 
        comm_tod.add_attribute(prefix + '/vsun','index', '[x, y, z]')
        comm_tod.add_attribute(prefix + '/vsun','coords','galactic')

        #satelite position
        posIndex = np.where(posArray[0] > exFile['Time/MJD'][pid_start])[0][0]
        comm_tod.add_field(prefix + '/satpos', [posArray[1][posIndex], posArray[2][posIndex], posArray[3][posIndex]])
        #add metadata
        comm_tod.add_attribute(prefix + '/satpos','index','X, Y, Z')
        comm_tod.add_attribute(prefix + '/satpos','coords','heliocentric')

        #open per freq npipe gains file if required
        if args.gains_dir is not None and "npipe" in args.gains_dir:#this is a shitty test
            gainsFile = fits.open(os.path.join(args.gains_dir, 'gains_0' + str(freq) + '_iter01.fits'))

        #per detector fields
        for horn in lfi.horns[freq]:
            fileName = h5py.File(os.path.join(args.planck_dir, 'L2Data', 'LFI_0' + str(freq) + '_' + str(horn) + '_L2_002_OD' + str(od).zfill(4) +'.h5'), 'r')
      
            if(not args.differenced_data):
                undiffFile = h5py.File(os.path.join(args.planck_dir, 'L1Data', 'LFI_0' + str(freq) + '_' + str(horn) + '_L1_OD' + str(od).zfill(4) + '.h5'), 'r')
 
            for hornType in lfi.hornTypes:
                #print(horn, hornType)
                prefix = str(pid).zfill(6) + '/' + str(horn) + hornType

                #get RIMO index
                #print(rimo[1].data.field('detector').flatten().shape, rimo[1].data.field('detector').flatten(), 'LFI' +str(horn) + hornType)
                rimo_i = np.where(rimo[1].data.field('detector').flatten() == 'LFI' + str(horn) + hornType)
                
                #make flag data
                flagArray = fileName[str(horn) + hornType + '/FLAG'][pid_start:pid_end]
                
                if (len(flagArray) > 0):
                    comm_tod.add_field(prefix + '/flag', flagArray, compArr)

                #make pixel number
                newTheta, newPhi = r(fileName[str(horn) + hornType + '/THETA'][pid_start:pid_end], fileName[str(horn) + hornType + '/PHI'][pid_start:pid_end])
                pixels = hp.pixelfunc.ang2pix(nside, newTheta, newPhi)

                if len(pixels > 0):
                    #compute average outer product
                    outAng = lfi.ring_outer_product(newTheta, newPhi)
                    comm_tod.add_field(prefix + '/outP', data=outAng)
                    if(args.no_compress):
                        comm_tod.add_field(prefix+'/theta', data=newTheta)
                        comm_tod.add_field(prefix+'/phi', data=newPhi) 
                    comm_tod.add_field(prefix + '/pix', pixels, compArr)


                #make pol angle
                psiArray = fileName[str(horn) + hornType + '/PSI'][pid_start:pid_end] + r.angle_ref(fileName[str(horn) + hornType + '/THETA'][pid_start:pid_end], fileName[str(horn) + hornType + '/PHI'][pid_start:pid_end]) + math.radians(rimo[1].data.field('psi_pol')[rimo_i])
                if(len(psiArray) > 0):
                    psiArray = np.where(psiArray < 0, 2*np.pi + psiArray, psiArray)
                    psiArray = np.where(psiArray >= 2*np.pi, psiArray - 2*np.pi, psiArray)
                    compArray = None
                    if not args.no_compress:
                        compArray = [lfi.psiDigitize, lfi.huffman]            
                    comm_tod.add_field(prefix + '/psi', psiArray, compArray)
                
                #scalars
                
                gain = 1
                #make gain
                if args.gains_dir is not None and "npipe" in args.gains_dir:#this is a shitty test
                    baseGain = fits.getdata(os.path.join(args.gains_dir, 'C0' + str(freq) + '-0000-DX11D-20150209_uniform.fits'),extname='LFI' + str(horn) + hornType)[0][0]
                    gainArr = gainsFile['LFI' + str(horn) + hornType].data.cumulative
                    obtArr = (1e-9 * pow(2,16)) * gainsFile[1].data.OBT
                    gainI = np.where(obtArr <= obt)[0][-1]

                    gain = np.array([1.0/(baseGain * gainArr[gainI])])


                elif args.gains_dir is not None:
                    gainFile = fits.open(os.path.join(args.gains_dir, 'LFI_0' + str(freq) + '_LFI' + str(horn) + hornType + '_001.fits'))
                    gain=1.0/gainFile[1].data.GAIN[np.where(gainFile[1].data.PID == pid)]
                    gainFile.close()
                #TODO: fix this
                if(type(gain) is int or gain.size == 0):
                    gain = [0.06]
         
                #make white noise
                sigma0 = rimo[1].data.field('net')[rimo_i] * math.sqrt(fsamp)

                #make f_knee
                fknee = rimo[1].data.field('f_knee')[rimo_i]

                #make 1/f noise exponent 
                alpha = rimo[1].data.field('alpha')[rimo_i]

                #print(gain, sigma0, fknee, alpha)
                comm_tod.add_field(prefix + '/scalars', np.array([gain, sigma0, fknee, alpha]).flatten())
                comm_tod.add_attribute(prefix + '/scalars','index','gain, sigma0, fknee, alpha')

                #make psd noise
               
                #make tod data
                if(args.differenced_data):
                    tod = fileName[str(horn) + hornType +'/SIGNAL'][pid_start:pid_end]
                    todSigma = lfi.todSigma
                    todSigma[1]['sigma0'] = sigma0*gain[0]
                    compArray = [lfi.todDtype, todSigma, lfi.huffTod]
                    if(args.no_compress or args.no_compress_tod):
                        compArray = [lfi.todDtype] 
                    comm_tod.add_field(prefix + '/tod', tod, compArray)
                else: #undifferenced data
               
                    diode_list = []
                    name_list = []
 
                    for diode in lfi.diodeTypes[hornType]:
                        ref = undiffFile[str(horn) + diode + '/REF'][pid_start:pid_end]
                        diode_list.append(ref)
                        name_list.append('ref'+diode)
                        sky = undiffFile[str(horn) + diode + '/SKY'][pid_start:pid_end] 
                        diode_list.append(sky)
                        name_list.append('sky'+diode)

                    compArray = [lfi.todDtype, lfi.huffTod]
                    if(args.no_compress or args.no_compress_tod):
                        compArray = [lfi.todDtype]
                    
                    comm_tod.add_matrix(prefix + '/diodes', diode_list, name_list, compArray)
                    

 
        comm_tod.finalize_chunk(pid, loadBalance=outAng)
    comm_tod.finalize_file()


if __name__ == '__main__':
    main()
