import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
import h5py

import huffman

from multiprocessing import Pool

from time import time
from tqdm import tqdm


from scipy import sparse


def cg_solve(A, b, Minv, imax=1000, eps=1e-6):
    x = np.zeros_like(b)
    i = 0
    r = b - A.dot(x)
    d = Minv.dot(r)
    delta_new = r.dot(d)
    delta_0 = r.dot(d)
    while ((i < imax) & (delta_new > eps**2*delta_0)):
        q = A.dot(d)
        alpha = delta_new/d.dot(q)
        x = x + alpha*d
        if (i % 50 == 0):
            r = b - A.dot(x)
        else:
            r = r - alpha*q
        s = Minv.dot(r)
        delta_old = np.copy(delta_new)
        delta_new = r.dot(s)
        beta = delta_new/delta_old
        d = s + beta*d
        i += 1
    return x

def cg_solve_map(A, b, Minv, imax=1000, eps=1e-6):
    '''
    Try to rewrite this so that A is an operator, not a matrix that needs to be
    held in memory.
    '''
    x = np.zeros_like(b)
    i = 0
    r = b - A.dot(x)
    d = Minv.dot(r)
    delta_new = r.dot(d)
    delta_0 = r.dot(d)
    while ((i < imax) & (delta_new > eps**2*delta_0)):
        q = A.dot(d)
        alpha = delta_new/d.dot(q)
        x = x + alpha*d
        if (i % 50 == 0):
            r = b - A.dot(x)
        else:
            r = r - alpha*q
        s = Minv.dot(r)
        delta_old = np.copy(delta_new)
        delta_new = r.dot(s)
        beta = delta_new/delta_old
        d = s + beta*d
        i += 1
    return x


def cg_test():
    A = np.array([[3,2],
                  [2,6]])
    b = np.array([2,-8])
    Minv = np.eye(2)

    x = cg_solve(A, b, Minv)
    assert np.allclose(A.dot(x), b), 'CG solution is not close enough'
    return

def get_data(fname, nside=256):
    npix = hp.nside2npix(nside)
    M = np.zeros(npix)
    b = np.zeros(npix)
    labels = ['K113', 'K114', 'K123', 'K124']
    f= h5py.File(fname, 'r')
    obsid = str(list(f.keys())[0])
    TOD0 = np.array(f[obsid + '/' + labels[0] + '/tod'])
    if len(TOD0) != 675000:
        print(f'{fname} has wrong length')
        return None
    
    
    DAs = [[], [], [], []]
    pixAs = []
    pixBs = []
    sigmas = []
    for num, label in enumerate(labels):
        TODs = np.array(f[obsid + '/' + label + '/tod'])
        gain, sigma0 = f[obsid + '/' + label + '/scalars'][0:2]
        DAs[num] = DAs[num] + TODs.tolist()
        sigmas.append(sigma0)
        if label == 'K113':
            pixA = np.array(f[obsid + '/' + label + '/pixA']).tolist()
            pixB = np.array(f[obsid + '/' + label + '/pixB']).tolist()

    DAs = np.array(DAs)
    sigma0 = sum(np.array(sigmas)**2)**0.5

    
    
    d1 = 0.5*(DAs[0] + DAs[1])
    d2 = 0.5*(DAs[2] + DAs[3])
    
    d = 0.5*(d1 + d2) # = i_A - i_B
    p = 0.5*(d1 - d2) # = q_A*cos(2*g_A) + u_A*sin(2*g_A) - q_B*cos(2*g_B) - u_B*sin(2*g_B)



    for t in range(len(d)):
        b[pixA[t]] += d[t]/sigma0
        b[pixB[t]] -= d[t]/sigma0

        M[pixA[t]] += sigma0**-2
        M[pixB[t]] += sigma0**-2

    return M, b, pixA, pixB, sigma0


def get_cg_K(nside=256):
    npix = hp.nside2npix(nside)
    b = np.zeros(hp.nside2npix(nside))
    M_diag = np.zeros(npix)

    from glob import glob
    fnames = glob('/mn/stornext/d16/cmbco/bp/wmap/data/wmap_K1_*v6.h5')
    fnames.sort()
    fnames = fnames[:1000]
    #fnames = [fnames[0]]
    pool = Pool(processes=120)
    x = [pool.apply_async(get_data, args=[fname]) for fname in fnames]
    pixA = []
    pixB = []
    sigma0s = []
    for res in x:
        r = res.get()
        if r is not None:
            M_i, b_i, pixA_i, pixB_i, sigma_i = r
            M_diag += M_i
            b += b_i
            pixA += pixA_i
            pixB += pixB_i
            sigma0s += [sigma_i]
    print('Loaded data')
    sigma0 = np.mean(sigma0s)


    times = np.arange(len(pixA))
    print('Creating the sparse matrices')
    t0 = time()
    P_A = sparse.csr_matrix((np.ones_like(times), (times, pixA)))
    P_B = sparse.csr_matrix((np.ones_like(times), (times, pixB)))
    print(f'sparse matrix construction takes {time()-t0} seconds')
    P = P_A - P_B

    print('Constructing the CG matrix A')
    t0 = time()
    A = P.T.dot(P)/sigma0**2
    print(f'Inner product takes {time()-t0} seconds')






    dts = []
    i = 0
    x = np.zeros_like(b)
    d = np.zeros_like(b)
    q = np.zeros_like(b)
    s = np.zeros_like(b)
    r = b
    p = (M_diag != 0)
    d[p] = r[p]/M_diag[p]
    delta_new = r.dot(d)
    delta_0 = np.copy(delta_new)
    i_max = npix
    i_max = 1000
    eps = 1e-7
    while ((i < i_max) & (delta_new > eps**2*delta_0)):
        #print(delta_new)
        t0 = time()
        q = A.dot(d)
        alpha = delta_new/d.dot(q)
        x = x + alpha*d
        if i % 50 == 0:
            #print('Divisible by 50')
            r = b - A.dot(x)
        else:
            r = r - alpha*q
        if i % 10 == 0:
            print(np.round(i/i_max,3),\
                    np.round(1/np.log10(delta_new/(delta_0*eps**2)),3),\
                    int(delta_new))
        s[p] = r[p]/M_diag[p]
        delta_old = np.copy(delta_new)
        delta_new = r.dot(s)
        beta = delta_new/delta_old
        d = s + beta*d
        i += 1
        dts.append(time()-t0)
    hp.write_map('cg_v3.fits', x, overwrite=True)

    print(f'Done with {i} iterations, delta is {delta_new}')
    print(f"Each iteration is {np.mean(dts)}\pm{np.std(dts)}")

    #hp.mollview(b, cmap='coolwarm', norm='hist', title='Noise-weighted average')
    #hp.mollview(M_diag, norm='hist', title='Preconditioner')
    #hp.mollview(x, norm='hist', title='Solution')
    #hp.mollview(x, min=-250, max=250, title='Solution')


    #plt.show()

    

    return


def check_hdf5(nside=256):
    # Take official W-band K-band, scan it with the same pointing matrix, divide
    # by gain, subtract from timestream, check the gain and pointing solution
    # directly. Should just be white noise.
    npix = hp.nside2npix(nside)
    b = np.zeros(hp.nside2npix(nside))
    M_diag = np.zeros(npix)

    from glob import glob
    fnames = glob('/mn/stornext/d16/cmbco/bp/wmap/data/wmap_K1_*v6.h5')
    fnames.sort()
    fname = fnames[500]
    f= h5py.File(fname, 'r')
    obsid = str(list(f.keys())[0])

    DAs = [[], [], [], []]
    pixAs = []
    pixBs = []
    sigmas = []
    labels = ['K113', 'K114', 'K123', 'K124']
    for num, label in enumerate(labels):
        TODs = np.array(f[obsid + '/' + label + '/tod'])
        gain, sigma0 = f[obsid + '/' + label + '/scalars'][0:2]
        DAs[num] = DAs[num] + TODs.tolist()
        sigmas.append(sigma0)
        if label == 'K113':
            pixA = np.array(f[obsid + '/' + label + '/pixA']).tolist()
            pixB = np.array(f[obsid + '/' + label + '/pixB']).tolist()

    DAs = np.array(DAs)
    sigma0 = sum(np.array(sigmas)**2)**0.5
    
    d1 = 0.5*(DAs[0] + DAs[1])
    d2 = 0.5*(DAs[2] + DAs[3])
    
    d = 0.5*(d1 + d2) # = i_A - i_B
    p = 0.5*(d1 - d2) # = q_A*cos(2*g_A) + u_A*sin(2*g_A) - q_B*cos(2*g_B) - u_B*sin(2*g_B)

    sol = hp.read_map('data/wmap_imap_r9_9yr_K1_v5.fits')
    sol = hp.ud_grade(sol, nside)
    #hp.mollview(sol, min=-250, max=250)

    d_sol = np.zeros(len(pixA))
    for t in range(len(pixA)):
        d_sol[t] = gain*(sol[pixA[t]] - sol[pixB[t]])
    plt.figure()
    plt.plot(d_sol)
    plt.plot(d)
    plt.show()

    return


if __name__ == '__main__':
    #cg_test()
    get_cg_K()
    #check_hdf5()


