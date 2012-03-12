#!usr/bin/env python
import numpy
from scipy import special
from numpy import *
import math
#-----------------------------------------------
class orth_polynomials:
    def __init__(self, n_root, n_dim):
        self.n_root = n_root
        self.n_dim = n_dim
        self.norm = []
        self.roots = None
        self.weights = None
        self.pols = []

class orth_legendre(orth_polynomials):
    pass

def set_legendre(n_root, n_dim):
    legP = orth_legendre(n_root, n_dim)
    # set normalization constant
    legP.norm.append( 2.0 )  #for i = 0
    for i in range(1, n_dim):
        legP.norm.append( 2.0 / (2.0*i + 1.0) )
    # set roots and weights
    legP.roots, legP.weights = special.orthogonal.p_roots(n_root)
    legP.roots = legP.roots.tolist()
    legP.weights = legP.weights.tolist() # convert into lists
    legP.roots.insert(0, 0.0)
    legP.weights.insert(0, 0.0)
    # generate legendre polynomial object
    for i in range(0, n_dim, 1):
        legP.pols.append( special.legendre(i) )
    return legP
#-----------------------------------------------
class Tmatrix:
    def __init__(self):
        self.val = None

def makeTmat(n_dvr, orth_pol):
    Tmat = Tmatrix()
    Tmat.val = numpy.zeros( (n_dvr+2)*(n_dvr+2) )
    Tmat.val.shape = (n_dvr+2, n_dvr+2) # shape into (n_dvr+1)x(n_dvr+1) matrix
    # check type
    if isinstance(orth_pol, orth_legendre):
        for n in range(1, n_dvr+1, 1):
            for i in range(1, n_dvr+1, 1):
                Tmat.val[n,i] = \
                    math.sqrt( orth_pol.weights[i] \
                                   / orth_pol.norm[n-1] ) \
                                   * orth_pol.pols[n-1](orth_pol.roots[i])
    else:
        print 'invalid class for 2nd argument'
    return Tmat
                        
def checkTmat(n_dvr, Tmat):  # n_dvr: number of dvr-basis
    out = numpy.zeros( (n_dvr+2)*(n_dvr+2) )
    out.shape = (n_dvr+2, n_dvr+2)
    if isinstance(Tmat, Tmatrix):
        TmatT = Tmatrix()
        TmatT.val = Tmat.val.transpose()
        for i in range(0, n_dvr+2, 1):
            for j in range(0, n_dvr+2, 1):
                for k in range(0, n_dvr+2, 1):
                    out[i,j] = out[i,j] + Tmat.val[i,k] * TmatT.val[k,j]
    print out

#-----------------------------------------------
class DVR_Leg:
    def __init__(self, n_dvr):
        self.n_dvr = n_dvr
        self.val = [0.0] # 0th argument is 0(not used).
        self.weight_f = 1.0  # weight function is always 1 for Legendre
        
# DVR basis with one argument
def DVR_Leg1arg(n_dvr, Tmat, orth_leg, x):
    if ( type(float(x)) != float ):
        print 'Bad DVR-argument in DVR_Leg1arg\n'
        sys.exit(1)

    dvrBasis = DVR_Leg(n_dvr)
    if isinstance(orth_leg, orth_legendre):
        midfunc = numpy.zeros( n_dvr + 1 )
        for n in range(1, n_dvr+1, 1):
            midfunc[n] = sqrt( dvrBasis.weight_f / orth_leg.norm[n-1] ) \
                * orth_leg.pols[n-1](x)
            sums = 0.0  # initialize
            for m in range(1, n_dvr+1, 1):
                sums = midfunc[m] * Tmat.val[m,n]
            dvrBasis.val.append( sums )
            sums = 0.0  # zeroclear
    return dvrBasis
#-----------------------------------------------
