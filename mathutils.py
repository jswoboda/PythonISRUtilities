#!/usr/bin/env python
"""
mathutils.py
Created on Thu Jun 12 17:43:43 2014
This module has various math function that are not included in scipy or numpy.
@author: John Swoboda
"""

import numpy as np
import scipy.fftpack as fftsy
import scipy.special
import matplotlib.pylab as plt
import pdb

def diric(x,M):
    """ This calculates the dirichete sinc function """
    if M % 1 != 0 or M <= 0:
        raise RuntimeError('n must be a strictly positive integer')

    y=np.sin(0.5*x)
    ilog = np.abs(y) < 1e-12
    nilog = np.logical_not(ilog)
    y[nilog]=np.sin((M/2)*x[nilog])/(M*y[nilog])
    y[ilog]=np.sign(np.cos(x[ilog]*((M+1.0)/2.0)))
    return y

def phys2array(az,el):
    """ This takes the physical angles of azimuth and elevation in degrees
    and brings them to the array space."""

    azt = (az)*np.pi/180.0
    elt = 90-el
    xout = elt*np.sin(azt)
    yout = elt*np.cos(azt)
    return (xout,yout)

def chirpz(Xn,A,W,M):
    """ chirpz(Xn,A,W,M)
        by John Swoboda
        This function calculates the chirpz transfrom for the numpy array Xn given the
        complex constants A and W along with the length of the final array M.
        Inputs
        Xn - The signal that the Chirp z transform will be calculated for.
        A - A complex constant used to determine the direction of the integration in
            Z
        W - Another complex constant that will determine the direction of the integration in Z.
        M - The length of the final chirpz transfrom.
        Output
        yk - The M length chirp z tranfrom given Xn and the complex constants.
        """
    N = Xn.shape[0]
    # Make an L length output so the circular convolution does not wrap. Added 1 extra sample to make coding easier.
    L = np.power(2,np.ceil(np.log2(N+M-1)))
    # chirpz numbering
    k = np.arange(M)
    #numbering for time axis
    n = np.arange(N)
    # Make complex arrays
    xn = np.zeros(L) +1j* np.zeros(L)
    yn = np.zeros(L) +1j* np.zeros(L)
    #complex constants raised to power
    An = A**(-n)
    Wn = W**(n**2/2.0)

    xn[:N] = Xn*An*Wn
    # Make the chirp kernal
    yn[:M] =  W**(-k**2/2.0)
    nn = np.arange(L-N+1,L)
    yn[L-N+1:] = W**(-(L-nn)**2/2.0)
    # perform the circular convolution and multiply by chirp
    gk =fftsy.ifft(fftsy.fft(xn)*fftsy.fft(yn))[:M]
    yk = gk*W**(k**2/2.0)

    return yk

def sommerfeldchirpz(func,N,M,dk,exparams=None,a=-1.0*1j,p=1.0,x0=None):

    k = np.arange(N)*dk
    if x0 is None:
        x0 = np.pi/dk
    if exparams is None:
        fk = func(k)
    else:
        fk = func(k,*exparams)

    wk = np.ones(N)
    wk[np.mod(k,2)==0] = 2.0/3.0
    wk[np.mod(k,2)==1] = 4.0/3.0
    wk[0] = 1.5/3.0
    wk[-1] = 1.5/3.0
    A_0 = np.exp(a*dk*x0)
    M2 = M - np.mod(M,2)
    W_0 = np.exp(a*2.0*p*np.pi/M2)
    Xk = chirpz(fk*wk,A_0,W_0,M)
    return Xk

def sommerfelderf(func,N,theta,a,b,exparams=None):


    nvec = np.arange(-N,N+1)

    h = np.log(1.05*np.sqrt(2*N))/N
    kn = 0.5*(b+a)+0.5*(b-a)*scipy.special.erf(np.sinh(nvec*h))

    An = np.cosh(nvec*h)*np.exp(-np.power(np.sinh(nvec*h),2))

    Xk3 = np.zeros_like(theta)*1.0j

    if exparams is None:
        fk = func(kn)
    else:
        fk = func(kn,*exparams)
        #XXX make matrix version
    for omegk, omeg in enumerate(theta):
        fkn = np.exp(-1j*kn*omeg) *fk
        Xk3[omegk] = np.sum(fkn*An)
    return Xk3*h*(b-a)/np.sqrt(np.pi)
