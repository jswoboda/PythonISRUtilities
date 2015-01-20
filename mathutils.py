#!/usr/bin/env python
"""
mathutils.py
Created on Thu Jun 12 17:43:43 2014
This module has various math function that are not included in scipy or numpy.
@author: John Swoboda
"""

import numpy as np
import scipy.fftpack as fftsy

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
    L = N+M
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
    yn[:M] =  W**(k**2/2.0)
    yn[M:] = W**((L-(n+M))**2/2.0)
    # perform the circular convolution and multiply by chirp
    gk =fftsy.ifft(fftsy.fft(xn)*fftsy.fft(yn))[:M]
    yk = gk*W**(k**2/2.0)

    return yk