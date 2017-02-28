#!/usr/bin/env python
"""
.. module:: mathutils
    :platform: Unix, Windows, Mac
    :synopsis: This module has various math function that are not included in scipy or numpy.

.. moduleauthor:: John Swoboda <swoboj@bu.edu>

"""

import numpy as np
import scipy.fftpack as fftsy
import scipy.special
#import matplotlib.pylab as plt
#import pdb
try:
    from numba import jit
except ImportError:
    print('Numba not installed, falling back to Numpy')

def diric(x,n):
    """ Based on octave-signal v. 1.10 diric. Jitted for up to 2x speedup in this case, and demo of Numba 0.17
    
    :Authors: Michael Hirsch
    """
    n = int(n)
    if n < 1:
        raise RuntimeError('n must be a strictly positive integer')
    return _diric(np.asanyarray(x),n)

@jit
def _diric(x,n):
    """ Based on octave-signal v. 1.10 diric. Jitted for up to 2x speedup in this case, and demo of Numba 0.17
    
    :Authors:Michael Hirsch
    """
    y = np.divide(np.sin(n*x/2),(n*np.sin(x/2)))
    #edge case
    badx = np.isnan(y)
    y[badx] = (-1)**((n-1)*x[badx]/(2*np.pi))
    return y
@jit
def jinc(t):
    """ This will output a jinc function.
    
    Args:
        t (:obj:`numpy array`): Time array in seconds.
        
    Returns:
        outdata (:obj:`numpy array`): Jinc(t)
    """
    t = np.asanyarray(t)

    outdata = np.ones_like(t)
    outdata[t==0.] = 0.5
    outdata[t!=0.0] =scipy.special.jn(1,np.pi*t[t!=0.0])/(2*t[t!=0.0])

    return outdata

def angles2xy(az,el):
    """Creates x and y coordinates from az and el arrays.
    
    This function will take a numpy array of az and elevation angles in degrees and
    create the x, y locations projected on to the z=0 plane. It is assumed that the
    elevation angle is mesured from the z=0 plane. The az angle is from x=0 and goes
    clockwise.
    
    Args:
        az (:obj:`numpy array`): Azimuth angles in degrees.
        el (:obj:`numpy array`): Elevation angles in degrees.
    
    Returns:
        (xout (:obj:`numpy array`), yout (:obj:`numpy array`)): x and y coordinates.
   """

    azt = (az)*np.pi/180.0
    elt = 90-el
    xout = elt*np.sin(azt)
    yout = elt*np.cos(azt)
    return (xout,yout)

def array2cart(Az,El):
    """ This function will turn azimuth and elevation angles to X, Y and Z coordinates
        assuming a unit sphere.
        
    Args:
        Az (:obj:`numpy array`): Azimuth angles in degrees.        
        El (:obj:`numpy array`): Elevation angles in degrees.
    
    Returns:
        (X (:obj:`numpy array`),  Y (:obj:`numpy array`), Z (:obj:`numpy array`)): x, y and z coordinates.        
    """
    d2r = np.pi/180.
    X = np.cos(Az*d2r)*np.cos(El*d2r)
    Y = np.sin(Az*d2r)*np.cos(El*d2r)
    Z = np.sin(El*d2r)
    return (X,Y,Z)

def cart2array(X,Y,Z):
    """ This function will turn the X, Y and Z coordinate to azimuth and elevation angles 
        assuming a unit sphere.
        
    Args:
        X (:obj:`numpy array`): x coordinates.
        Y (:obj:`numpy array`): y coordinates.
        Z (:obj:`numpy array`): z coordinates.
    
    Returns:
        (Az (:obj:`numpy array`), El (:obj:`numpy array`)): Azimuth and elevation angles in degrees.
    """
    r2d = 180./np.pi
    Az = np.arctan2(Y,X)*r2d
    El = np.arcsin(Z)*r2d
    return (Az,El)
def rotmatrix(Az_0,El_0):
    """ Makes a rotation matrix.
    
    This creates a rotation matrix for the rotcoords function. First rotate 
    about the z axis and then rotate about the new y axis
    http://www.agi.com/resources/help/online/stk/11.0/index.htm#comm/CommRadar03-03.htm
    
    Args:
        Az_0 (float): The azimuth rotation angle in degrees. 
        El_0 (float): The elevation rotation angle in degrees.
        
    Return:
        rotmat (:obj:`numpy array`): A 3x3 rotation matrix.
    """
    d2r = np.pi/180.
    Az_0= Az_0*d2r
    El_0 = El_0*d2r
    R_Az = np.array([[np.cos(Az_0),-np.sin(Az_0),0.],[np.sin(Az_0),np.cos(Az_0),0.],[0.,0.,1.]])
    R_El = np.array([[np.cos(El_0),0.,np.sin(El_0)],[0.,1.,0.],[-np.sin(El_0),0.,np.cos(El_0)]])
    return R_El.dot(R_Az)

def rotcoords(Az,El,Az_0,El_0):
    """ Applies rotation matrix to Az and El arrays.
    
    This function will rotate the Az and Elevation cordinates given offset
    angles. This will use a rotation matrix after the angles have been 
    changed to Cartisian coordinates assuming a unit sphere.
    
    Args:
        Az (:obj:`numpy array`): Azimuth angles in degrees.
        El (:obj:`numpy array`): Elevation angles in degrees.
        Az_0 (float): The azimuth rotation angle in degrees. 
        El_0 (float): The elevation rotation angle in degrees.
        
    Returns:
        (Az_out (:obj:`numpy array`), El_out (:obj:`numpy array`)): Rotated azimuth and elevation angles in degrees.
   """
    cartcords = array2cart(Az,El)
    cartmat = np.column_stack(cartcords).transpose()
    rotmat = rotmatrix(Az_0,El_0)

    rotcart = np.dot(rotmat,cartmat)
    return cart2array(rotcart[0],rotcart[1],rotcart[2])
def chirpz(Xn,A,W,M):
    """ Chirpz calculation for a single array.
    
    This function calculates the chirpz transfrom for the numpy array Xn given the
    complex constants A and W along with the length of the final array M.
    
    Args:
        Xn (:obj:`numpy array`): The signal that the Chirp z transform will be calculated for.
        A (:obj:`complex`): A complex constant used to determine the direction of the integration in Z.
        W (:obj:`complex`): Another complex constant that will determine the direction of the integration in Z.
        M (int): The length of the final chirpz transfrom.
    
    Returns: 
        yk (:obj:`numpy array`): The M length chirp z tranfrom given Xn and the complex constants.
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

def sommerfeldchirpz(func,N,M,dk,Lmax=1,errF=0.1,a=-1.0j,p=1.0,x0=None,exparams=()):
    """ Numerically integrate Sommerfeld like integral using chirpz.
    
    This function will numerically integrate a Sommerfeld like integral, int(exp(awk)f(k),t=0..inf)
    using at most Lmax number of N length chirpz transforms to make an M length
    array. If the normalized difference between the previous estimate of the output Xk is
    less then the parameter errF then the loop stops and a flag that represents convergence is set to true.
    A number of repeats are also output as well. This technique is based off the article
    Li et. al Adaptive evaluation of the Sommerfeld-type integral using the chirp z-transform, 1991.
    This function also uses a modified Simpsons rule for integration found in Milla PhD Thesis (2010).

    Args:
        func (func): A function that is used to create f(k).
        N (int): The length of the chirp z transform used.
        M (int): The length of the output array.
        dk (float): The sample period of k.
        Lmax (:obj:`int`, optional): default 1, The maximum number of repeats of the integration before the loop finishes.
        errF (:obj:`float`, optional): default .1, The threshold of the normalized difference between the new iteration and the old to stop to stop the iteration.
        a (:obj:`complex`, optional): default -1.0*1j, A complex number that determines the trejectory of the integration on the z plane.
        p (:obj:`float`, optional): default 1.0, A real number that helps to define the spacing on the omega plane.
        x0(:obj:`complex`, optional): default p*np.pi/dk, The starting point on the omega plane.
        exparams(:obj:`tuple`, optional): default (), Any extra params other then k to create f(k).
        
    Returns:
        (Xk (:obj:`numpy array`),flag_c (bool),irep (int)): The integrated data that is of length M.
        A convergence flag. The number of repetitions until convergence.
    """

    k = np.arange(N)*dk
    if x0 is None:
        x0 = p*np.pi/dk

    #set up simpson rule
    wk = np.ones(N)
    wk[np.mod(k,2)==0] = 2.0/3.0
    wk[np.mod(k,2)==1] = 4.0/3.0
    wk[0] = 1.5/3.0
    wk[-1] = 1.5/3.0
    # Make A_0 and W_0
    A_0 = np.exp(a*dk*x0)
    M2 = M - np.mod(M,2)
    W_0 = np.exp(a*2.0*p*np.pi/M2)
    # Make frequency array for phase offset
    freqm = np.arange(-np.ceil((M-1)/2.0),np.floor((M-1)/2.0)+1)
    Xk = np.zeros(M)*1j
    flag_c = False


    for irep in range(Lmax):

        fk = func(k+N*dk*irep,*exparams)
        Xkold = Xk
        Xk = chirpz(fk*wk,A_0,W_0,M)*np.power(W_0,N*dk*irep*freqm)+Xk

        Xkdiff = np.sqrt(np.sum(np.power(np.abs(Xk-Xkold),2.0)))
        Xkpow = np.sqrt(np.sum(np.power(np.abs(Xk),2.0)))
        outrep = irep+1
        # check for convergence
        if Xkdiff/Xkpow<errF:
            flag_c = True

            break

    return (Xk,flag_c,outrep)


def sommerfelderfrep(func,N,omega,b1,Lmax=1,errF=0.1,exparams=()):
    """  Numerically integrate Sommerfeld like integral using erf transform function loop.
    
    This function will numerically integrate a Sommerfeld like integral, int(exp(-jwk)f(k),k=0..inf)
    using the ERF transform and 2N+1 samples and at most Lmax loops. If the normalized difference
    between the previous estimate of the output Xk is less then the parameter errF then the loop stops
    and a flag that represents convergence is set to true. A number of loops is also output as well.
    This function uses sommerfelderf to do the integration
    
    Args:
        func (func): A function that is used to create f(k).
        N (int): The integration uses 2N+1 samples to do the integration.
        omega (float): The Frequency array in radians per second.
        b1 (int): The inital bounds of the first try and then step size for each subsiquent integral.
        Lmax (:obj:`int`, optional): default 1, The maximum number of repeats of the integration before the loop finishes.
        errF (:obj:`float`, optional): default .1, The threshold of the normalized difference between the new
            iteration and the old to stop to stop the iteration.
        exparams(:obj:`tuple`, optional): default (), Any extra params other then k to create f(k).
        
    Returns:
        (Xk (:obj:`numpy array`),flag_c (bool),irep (int)): The integrated data that is of length M.
        A convergence flag. The number of repetitions until convergence.
    """
    Xk =np.zeros_like(omega)*(1+1j)
    flag_c=False
    for irep in range(Lmax):

        Xktemp = sommerfelderf(func,N,omega,b1*irep,b1*(irep+1),exparams)
        Xkdiff = Xktemp.real**2+Xktemp.imag**2
        Xk = Xk+Xktemp
        Xkpow = Xk.real**2 +Xk.imag**2
#        Xkdiff = np.sqrt(np.sum(np.power(np.abs(Xktemp),2.0)))
#        Xkpow = np.sqrt(np.sum(np.power(np.abs(Xk+Xktemp),2.0)))
#        Xk = Xk+Xktemp
        outrep = irep+1
        # check for convergence
        if np.sum(Xkdiff/Xkpow)<errF:
            flag_c = True
            outrep = irep
            break
    return (Xk,flag_c,outrep)
def sommerfelderf(func,N,omega,a,b,exparams=()):
    """ Integrate somerfeld integral using ERF transform for single portion.

    This function will numerically integrate a Sommerfeld like integral, int(exp(-jwk)f(k),k=a..b)
    using the ERF transform and 2N+1 samples. This technique is from the paper B. L. Ooi 2007.
    
    Args:
        func (func): A function that is used to create f(k).
        N (int): The integration uses 2N+1 samples to do the integration.
        omega: The Frequency array in radians per second.
        a (float): Lower bound of the integral.
        b (float): Upper bound of teh integral.
        exparams (:obj:`tuple`, optional): default (), Any extra params other then k to create f(k).
        
    Returns:
        Xk (:obj:`numpy array`): The integrated data that is the same length as omega.
    """

    nvec = np.arange(-N,N+1)

    h = np.log(1.05*np.sqrt(2*N))/N
    kn = 0.5*(b+a)+0.5*(b-a)*scipy.special.erf(np.sinh(nvec*h))

    An = np.cosh(nvec*h)*np.exp(-np.power(np.sinh(nvec*h),2))

    fk = func(kn,*exparams)
    kmat = np.tile(kn[:,np.newaxis],(1,len(omega)))
    omegamat = np.tile(omega[np.newaxis,:],(len(kn),1))
    Xk3 = np.dot(An*fk,np.exp(-1j*kmat*omegamat));

    return Xk3*h*(b-a)/np.sqrt(np.pi)
