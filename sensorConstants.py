#!/usr/bin/env python
"""
Created on Tue Nov 26 12:42:11 2013

@author: John Swoboda
These are system constants for various sensors
"""
import pdb
import os

import tables
import numpy as np
from scipy.interpolate import griddata
import scipy as sp
from mathutils import diric, angles2xy, jinc, rotcoords
from physConstants import v_C_0
## Parameters for Sensor
#AMISR = {'Name':'AMISR','Pt':2e6,'k':9.4,'G':10**4.3,'lamb':0.6677,'fc':449e6,'fs':50e3,\
#    'taurg':14,'Tsys':120,'BeamWidth':(2,2)}
#AMISR['t_s'] = 1/AMISR['fs']


def getConst(typestr,angles = None):
    """Get the constants associated with a specific radar system. This will fill
    out a dictionary with all of the parameters."""
    dirname, filename = os.path.split(os.path.abspath(__file__))
    if typestr.lower() =='risr':
        arrayfunc = AMISR_Patternadj
        h5filename = os.path.join(dirname,'RISR_PARAMS.h5')

    elif typestr.lower() =='pfisr':
        arrayfunc = AMISR_Patternadj
        h5filename = os.path.join(dirname,'PFISR_PARAMS.h5')
    elif typestr.lower() =='millstone':
        arrayfunc = Millstone_Pattern
        h5filename = os.path.join(dirname,'Millstone_PARAMS.h5')
    elif typestr.lower() =='sondrestrom':
        arrayfunc = Sond_Pattern
        h5filename = os.path.join(dirname,'Sondrestrom_PARAMS.h5')
    h5file = tables.open_file(h5filename)
    kmat = h5file.root.Params.Kmat.read()
    freq = float(h5file.root.Params.Frequency.read())
    P_r = float(h5file.root.Params.Power.read())
    bandwidth = h5file.getNode('/Params/Bandwidth').read()
    ts = h5file.getNode('/Params/Sampletime').read()
    systemp = h5file.getNode('/Params/Systemp').read()
    Ang_off = h5file.root.Params.Angleoffset.read()
    h5file.close()
    Ksens = freq*2*np.pi/v_C_0
    lamb=Ksens/2.0/np.pi
    az = kmat[:,1]
    el = kmat[:,2]
    ksys = kmat[:,3]

    (xin,yin) = angles2xy(az,el)
    points = sp.column_stack((xin,yin))
    if angles is not None:
        (xvec,yvec) = angles2xy(angles[:,0],angles[:,1])
        ksysout = griddata(points, ksys, (xvec, yvec), method='nearest')
    else:
        ksysout = None

    #'G':10**4.3, This doesn't get used anymore it seems
    sensdict = {'Name':typestr,'Pt':P_r,'k':Ksens,'lamb':lamb,'fc':freq,'fs':1/ts,\
    'taurg':14,'Tsys':systemp,'BeamWidth':(2,2),'Ksys':ksysout,'BandWidth':bandwidth,\
    'Angleoffset':Ang_off,'ArrayFunc':arrayfunc}
    sensdict['t_s'] = ts
    return sensdict



def phys2array(az,el):
    """ This takes the physical angles of azimuth and elevation in degrees
    and brings them to the array space."""

    azt = (az)*np.pi/180.0
    elt = 90-el
    xout = elt*np.sin(azt)
    yout = elt*np.cos(azt)
    return (xout,yout)

def AMISR_Patternadj(Az,El,Az0,El0,Angleoffset):
    """
    AMISR_Patternadj
    by John Swoboda
    This function will call """
    d2r= np.pi/180.0
    Azs = np.mod(Az-Angleoffset[0],360.0)
    Az0s = np.mod(Az-Angleoffset[0],360.0)

    Els = El+Angleoffset[1]
    elg90 = Els>90.0
    Els[elg90] = 180.0-Els[elg90]
    Azs[elg90] = np.mod(Azs[elg90]+180.0,360.0)

    El0s = El0+Angleoffset[1]
    if El0s>90.0:
        El0s= 180.0-El0s
        Az0s = np.mod(Az0s+180.0,360.0)


    Elr = (90.0-Els)*d2r
    El0r = (90-El0s)*d2r
    Azr = Azs*d2r
    Az0r = Az0s*d2r

    return AMISR_Pattern(Azr,Elr,Az0r,El0r)
def Sond_Pattern(Az,El,Az0,El0,Angleoffset):



    d2r= np.pi/180.0
    r = 30.
    lamb = v_C_0/1.2e9

    Azadj,Eladj = rotcoords(Az,El,-Az0,El0-90.)
    Elr = (90.0-Eladj)*d2r
    return Circ_Ant_Pattern(Elr,r,lamb)

def Millstone_Pattern(Az,El,Az0,El0,Angleoffset):
    d2r= np.pi/180.0
    r = 34.
    lamb = v_C_0/4.4e8
    Azadj,Eladj = rotcoords(Az,El,-Az0,El0-90.)
    Elr = (90.0-Eladj)*d2r
    return Circ_Ant_Pattern(Elr,r,lamb)

def Circ_Ant_Pattern(EL,r,lamb):
    """
    Circ_Ant_Pattern
    by John Swoboda
    This function will create an idealized antenna pattern for a circular antenna
    array. The pattern is not normalized.
    The antenna is assumed to made of a grid of ideal cross dipole
    elements. In the array every other column is shifted by 1/2 dy. The
    parameters are taken from the AMISR spec and the method for calculating
    the field is derived from a report by Adam R. Wichman.
    The inputs for the az and el coordinates can be either an array or
    scalar. If both are arrays they must be the same shape.
    ###########################################################################
    Inputs

    EL - An array or scalar holding the elevation coordinates in radians.
       Also vertical is at zero radians.
    r - Radius of the antenna in meters
    lamb - wavelength of radiation in meters
    ###########################################################################
    Outputs
    Patout - The normalized radiation density.
    ###########################################################################"""

    Patout = (2.*r/lamb)**2* np.abs(jinc((2.*r/lamb)*np.sin(EL)))
    normfactor = (2.*r/lamb)**2* jinc(0.)
    return Patout/normfactor



def AMISR_Pattern(AZ,EL,Az0,El0):
    """
    AMISR_Pattern
    by John Swoboda
    This function will create an idealized antenna pattern for the AMISR
    array. The pattern is not normalized.
    The antenna is assumed to made of a grid of ideal cross dipole
    elements. In the array every other column is shifted by 1/2 dy. The
    parameters are taken from the AMISR spec and the method for calculating
    the field is derived from a report by Adam R. Wichman.
    The inputs for the az and el coordinates can be either an array or
    scalar. If both are arrays they must be the same shape.
    ###########################################################################
    Inputs
    Az - An array or scalar holding the azimuth coordinates in radians.
    EL - An array or scalar holding the elevation coordinates in radians.
       Also vertical is at zero radians.
    Az0 - A scalar that determines the azimuth pointing angle of the antenna.
    El0 - A scalar that determines the elevation pointing angle of the
    antenna.
    ###########################################################################
    Outputs
    Patout - The normalized radiation density.
    ###########################################################################"""
    f0=440e6 # frequency of AMISR in Hz
    lam0=v_C_0/f0 # wavelength in m
    k0=2*np.pi/lam0 # wavenumber in rad/m

    dx=0.4343 # x spacing[m]
    dy=0.4958 # y spacing[m]
    # element pattern from an ideal cross dipole array.
    elementpower=(1.0/2.0)*(1.0+(np.cos(EL))**2)
#    pdb.set_trace()
    m=8.0;# number of pannels in the x direction
    mtot = 8.0*m;# number of elements times panels in x direction

    n = 16.0;# number of pannels in the y direction
    ntot = n*4.0;# number of elements times panels in y direction
    # relative phase between the x elements
    phix = k0*dx*(np.sin(EL)*np.cos(AZ)-np.sin(El0)*np.cos(Az0))
    # relative phase between the y elements
    phiy = k0*dy*(np.sin(EL)*np.sin(AZ)-np.sin(El0)*np.sin(Az0))

    AF = (1.0/2.0)*(1.0+np.exp(1j*((1/2)*phiy+phix)))*diric(2.0*phix,mtot/2.0)*diric(phiy,ntot);
    arrayfac = abs(AF)**2;
    Patout = elementpower*arrayfac
    return Patout