#!/usr/bin/env python
"""
Created on Sat Mar 15 10:21:20 2014

@author: John Swoboda

This script is here to create the h5 files for the ISR systems
"""

import os
import numpy as np
from isrutilities.mathutils import angles2xy
from scipy.interpolate import griddata

import tables
import glob


if __name__ == "__main__":

     # Get only .dt0 files
    flist1 = glob.glob(os.path.expanduser('~')+'/Documents/ISRData/20140224.002/'+'*.dt0.h5')
    flist2 = glob.glob(os.path.expanduser('~')+'/Documents/ISRData/20111215.022/'+'*.dt0.h5')
    #%% Read RISR Params
    curfile = flist1[0]
    filepath, fext =os.path.splitext(curfile)
    # open file
    with tables.open_file(curfile) as f:
        # get the tx power for RISR
        txpowR = f.get_node('/Tx/Power').read()[0,0]
        # get freq 1 for RISR
        txfreqR = f.get_Node('/Tx/Frequency').read()[0,0]
        # get cal temp
        caltempR=120
    #    caltempR = f.getNode('/Rx/CalTemp').read()
        # get bandwidth
        bandwidthR = f.get_node('/Rx/Bandwidth').read()
        # sample time
        sampletimeR = f.get_node('/Rx/SampleTime').read()
        # angle offset
        angoffR = np.array([26.0,45.0])

        beamcodemapR = f.get_node('/Setup/BeamcodeMap').read()
    beamcodemapR = np.array(beamcodemapR)

    posel = np.where(beamcodemapR[:,2]>0)
    beamcodemapR = beamcodemapR[posel]
    #%% Write RISR parameters
    with tables.open_file('RISR_PARAMS.h5',mode='w',title='RISR Parameters') as f:
        fgroup = f.create_group('/','Params','Parameters')

        f.create_array(fgroup, 'Frequency', txfreqR)
        f.create_array(fgroup, 'Power', txpowR)
        f.create_array(fgroup, 'Kmat', beamcodemapR)
        f.create_array(fgroup,'Systemp',caltempR)
        f.create_array(fgroup,'Bandwidth',bandwidthR)
        f.create_array(fgroup,'Sampletime',sampletimeR)
        f.create_array(fgroup,'Angleoffset',angoffR)


    #%% Set up PFISR parameters
    curfile = flist2[0]
    filepath, fext =os.path.splitext(curfile)
    # open file
    h5file = tables.open_file(curfile)
    # get the tx power for RISR
    txpowP = h5file.getNode('/Tx/Power').read()[0,0]
    # get freq 1 for RISR
    txfreqP = h5file.getNode('/Tx/Frequency').read()[0,0]
    # get cal temp
    #hard code in because PFISR example is 3x to high
    caltempP = 120
#    caltempP = h5file.getNode('/Rx/CalTemp').read()
    # get bandwidth
    bandwidthP = h5file.getNode('/Rx/Bandwidth').read()
    # sample time
    sampletimeP = h5file.getNode('/Rx/SampleTime').read()
    # angle offset
    angoffP = np.array([15.0,16.0])

    beamcodemapP = h5file.getNode('/Setup/BeamcodeMap').read()
    beamcodemapP = np.array(beamcodemapP)
    h5file.close()
    posel = np.where(beamcodemapP[:,2]>0)
    beamcodemapP = beamcodemapP[posel]
    #%% Write PFISR parameters
    h5fileout = tables.open_file('PFISR_PARAMS.h5',mode='w',title='PFISR Parameters')
    fgroup = h5fileout.create_group('/','Params','Parameters')

    h5fileout.create_array(fgroup, 'Frequency', txfreqP)
    h5fileout.create_array(fgroup, 'Power', txpowP)
    h5fileout.create_array(fgroup, 'Kmat', beamcodemapP)
    h5fileout.create_array(fgroup,'Systemp',caltempP)
    h5fileout.create_array(fgroup,'Bandwidth',bandwidthP)
    h5fileout.create_array(fgroup,'Sampletime',sampletimeP)
    h5fileout.create_array(fgroup,'Angleoffset',angoffP)
    h5fileout.close()


    #%% Set up Sondrestrom parameters

    # get the tx power for Sondrestrom
    txpowP = 3.5e6
    # get freq 1 for Sondrestrom
    txfreqP = 1.290e9
    # get cal temp
    #hard code in because Sondrestrom example is 3x to high
    caltempP = 85.
#    caltempP = h5file.getNode('/Rx/CalTemp').read()
    # get bandwidth
    bandwidthP = 100e3
    # sample time
    sampletimeP = 1./100e3
    # angle offset
    angoffP = np.array([0.0,0.0])

    azvec = np.linspace(0.,360.,361)
    elvec = np.linspace(25.,89,65)
    (azmat,elmat) = np.meshgrid(azvec,elvec)

    (az,el,ksys)=beamcodemapR.transpose()[1:]
    (xin,yin) = angles2xy(az,el)
    points = np.column_stack((xin,yin))
    (xvec,yvec) = angles2xy(azmat.flatten(),elmat.flatten())
    ksysout = griddata(points, ksys, (xvec, yvec), method='nearest')
    # size ksys to deal with the different wavelength and antenna gain from Sondrestrom
    # makes the returns go down by about half
    ksysout = ksysout*np.power(10.,6./10.)*(.2323/.6677)**2
    beamcodemapP = np.column_stack((np.arange(azmat.size),azmat.flatten(),elmat.flatten(),ksysout))


    #%% Write Sondrestrom parameters
    h5fileout = tables.open_file('Sondrestrom_PARAMS.h5',mode='w',title='Sondrestrom Parameters')
    fgroup = h5fileout.create_group('/','Params','Parameters')

    h5fileout.create_array(fgroup, 'Frequency', txfreqP)
    h5fileout.create_array(fgroup, 'Power', txpowP)
    h5fileout.create_array(fgroup, 'Kmat', beamcodemapP)
    h5fileout.create_array(fgroup,'Systemp',caltempP)
    h5fileout.create_array(fgroup,'Bandwidth',bandwidthP)
    h5fileout.create_array(fgroup,'Sampletime',sampletimeP)
    h5fileout.create_array(fgroup,'Angleoffset',angoffP)
    h5fileout.close()

    #%% Set up Millstone parameters

    # get the tx power for Sondrestrom
    txpowP = 2.5e6
    # get freq 1 for Sondrestrom
    txfreqP = 4.4e8
    # get cal temp
    #hard code in because Sondrestrom example is 3x to high
    caltempP = 120.
#    caltempP = h5file.getNode('/Rx/CalTemp').read()
    # get bandwidth
    bandwidthP = 50e3
    # sample time
    sampletimeP = 1./50e3
    # angle offset
    angoffP = np.array([0.0,0.0])

    azvec = np.linspace(0.,360.,361)
    elvec = np.linspace(25.,89,65)
    (azmat,elmat) = np.meshgrid(azvec,elvec)

    (az,el,ksys)=beamcodemapR.transpose()[1:]
    (xin,yin) = angles2xy(az,el)
    points = np.column_stack((xin,yin))
    (xvec,yvec) = angles2xy(azmat.flatten(),elmat.flatten())
    ksysout = griddata(points, ksys, (xvec, yvec), method='nearest')
    # size ksys to deal with the different wavelength and antenna gain from Millstone
    # makes the returns go up by 2.5 dB
    ksysout = ksysout*np.power(10.,2.5/10.)
    beamcodemapP = np.column_stack((np.arange(azmat.size),azmat.flatten(),elmat.flatten(),ksysout))


    #%% Write Millstone parameters
    h5fileout = tables.open_file('Millstone_PARAMS.h5',mode='w',title='Millstone Parameters')
    fgroup = h5fileout.create_group('/','Params','Parameters')

    h5fileout.create_array(fgroup, 'Frequency', txfreqP)
    h5fileout.create_array(fgroup, 'Power', txpowP)
    h5fileout.create_array(fgroup, 'Kmat', beamcodemapP)
    h5fileout.create_array(fgroup,'Systemp',caltempP)
    h5fileout.create_array(fgroup,'Bandwidth',bandwidthP)
    h5fileout.create_array(fgroup,'Sampletime',sampletimeP)
    h5fileout.create_array(fgroup,'Angleoffset',angoffP)
    h5fileout.close()

