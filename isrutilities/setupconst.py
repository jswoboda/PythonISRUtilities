#!/usr/bin/env python
"""
Created on Sat Mar 15 10:21:20 2014

@author: John Swoboda

This script is here to create the h5 files for the ISR systems
"""

import os
import numpy as np
import yaml
from isrutilities import Path
from mathutils import angles2xy
import physConstants as Const
from scipy.interpolate import griddata


import tables
import glob


if __name__ == "__main__":
    dirname = Path(__file__).expanduser().parent.parent
    # Get only .dt0 files
    flist1 = glob.glob(os.path.expanduser('~')+'/DATA/20140224.002/'+'*.dt0.h5')
    flist2 = glob.glob(os.path.expanduser('~')+'/DATA/20111215.022/'+'*.dt0.h5')
    #%% Read RISR Params
    # open file
    caltempR = 120.
    angoffR = np.array([26.0,45.0])
    try:

        curfile = flist1[0]
        filepath, fext = os.path.splitext(curfile)
        with tables.open_file(curfile) as f:
            # get the tx power for RISR
            txpowR = float(f.get_node('/Tx/Power').read()[0,0])
            # get freq 1 for RISR
            txfreqR = float(f.get_node('/Tx/Frequency').read()[0,0])
            # get cal temp
        #    caltempR = f.getNode('/Rx/CalTemp').read()
            # get bandwidth
            bandwidthR = float(f.get_node('/Rx/Bandwidth').read())
            # sample time
            sampletimeR = float(f.get_node('/Rx/SampleTime').read())
            # angle offset

            beamcodemapR = f.get_node('/Setup/BeamcodeMap').read()
    except EnvironmentError:

        newfile = dirname.joinpath('AMISR_Params.h5')
        with tables.open_file(str(newfile)) as f:

            txpowR = 1.72e6
            # get freq 1 for RISR
            txfreqR = 440e6

        #    caltempR = f.getNode('/Rx/CalTemp').read()
            # get bandwidth
            bandwidthR = 22970.26336932706
            # sample time
            sampletimeR = 2e-5
            # angle offset

            beamcodemapR = f.get_node('/RISRNKsys').read()

    beamcodemapR = np.array(beamcodemapR)

    posel = np.where(beamcodemapR[:,2]>0)
    beamcodemapR = beamcodemapR[posel]
    #%% Write RISR parameters
    RISR_file = str(dirname.joinpath('RISR_PARAMS.h5'))
    with tables.open_file(RISR_file, mode='w', title='RISR Parameters') as f:
        fgroup = f.create_group('/', 'Params', 'Parameters')

        f.create_array(fgroup, 'Frequency', txfreqR)
        f.create_array(fgroup, 'Power', txpowR)
        f.create_array(fgroup, 'Kmat', beamcodemapR)
        f.create_array(fgroup, 'Systemp', caltempR)
        f.create_array(fgroup, 'Bandwidth', bandwidthR)
        f.create_array(fgroup, 'Sampletime', sampletimeR)
        f.create_array(fgroup, 'Angleoffset', angoffR)


    #%% Set up PFISR parameters
    #hard code in because PFISR example is 3x to high
    caltempP = 120
    angoffP = np.array([15.0,16.0])

    try:
        curfile = flist2[0]
        filepath, fext = os.path.splitext(curfile)
        # open file
        with tables.open_file(curfile) as f:
        # get the tx power for RISR
            txpowP = f.get_node('/Tx/Power').read()[0,0]
            # get freq 1 for RISR
            txfreqP = f.get_node('/Tx/Frequency').read()[0,0]
            # get cal temp

        #    caltempP = h5file.getNode('/Rx/CalTemp').read()
            # get bandwidth
            bandwidthP = f.get_node('/Rx/Bandwidth').read()
            # sample time
            sampletimeP = f.get_node('/Rx/SampleTime').read()
            # angle offset

            beamcodemapP = f.get_node('/Setup/BeamcodeMap').read()


    except EnvironmentError:
        dirname = Path(__file__).expanduser().parent.parent
        newfile = dirname.joinpath('AMISR_Params.h5')
        with tables.open_file(str(newfile)) as f:

            txpowP = 1.72e6
            # get freq 1 for RISR
            txfreqP = 440e6
            # get bandwidth
            bandwidthP = 22970.26336932706
            # sample time
            sampletimeP = 2e-5
            # angle offset
            beamcodemapP= f.get_node('/PFISRKsys').read()
    beamcodemapP = np.array(beamcodemapP)
    posel = np.where(beamcodemapP[:,2]>0)
    beamcodemapP = beamcodemapP[posel]
    #%% Write PFISR parameters
    PFISR_file = str(dirname.joinpath('PFISR_PARAMS.h5'))
    with tables.open_file(PFISR_file, mode='w', title='PFISR Parameters') as h5fileout:
        fgroup = h5fileout.create_group('/','Params','Parameters')

        h5fileout.create_array(fgroup, 'Frequency', txfreqP)
        h5fileout.create_array(fgroup, 'Power', txpowP)
        h5fileout.create_array(fgroup, 'Kmat', beamcodemapP)
        h5fileout.create_array(fgroup, 'Systemp', caltempP)
        h5fileout.create_array(fgroup, 'Bandwidth', bandwidthP)
        h5fileout.create_array(fgroup, 'Sampletime', sampletimeP)
        h5fileout.create_array(fgroup, 'Angleoffset', angoffP)


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

    # Gain is 49 dB according to data sheet
    G = np.power(10,4.9)
    lam = .2323
    ksysc = G*lam**2*Const.v_C_0*Const.v_electronradius**2/2./4**2/np.pi**2
    # angle offset
    angoffP = np.array([0.0,0.0])

    azvec = np.linspace(0.,360.,361)
    elvec = np.linspace(25.,89,65)
    (azmat,elmat) = np.meshgrid(azvec,elvec)
    azflat = np.append(azmat.flatten(),0.)
    elflat = np.append(elmat.flatten(),90.)
    ksysout = ksysc*np.ones_like(azflat)

    beamcodemapP = np.column_stack((np.arange(azflat.size),azflat,elflat,ksysout))


    #%% Write Sondrestrom parameters
    sondfile = str(dirname.joinpath('Sondrestrom_PARAMS.h5'))
    with tables.open_file(sondfile, mode='w', title='Sondrestrom Parameters') as h5fileout:
        fgroup = h5fileout.create_group('/', 'Params', 'Parameters')

        h5fileout.create_array(fgroup, 'Frequency', txfreqP)
        h5fileout.create_array(fgroup, 'Power', txpowP)
        h5fileout.create_array(fgroup, 'Kmat', beamcodemapP)
        h5fileout.create_array(fgroup, 'Systemp', caltempP)
        h5fileout.create_array(fgroup, 'Bandwidth', bandwidthP)
        h5fileout.create_array(fgroup, 'Sampletime', sampletimeP)
        h5fileout.create_array(fgroup, 'Angleoffset', angoffP)

    #%% Set up Millstone parameters

    # set the tx power for Millstone
    txpowP = 2.5e6
    # set freq 1 for Millstone
    txfreqP = 4.4e8
    # set cal temp
    caltempP = 120.
    # get bandwidth
    bandwidthP = 50e3
    # sample time
    sampletimeP = 1./50e3
    # angle offset

    #MISA gain is 42.5 dB and Zenith is 45.5
    G = np.power(10,4.25)
    lam =  0.6813464954545455
    ksysc = G*lam**2*Const.v_C_0*Const.v_electronradius**2/2./4**2/np.pi**2
    # angle offset
    angoffP = np.array([0.0,0.0])

    azvec = np.linspace(0.,360.,361)
    elvec = np.linspace(25.,89,65)
    (azmat,elmat) = np.meshgrid(azvec,elvec)
    azflat = np.append(azmat.flatten(),0.)
    elflat = np.append(elmat.flatten(),90.)
    ksysout = ksysc*np.ones_like(azflat)

    beamcodemapP = np.column_stack((np.arange(azflat.size),azflat,elflat,ksysout))


    #%% Write Millstone parameters
    millfile = str(dirname.joinpath('Millstone_PARAMS.h5'))
    with tables.open_file(millfile, mode='w', title='Millstone Parameters') as h5fileout:
        fgroup = h5fileout.create_group('/','Params','Parameters')

        h5fileout.create_array(fgroup, 'Frequency', txfreqP)
        h5fileout.create_array(fgroup, 'Power', txpowP)
        h5fileout.create_array(fgroup, 'Kmat', beamcodemapP)
        h5fileout.create_array(fgroup,'Systemp',caltempP)
        h5fileout.create_array(fgroup,'Bandwidth',bandwidthP)
        h5fileout.create_array(fgroup,'Sampletime',sampletimeP)
        h5fileout.create_array(fgroup,'Angleoffset',angoffP)

    #%% Set up Millstone parameters

    # set the tx power for Millstone
    txpowP = 2.5e6
    # set freq 1 for Millstone
    txfreqP = 4.4e8
    # set cal temp
    caltempP = 120.
    # get bandwidth
    bandwidthP = 50e3
    # sample time
    sampletimeP = 1./50e3
    # angle offset

    #MISA gain is 42.5 dB and Zenith is 45.5
    G = np.power(10,4.55)
    lam =  0.6813464954545455
    ksysc = G*lam**2*Const.v_C_0*Const.v_electronradius**2/2./4**2/np.pi**2
    # angle offset
    angoffP = np.array([0.0,0.0])

    azvec = np.linspace(0.,360.,361)
    elvec = np.linspace(25.,89,65)
    (azmat,elmat) = np.meshgrid(azvec,elvec)
    azflat = np.append(azmat.flatten(),0.)
    elflat = np.append(elmat.flatten(),90.)
    ksysout = ksysc*np.ones_like(azflat)

    beamcodemapP = np.column_stack((np.arange(azflat.size),azflat,elflat,ksysout))


    #%% Write Millstone parameters
    millfile = str(dirname.joinpath('Millstonez_PARAMS.h5'))
    with tables.open_file(millfile, mode='w', title='Millstone Parameters') as h5fileout:
        fgroup = h5fileout.create_group('/','Params','Parameters')

        h5fileout.create_array(fgroup, 'Frequency', txfreqP)
        h5fileout.create_array(fgroup, 'Power', txpowP)
        h5fileout.create_array(fgroup, 'Kmat', beamcodemapP)
        h5fileout.create_array(fgroup,'Systemp',caltempP)
        h5fileout.create_array(fgroup,'Bandwidth',bandwidthP)
        h5fileout.create_array(fgroup,'Sampletime',sampletimeP)
        h5fileout.create_array(fgroup,'Angleoffset',angoffP)
