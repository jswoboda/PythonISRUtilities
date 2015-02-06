#!/usr/bin/env python
"""
Created on Sat Mar 15 10:21:20 2014

@author: John Swoboda

This script is here to create the h5 files for the ISR systems
"""

import os
import numpy as np
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
    h5file = tables.open_file(curfile)
    # get the tx power for RISR
    txpowR = h5file.getNode('/Tx/Power').read()[0,0]
    # get freq 1 for RISR
    txfreqR = h5file.getNode('/Tx/Frequency').read()[0,0]
    # get cal temp
    caltempR = h5file.getNode('/Rx/CalTemp').read()
    # get bandwidth
    bandwidthR = h5file.getNode('/Rx/Bandwidth').read()
    # sample time
    sampletimeR = h5file.getNode('/Rx/SampleTime').read()
    # angle offset
    angoffR = np.array([26.0,45.0])

    beamcodemapR = h5file.getNode('/Setup/BeamcodeMap').read()
    beamcodemapR = np.array(beamcodemapR)

    h5file.close()
    #%% Write RISR parameters
    h5fileout = tables.open_file('RISR_PARAMS.h5',mode='w',title='RISR Parameters')
    fgroup = h5fileout.create_group('/','Params','Parameters')

    h5fileout.create_array(fgroup, 'Frequency', txfreqR)
    h5fileout.create_array(fgroup, 'Power', txpowR)
    h5fileout.create_array(fgroup, 'Kmat', beamcodemapR)
    h5fileout.create_array(fgroup,'Systemp',caltempR)
    h5fileout.create_array(fgroup,'Bandwidth',bandwidthR)
    h5fileout.create_array(fgroup,'Sampletime',sampletimeR)
    h5fileout.create_array(fgroup,'Angleoffset',angoffR)

    h5fileout.close()
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
    caltempP = h5file.getNode('/Rx/CalTemp').read()
    # get bandwidth
    bandwidthP = h5file.getNode('/Rx/Bandwidth').read()
    # sample time
    sampletimeP = h5file.getNode('/Rx/SampleTime').read()
    # angle offset
    angoffP = np.array([15.0,16.0])

    beamcodemapP = h5file.getNode('/Setup/BeamcodeMap').read()
    beamcodemapP = np.array(beamcodemapR)
    h5file.close()
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

