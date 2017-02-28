#!/usr/bin/env python

"""
:platform: Unix, Windows, Mac
:synopsis: A set of physical constants that can be used by a bunch of different programs.
    table:: Variable Descriptions

    =================  ==============
    Variable Name      Description
    =================  ==============
    v_C_0  	          Speed of light m/s.
    v_Boltz            Boltzmans Const. J/K.
    v_me 	          Mass of electron in kg
    v_amu              Mass of proton/neutron in kg.
    v_electronradius   Electron radius in m.
    v_electron_rcs     Cross section of electon m^2.
    v_epsilon0         Vacuum permittivity.
    v_mu0              Vacuum permeability.
    v_elemcharge       Elementary charge in mks units
    v_elemcharge_cgs   Elementary charge in cgs units.
    R_EARTH            Radius of earth in km.
    ==============  ==============
        
.. moduleauthor:: John Swoboda <swoboj@bu.edu>
"""

import math

pi=math.pi

# some natural constants
v_C_0=299792458# speed of light m/s
v_Boltz=1.380658e-23 # Boltzmans Const. J/K
v_me=9.1093897e-31# mass of electron in kg
v_amu=1.6605402e-27# mass of proton/neutron in kg
v_electronradius=2.81794092e-15# electron radius in m
v_electron_rcs = 4*pi*v_electronradius**2
v_epsilon0=8.854188E-12#vacuum permittivity
v_mu0 = 4*pi * 1e-7;# H/m, m kg s^-2 A^-2 vacuum permeability
v_elemcharge=1.602177E-19 # elementary charge in mks
v_elemcharge_cgs=4.8032068e-10 #elementary charge in cgs
R_EARTH=6371.2# radius of earth in km
