#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

def CdSe():
    """Call the CdSe Physical properties (lattice parameter, C11, volumic density
    , deduced longitudinal sound velocity)"""
    lattice_para = 6.08    #Lattice parameter Angstrom
    elastic_const = 88.    #Elastic constant C11 (GPa)
    density = 5655.    #Volumic density kg.m-3
    L_Sound_Velocity = np.sqrt((elastic_const*1E9)/(density))
    mat_datas = (lattice_para,elastic_const, density, L_Sound_Velocity)
    return mat_datas

def CdS():
    """Call the CdS Physical properties (lattice parameter, C11, volumic density
    , deduced longitudinal sound velocity)"""
    lattice_para = 5.82    #Lattice parameter Angstrom
    elastic_const = 98.    #Elastic constant C11 (GPa)
    density = 4870.    #Volumic density kg.m-3
    L_Sound_Velocity = np.sqrt((elastic_const*1E9)/(density))
    mat_datas = (lattice_para,elastic_const, density, L_Sound_Velocity)
    return mat_datas

def h_thick(layer_number,lattice_para):
    """Deduce the thickness in meter from the number of layers"""
    #ajouter une ligne pour demander le nombre de layer
    return layer_number*((lattice_para*1e-10)/2)

def naked_vib(h,V_L):
    """Calculate the vibration of a 'free' nanoplatelets of thickness h in wavenumber(cm-1)"""
    nu_naked = ((1/(2*h))*V_L)/(30*1e9)
    return nu_naked

#def mass_cos_func(omega,h,)
