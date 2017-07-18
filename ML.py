#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import scipy.constants as constant
import scipy.optimize as optimize
import matplotlib.pyplot as pyplot

pyplot.close("all")

def CdSe():
    """Call the CdSe Physical properties (lattice parameter, C11, volumic density
    , deduced longitudinal sound velocity)"""
    lattice_para = 6.08    #Lattice parameter Angstrom
    elastic_const = 88.    #Elastic constant C11 (GPa)
    density = 5655.    #Volumic density kg.m-3
    L_Sound_Velocity = np.sqrt((elastic_const*1E9)/(density))
    mat_datas = {'a_lat' : lattice_para, 'C_11' : elastic_const, 'rho' : density, 'v_long':L_Sound_Velocity}
    return mat_datas

def CdS():
    """Call the CdS Physical properties (lattice parameter, C11, volumic density
    , deduced longitudinal sound velocity)"""
    lattice_para = 5.82    #Lattice parameter Angstrom
    elastic_const = 98.    #Elastic constant C11 (GPa)
    density = 4870.    #Volumic density kg.m-3
    L_Sound_Velocity = np.sqrt((elastic_const*1E9)/(density))
    mat_datas = {'a_lat' : lattice_para, 'C_11' : elastic_const, 'rho' : density, 'v_long':L_Sound_Velocity}
    return mat_datas

def h_thick(layer_number):
    """Deduce the thickness in meter from the number of layers"""
    #ajouter une ligne pour demander le nombre de layer
    return layer_number*((CdS().get('a_lat')*1e-10)/2)

def naked_vib(n_layer):
    """Calculate the vibration of a 'free' nanoplatelets of thickness h in wavenumber(cm-1)"""
    nu_naked = (1/(2*h_thick(n_layer)))*(CdS().get('v_long'))/(30*1e9)
    return nu_naked

n_layer=int(input('entrer le nombre de layer:'))


print("frequence de vibration d'une plaquette nue avec {} layers".format(n_layer), naked_vib(n_layer))

def sigma_molecule(mass_mol):
    """Calculate the number density of ligands as a function of the Molecular mass (g/mol) and the lattice parameter (angstrom). Result express in g.m^-2"""
    return (2/((CdS().get('a_lat')*1e-10)**2))*(mass_mol*1e-3/constant.Avogadro)

def cos_mload(freq):
    """Calculate the cosinus term"""
    return np.cos( (2 * constant.pi * freq * h_thick(n_layer)) / (CdS().get('v_long') * 2))

def sin_mload(freq,mass_mol):
    """Calculate the sinusoidal term"""
    return ((sigma_molecule(mass_mol)*2 * constant.pi * freq)/(CdS().get('rho')*CdS().get('v_long'))) * np.sin( (2 * constant.pi * freq * h_thick(n_layer)) / (CdS().get('v_long') * 2))

n=301
result=[0 for j in range(n)]
mass_mol = np.linspace(0,300,301) #masse molaire variant de 0 à 300g/mol par pas de 1
for i in range (0,301): #dans le même range que la masse molaire calcul de la frequence de vibration pour chaque masse
    def zero_function(freq):
        """Calculate the zero of the function and gives the frequency of vibration"""
        return cos_mload(freq*30e9)-sin_mload(freq*30e9,i)
    result[i] = optimize.brentq(zero_function,10,50)
    print("frequence de vibration d'une plaquette avec des ligands en surface avec {} layers".format(n_layer),result)

wave_number = np.linspace (10,60,200)


#Plotting
#fig = P.figure()
#for i in range(1, 5):
#    ax = fig.add_subplot(2, 2, i)
#    j=2*i-2
#    ax.plot(data[:,j],data[:,j+1], 'ro')
#    ax.plot(x,(calc(j)[5]*x)+calc(j)[6],ls=':')
#    ax.set_title("Set({})".format(i))

fig, ax = pyplot.subplots()  # Création d'une figure contenant un seul système d'axes
ax.plot(mass_mol, result, c='b', ls='-', label="Sinus")    # Courbe y = sin(x)

#ax.plot(wave_number, sin_mload(wave_number*30e9), c='b', ls='-', label="Sinus")    # Courbe y = sin(x)
#ax.plot(wave_number, cos_mload(wave_number*30e9), c='r', ls=':', label="Cosinus")
#ax.plot(wave_number, zero_function(wave_number*30e9), c='g', ls='-', label="Cosinus")  # Courbe y = cos(x)
ax.set_xlabel("molecular wheight (g.mol$^{-1}$)")          # Nom de l'axe des x
ax.set_ylabel("frequency shift (cm$^{-1}$)")                # Nom de l'axe des y
#ax.set_title("Sinus et Cosinus")  # Titre de la figure
#ax.legend()                       # Légende
#fig.savefig("simple.png")         # Sauvegarde en PNG
pyplot.show()

data= mass_mol, result

#Sauvegarde du fichier représentant la fréquence en cm^-1 en fonction de la masse molaire des ligands
np.savetxt('CdS{}.txt'.format(n_layer),data, delimiter="    ")
