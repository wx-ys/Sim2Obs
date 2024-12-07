import os

import numpy as np



def make_header(particle,columns):
    trans = {'x':'Column $: position x (kpc)\n',
             'y':'Column $: position y (kpc)\n',
             'z':'Column $: position z (kpc)\n',
             'vx':'Column $: velocity vx (km/s)\n',
             'vy':'Column $: velocity vy (km/s)\n',
             'vz':'Column $: velocity vz (km/s)\n',
             'h':'Column $: smoothing length (kpc)\n',
             'mass':'Column $: mass (Msun)\n',
             'metals':'Column $: metallicity (1)\n',
             'temp':'Column $: temperature (K)\n',
             'mbirth':'Column $: initial mass (Msun)\n',
             'age':'Column $: age (Gyr)\n',
             'sfr':'Column $: star formation rate (Msun/yr)\n',
             'P':'Column $: pressure (J/m3)\n',
             'logC':'Column $: compactness (1)\n',
             'Fc':'Column $: coveringFactor (1)\n',
             'ncloud':'Column $: cloud number density (1/cm3)\n',
             'sfe':'Column $: star formation efficiency (1)\n',
             'electrons':'Column $: nr of electrons (1)\n',}
    header_con= ''
    for i in range(len(columns)):
        if i ==(len(columns)-1):
            header_con = header_con+trans[columns[i]].replace('$',str(i+1))[:-1]
        else:
            header_con = header_con+trans[columns[i]].replace('$',str(i+1))
    
    header = particle+'.txt: import file for '+particle+'\n' +header_con
    return header

def make_fmt(particle,columns):
    trans = {'x':'%.6f',
             'y':'%.6f',
             'z':'%.6f',
             'vx':'%.6f',
             'vy':'%.6f',
             'vz':'%.6f',
             'h':'%.6f',
             'mass':'%.1f',
             'metals':'%.8e',
             'temp': '%.1f',
             'mbirth':'%.1f',
             'age':'%.6f',
             'sfr':'%.6f',
             'P':'%.14e',
             'logC':'%.6f',
             'Fc':'%.4f',
             'ncloud':'%.6f',
             'sfe':'%.6f',
             'electrons':'%.2f',}
    fmt = [trans[i] for i in columns]
    return fmt

def readableFileSize(size: float) -> str:
    unit = ['Byte', 'KB', 'MB', 'GB', 'TB']
    for i in range(len(unit)):
        if size < 1024:
            return f'{size:.2f} {unit[i]}'
        size = size/1024.0
    return f'{size*1024:.2f} {unit[-1]}'

