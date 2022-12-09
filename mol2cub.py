"""
mol2cub
A tool to generate an ESP cube file from a mol2 file containing charged atoms.

Usage:
1) 
$ python mol2cub.py input.mol2 output.cube

In this default mode, a box is generated that bounds the molecule plus a margin of 10 Bohr. 
The resolution of points is 0.5 Bohr along each axis.

2)
$ python mol2cub.py input.mol2 output.cube head

head - the head section of a cube file (at least the first 6 lines) specifying 
cooridnates of the box origin, number of points and spacing along each direction 


Last update: 2 December 2022
"""

import sys, re, numpy, scipy, scipy.spatial
from decimal import Decimal

bohr = 1.88972599
read_head = False # the default mode - without reading the head file

# If needed add missing elements and their atomic numbers to this dictionary
atomicnumbersdict = {
    "H": 1, "He": 2,
    "Li": 3, "Be":4, "B":5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10,
    "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18,
    "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27,
    "Ni": 28, "Cu": 29, "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36,
    "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46,
    "Ag": 47, "Cd": 48, "In": 49, "Sn": 50, "Sb": 51, "Te": 52, "I": 53, "Xe": 54,
    "Cs": 55, "Ba": 56, "Hf": 72, "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, 
    "Au": 79, "Hg": 80, "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84, "At": 85, "Rn": 86,
    "Fr": 87, "Ra": 88, "E": 0
}

# Some function definitions
def frange(x, y, jump):
  while x < y:
    yield x
    x += jump

def frange2(x, N, jump):
    n = 0
    while n < N:
        yield x
        x += jump
        n += 1

def inv(x):
    if x != 0.000:
        return (1.0 / x)
    else:
        return float("NaN")

inv_v = numpy.vectorize(inv)


# File names read from command line
mol2_f_name = sys.argv[1]
output_f_name = sys.argv[2]
if len(sys.argv)>3:
    head_inp_file = sys.argv[3]
    read_head = True

if read_head:
    head_file = []
    box_orig = []
    with open (head_inp_file,'r') as file:
        for line in file:
            head_file.append(line)
    Natoms_head = eval(head_file[2].split()[0])
    for i in range(3):
        box_orig.append( eval(head_file[2].split()[i+1]) )
    Nx = eval(head_file[3].split()[0])
    Ny = eval(head_file[4].split()[0])
    Nz = eval(head_file[5].split()[0])
    #
    xres = eval(head_file[3].split()[1])
    yres = eval(head_file[4].split()[2])
    zres = eval(head_file[5].split()[3])
    #
    xbox = list(frange2(box_orig[0], Nx, xres))
    ybox = list(frange2(box_orig[1], Ny, yres))
    zbox = list(frange2(box_orig[2], Nz, zres))
else:
# Default Resolution of the resulting cubefile (in bohr)
    res = 0.5
    xres = res
    yres = res
    zres = res

mol2file = []
atoms = []

# Open the mol2 file as a list of strings corresponding to each line
with open (mol2_f_name,'r') as file:
    for line in file:
        mol2file.append(line)

# Find the beginning of the molecule specification in the mol2 file
molspecline = [i for i, item in enumerate(mol2file) if re.search('@<TRIPOS>MOLECULE', item)][0]

# Find the beginning of the atoms specification in the mol2 file
atomspecline = [i for i, item in enumerate(mol2file) if re.search('@<TRIPOS>ATOM', item)][0]

# The number of atoms should be the first integer 2 lines after the MOLECULE spec line
natoms = int(mol2file[molspecline + 2 ].split()[0])

atomnames = []
atomicnumbers = []
charges = []
xcoords = []
ycoords = []
zcoords = []

# Turn the atoms into a list of lists
# Also extracts the names, coordinates and charges
for i in range(0,natoms):
    atoms.append(mol2file[atomspecline + i + 1].split())
    atomnames.append(atoms[i][1])
    xcoords.append(float(atoms[i][2])*bohr)
    ycoords.append(float(atoms[i][3])*bohr)
    zcoords.append(float(atoms[i][4])*bohr)
    charges.append(float(atoms[i][8]))

# Find the atomic numbers for the atoms. New atom types can be added to the atomicnumbersdict dictionary
for name in atomnames:
    symbol = ''.join([i for i in name if not i.isdigit()])
    atomicnumbers.append(atomicnumbersdict[symbol])

# Write a little header for the cubefile output
with open (output_f_name,'w') as file:
    file.write("Potential from mol2\nMade using python\n")

# Determine the extents of the box
if not read_head:
    xbox = list(frange(min(xcoords)-10, max(xcoords)+10, xres))
    ybox = list(frange(min(ycoords)-10, max(ycoords)+10, yres))
    zbox = list(frange(min(zcoords)-10, max(zcoords)+10, zres))

# Write the extent of the box to the cubefile
with open (output_f_name,'a') as file:
    # Write the origin
    file.write("  " + str(natoms) 
    + "  " + str(xbox[0]) 
    + "  " + str(ybox[0])
    + "  " + str(zbox[0]) + "\n")    
    # Write the resolution and extents
    file.write("  " + str(len(xbox)) + " " + str(xres) + " 0.0 0.0 \n")
    file.write("  " + str(len(ybox)) + " 0.0 " + str(yres) + " 0.0 \n")
    file.write("  " + str(len(zbox)) + " 0.0 0.0 " + str(zres) + " \n")

# Write the atomic coordinates
with open (output_f_name,'a') as file:
    for i in range(0,natoms):
        file.write("  " + str(atomicnumbers[i]) + "  0.0  " + str(round(xcoords[i],6)) + \
                   "  " + str(round(ycoords[i],6)) + "  " + str(round(zcoords[i],6)) + " \n")
  

# generate an array of vectors - points of the box
temp = []
for x in xbox:
    for y in ybox:
        for z in zbox:
            temp.append([x, y, z])
points_array = numpy.array(temp, dtype=float)


# calculate esp values at the points of the box
esp_array = numpy.zeros(len(points_array), dtype=float)

for i in range(0,natoms):
    at_coords = numpy.array([xcoords[i], ycoords[i], zcoords[i]], dtype=float)
    at_coords_ex = numpy.expand_dims(at_coords, 0)
    dist = scipy.spatial.distance.cdist(at_coords_ex, points_array)
    dist_1 = inv_v(dist[0])
    esp_array +=  charges[i] * dist_1


# Write the calculated ESP values to the cubefile
with open (output_f_name,'a') as file:
     n=0
     c=0
     for esp in esp_array:
         if n > 5:
             file.write('\n')
             n=0
         file.write(' ' + str('%.5E' % Decimal(esp)) + ' ')
         n=n+1
         c=c+1
     file.write('\n\n') 

