import sys,re,numpy,math
from decimal import Decimal

# Resolution of the resulting cubefile (in bohr)
res = 0.5

xres = res
yres = res
zres = res

bohr = 1.8897161646320724

# Add other elements and their atomic numbers to this list
atomicnumbersdict = {
    "H": 1,
    "C": 6,
    "N": 7,
    "O": 8,
    "Se": 34,
    "E": 0
}

# Some function definitions
def frange(x, y, jump):
  while x < y:
    yield x
    x += jump

def distance(x0,y0,z0,x1,y1,z1):
    return math.sqrt((x0-x1)**2 + (y0-y1)**2 + (z0-z1)**2)

mol2file = []
atoms = []

# Open the mol2 file as a list of strings corresponding to each line
with open (sys.argv[1],'r') as file:
    for line in file:
        mol2file.append(line)

# Write a little header for the cubefile output
with open (sys.argv[2],'w') as file:
    file.write("Potential from mol2\nMade using python\n")
    
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

# Determine the extents of the box
boxpoints = []
xbox = list(frange(min(xcoords)-10, max(xcoords)+10, xres))
ybox = list(frange(min(ycoords)-10, max(ycoords)+10, yres))
zbox = list(frange(min(zcoords)-10, max(zcoords)+10, zres))

# Write the extent of the box to the cubefile
with open (sys.argv[2],'a') as file:
    # Write the origin
    file.write("  " + str(natoms) 
    + "  " + str(min(xcoords)-10) 
    + "  " + str(min(ycoords)-10)
    + "  " + str(min(zcoords)-10) + "\n")
    # Write the resolution and extents
    file.write("  " + str(len(xbox)) + " " + str(xres) + " 0.0 0.0 \n")
    file.write("  " + str(len(ybox)) + " 0.0 " + str(yres) + " 0.0 \n")
    file.write("  " + str(len(zbox)) + " 0.0 0.0 " + str(zres) + " \n")

# Write the atomic coordinates
with open (sys.argv[2],'a') as file:
    for i in range(0,natoms):
        file.write("  " + str(atomicnumbers[i]) + "  0.0  " + str(round(xcoords[i],6)) + "  " + str(round(ycoords[i],6)) + "  " + str(round(zcoords[i],6)) + " \n")
  

for x in xbox:
    for y in ybox:
        for z in zbox:
            #esp calculating code goes here
            esp = 0
            for i in range(0,natoms):
                esp = esp + charges[i]/distance(x,y,z,xcoords[i],ycoords[i],zcoords[i])
            boxpoints.append([x,y,z,esp])

# Write the calculated ESP values to the cubefile
with open (sys.argv[2],'a') as file:
    n=0
    c=0
    for i in boxpoints:
        if n > 5:
            file.write('\n')
            n=0
        file.write(' ' + str('%.5E' % Decimal(i[3])) + ' ')
        n=n+1
        c=c+1
    file.write('\n\n')

