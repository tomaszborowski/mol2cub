# mol2cub
A tool to generate an ESP cube file from a mol2 file containing charged atoms.

Usage:
1) 
$ python mol2cub.py input.mol2 output.cube

In this default mode, a box is generated that bounds the molecule plus a margin of 10 Bohr. 
The resolution of points is 0.5 Bohr along each axis.

2)
$ python mol2cub.py input.mol2 output.cube head

head - the head section of a cube file (at least the first 6 lines) specifying 
cooridnates of the box origin, number of points and spacing along each direction. 

The resulting cube file can be visualised using standard software.
