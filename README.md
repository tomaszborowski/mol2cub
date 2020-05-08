# mol2cub
A tool to generate an ESP cube file from a mol2 file containing charged atoms.

Usage: $ python mol2cub.py input.mol2 output.cube

The resulting cube file can be visualised using standard software.

By default, a box is generated that bounds the molecule plus a margin of 10 Bohr. The resolution of points is 0.5 Bohr along each axis.
