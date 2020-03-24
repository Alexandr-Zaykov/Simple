New version of Simput (courtesy of Eric Buchanan)
Patches (by Alexandr Zaykov):
- Added cutoff possibilities to get rid of negligible orbitals - px, py in planar parts of the structure
    = Decent speedup (see Non-Planar test)
- Added "fidelity" patch for Simple 2.0+ to better reflect the actual structure of input files
- Changed the Pople's 6-311+G basis set to 6-311G due to diffuse function problems (artifacts)
- Patched wrong method format in output file
