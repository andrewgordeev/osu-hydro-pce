# OSU hydro, PCE version

This is a modified version of the OSU hydro code available at https://github.com/jbernhard/osu-hydro

Usage is mostly idential to that code, so in lieu of a comprehensive README I have included the original README as READMEold.md in this folder.

There are a few key differences in usage:

- Following the standard CMake installation sequences produces a compiled binary osu-hydro-pce in <prefix>/bin/, instead of osu-hydro. The configuration file is similarly renamed osu-hydro-pce.conf. The command to run this program is thus osu-hydro-pce, not osu-hydro.
- Two equations of state, eos.dat and eos-gluon.dat, are produced at build time. The code interpolates between these with a fugacity factor.
- osu-hydro-pce.conf contains 2 new options: 
	1. Teq sets the chemical equilibration time used in determining the fugacity.
	2. Tfinal sets the final time at which the program will terminate regardless of the energy density (for testing purposes). By default this is 40 fm, the same value hard coded into osu-hydro.
- In this version the code outputs data for _all_ volume elements, not just the freezeout surface. Note that this greatly increases memory usage. This can be toggled by uncommenting the "If (intersect) THEN" line in FreezeoutPro10 in src/osu_hydro.for.
- The code outputs a Nx18 array of doubles, where N is the number of volume elements. Columns 0-15 are as described in READMEold.md, but the program additionally outputs the temperature of the volume elements in column 16 and the energy density in column 17.
