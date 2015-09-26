# Molden2AIM
Molden2AIM is a utility program which can be used to creat AIM-WFN and NBO-47 files from a Molden file.

## Latest Version
Version 3.1.0 (02/25/2015). Significant updates since version 3.0.0.

1. Supports the MOLDEN files generated by NWChem, BDF, PSI4 (spherical functions only), CADPAC, and MRCC.
2. Bug fix for NBO6.
3. ReOrdAtm.f90 has been updated for CFour.
4. Check the normalization of NBO's *.47 file.

## Features

1. It converts the data format from Molden to AIM's WFN. The latter format can be read by AIMPAC, AIMPAC2, AIM2000, AIMALL, AIM-UC, DGrid, MORPHY98, Multiwfn, PAMoC, ProMolden, TopChem, TopMoD, Xaim, and so on.
2. It saves NBO's *.47 data file. One can do NBO analysis using the stand-alone GENNBO program.

The WFX format will be supported in future.

## About the Molden file

See readme.html for details.
