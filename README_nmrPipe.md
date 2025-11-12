# Confidentiality Notice

The software and source code contained in this repository are provided solely for the purpose of peer review and evaluation in connection with our submitted manuscript. 

Redistribution, modification, or public dissemination of the code, in whole or in part, is strictly prohibited until the corresponding manuscript has been formally published. Unauthorized use of the code outside of the review process is not permitted.

Upon publication of the manuscript, the authors may release the code under an open-source license. Until then, all rights are reserved by the authors.

© 2025 Ze Fang, Yida Chen, Yao Luo, Yu Yang, Enping Lin*, and Zhong Chen. All rights reserved.

# Tutorial

All NMR source files used by the following program need to be downloaded from the BMRB. For specific addresses, please refer to the "Experimental Data" section of the README.md file. To avoid errors during program execution due to missing corresponding files, each script should be run in its designated location.

## The following steps are performed in the Linux system:

1. Open the terminal in the folder containing the original data.

2. Enter "csh".

3. Enter the manufacturer name corresponding to the original data file. If the original data file is fid, enter "varian"; if it is ser, enter "bruker".

4. In the pop-up window, click "Read Parameters" and "Save Script", then click "Quit".

5. Enter "fid.com".

6. Enter "process_direct" — this script performs Fourier transformation, phase correction, and removes the imaginary part along the direct dimension.

7. Enter "process_indirect.com" to obtain the projection spectrum of the original data.

## The following steps are performed in the Windows system:

8. Run "nmrPipe_to_mat_FID.py" — this script converts the FID data file from nmrPipe format to MATLAB format. Then run "nmrPipe_to_mat.py" — this script converts the projection spectrum data file from nmrPipe format to MATLAB format.

9. Run the corresponding reconstruction script.
10. Run "FID_to_nmrpipe_HERO.py" — this script converts the reconstructed data from MATLAB format back to nmrPipe format (the reconstructed data obtained using SVD and MF can be processed using their respective programs, and the same applies below).

## The following steps are performed in the Linux system:

11. In the terminal, enter "process_indirect_HERO.com" to obtain the projection spectrum.

## The following steps are performed in the Windows system:

12. Run "nmrPipe_to_mat_HERO.py" to obtain the projection spectrum data in MATLAB format. Typically, a single set of data has multiple projection files from different orientations. Users can modify the script to read the specific projection file as needed.



# Contact

If you have any questions or concerns about the code, please contact the author via email enpinglin@xmu.edu.cn. 