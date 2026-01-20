

# HERO: A Memory-Efficient and Fast Structured Low-Rank Reconstruction Algorithm for Accelerated Hypercomplex NMR

**HERO** (a **H**ankel-based **E**fficient **R**econstructions with **O**ptimization) is a fast and memory-saving structured low-rank algorithm for reconstructing hypercomplex NMR spectra. 

**This code is attached for the review of our related manuscript. We kindly request that it not be distributed or shared publicly until the paper is officially published.** 

# Content

data (Five synthetic hypercomplex data of different sizes and the corresponding non-uniform sampling schedules)



function 

├─ AhHX_Generate_pro.m (Optimization script for Hankel matrix right multiplication operation) 

├─ HhAB_Generate.m (Optimization script for Hankel matrix adjoint operations)

├─ HXBh_Generate_pro.m (Optimization script for Hankel matrix left multiplication operation) 

├─ NUS2D_HERO.m (The running script of the HERO algorithm)

├─ NUS2D_MF.m (The running script of the MF algorithm)

└─ NUS2D_SVD.m (The running script of the SVD algorithm)



LowRank_Toolbox

├─ BHankel2Matrix.m (The running script of the adjoint operator of the block-Hankel transformation)

├─ Hankel2vec.m (The running script of the adjoint operator of the Hankel transformation)

├─ Matrix2BHankel.m (The running script of the block-Hankel transformation)

├─ SVT.m (The script of Singular Value Thresholding)

└─ vec2Hankel.m  (The running script of the Hankel transformation)



README.md



test_2DNUS_simulated_slice.m (The demo is for running the program that reconstructs synthetic data. It effectively demonstrates the performance capabilities of HERO.)



# Instructions

Run the script "test_2DNUS_simulated_slice.m" to reconstruct the synthetic data using different algorithms (SVD, MF, HERO).



# Experimental Data

The experimental data reconstruction demo is not uploaded here due to the redistribution issue of experimental data. However, we believe that the simulated data experiments clearly demonstrate and validate the computational efficiency of the proposed algorithm. Interested users can download the data from the original BMRB data library (see below) and use the simulated data reconstruction demo as a reference. 

1. CPS_2611: https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr15317/timedomain_data/nesgCsR4_bmrb15317/csr4_n15_noesy_2_24_07.fid/
2. acpE: https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr17031/timedomain_data/bmrb17031td/sgr209c.013_d2o_hcch_tocsy_3_19_10.fid/
3. KPC1/RNF123: https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr19329/timedomain_data/bmrb19329td/HR8700A.005_152_D2O_Cnoesy_9_11_12.850mhz/
4. HR969: https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr6052/timedomain_data/hr969_4d_ccnoesy_7_21_03.fid/
5. IR24: https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr5842/timedomain_data/1Q48_bmrb5842.fids/



# Contact


If you have any questions or concerns about the code, please contact the author via email enpinglin@xmu.edu.cn. 




@article{Fang2026HERO,
  title   = {HERO: A Memory-Efficient and Fast Structured Low-Rank Reconstruction Algorithm for Accelerated Hypercomplex NMR},
  author  = {Fang, Ze and Chen, Yida and Luo, Yao and Yang, Yu and Lin, Enping and Chen, Zhong},
  journal = {Analytical Chemistry},
  volume  = {98},
  number  = {2},
  pages   = {1224--1232},
  year    = {2026},
  doi     = {10.1021/acs.analchem.5c05054},
  publisher = {American Chemical Society}
}


