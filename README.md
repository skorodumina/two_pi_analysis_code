# two_pi_analysis_code

Here is the analysis code that is used to analyse two pi data. As an input, files after the converter are used. As an output, a root file with multi dim hists is produced.

17feb17 - new electron selection, new files norm_nphe_27jan17_aft_cc_match_gt50.root and phel_integr_fract_22jan17.txt. New proton energy loss is added for data and sim. Electron correction is added for sim.

2018 -- new nphe files are produced after the cc plane match: norm_nphe_5july18_aft_cc_match_gt50.root and phel_integr_fract_12July18.txt.

13apr17 - new hadron id + new fid cuts for all particles, including th_vs_p + z-vertex cuts

15mar19 - adjusted pim fid and exclusivity cuts for both topologies in order to avoid accidental fluctuations in the single-differential distributions.

2021 -- this version of the code calculates the final version of the two pion cross sections.

Simulation output files that correspond to this code version are stored at

/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/aft_2pi_mymain_test1Jul2021_sig2/fin_root1jul2021test_*
/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/aft_2pi_mymain_test1Jul2021_evt/fin_root1jul2021test_evt_*

(1325 files in each directory)
-------------------------------


NOTE: for the empty target runs, the quality check cut is NOT needed.

-----------------------------

new_jsub/ folder contains script for batch farm job submission of the two pion code.
