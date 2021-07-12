#!/bin/tcsh -f
setenv size 0
foreach k (`seq 1 25 33125`)


@ z = $k / 25 + 1
setenv size 0


echo "PROJECT: e1e" > jsub_new
echo "TRACK: simulation" >>jsub_new
echo "JOBNAME: two_pi_evt_${k}" >>jsub_new
echo "MAIL: skorodum@jlab.org" >>jsub_new
echo "SINGLE_JOB: TRUE" >>jsub_new
echo "MEMORY: 9500 MB" >>jsub_new
echo "DISK_SPACE: 30 GB" >>jsub_new
@ k1 = $k + 1
@ k2 = $k + 2
@ k3 = $k + 3
@ k4 = $k + 4
@ k5 = $k + 5
@ k6 = $k + 6
@ k7 = $k + 7
@ k8 = $k + 8
@ k9 = $k + 9
@ k10 = $k + 10
@ k11 = $k + 11
@ k12 = $k + 12
@ k13 = $k + 13
@ k14 = $k + 14
@ k15 = $k + 15
@ k16 = $k + 16
@ k17 = $k + 17
@ k18 = $k + 18
@ k19 = $k + 19
@ k20 = $k + 20
@ k21 = $k + 21
@ k22 = $k + 22
@ k23 = $k + 23
@ k24 = $k + 24
echo "INPUT_FILES: /mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/converted_July2018_cc_ok/out_newFrad_July2018_${k}.root  /mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/converted_July2018_cc_ok/out_newFrad_July2018_${k1}.root /mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/converted_July2018_cc_ok/out_newFrad_July2018_${k2}.root /mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/converted_July2018_cc_ok/out_newFrad_July2018_${k3}.root /mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/converted_July2018_cc_ok/out_newFrad_July2018_${k4}.root /mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/converted_July2018_cc_ok/out_newFrad_July2018_${k5}.root /mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/converted_July2018_cc_ok/out_newFrad_July2018_${k6}.root /mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/converted_July2018_cc_ok/out_newFrad_July2018_${k7}.root /mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/converted_July2018_cc_ok/out_newFrad_July2018_${k8}.root /mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/converted_July2018_cc_ok/out_newFrad_July2018_${k9}.root /mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/converted_July2018_cc_ok/out_newFrad_July2018_${k10}.root /mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/converted_July2018_cc_ok/out_newFrad_July2018_${k11}.root /mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/converted_July2018_cc_ok/out_newFrad_July2018_${k12}.root /mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/converted_July2018_cc_ok/out_newFrad_July2018_${k13}.root /mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/converted_July2018_cc_ok/out_newFrad_July2018_${k14}.root /mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/converted_July2018_cc_ok/out_newFrad_July2018_${k15}.root /mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/converted_July2018_cc_ok/out_newFrad_July2018_${k16}.root /mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/converted_July2018_cc_ok/out_newFrad_July2018_${k17}.root /mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/converted_July2018_cc_ok/out_newFrad_July2018_${k18}.root /mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/converted_July2018_cc_ok/out_newFrad_July2018_${k19}.root /mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/converted_July2018_cc_ok/out_newFrad_July2018_${k20}.root /mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/converted_July2018_cc_ok/out_newFrad_July2018_${k21}.root  /mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/converted_July2018_cc_ok/out_newFrad_July2018_${k22}.root  /mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/converted_July2018_cc_ok/out_newFrad_July2018_${k23}.root  /mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/converted_July2018_cc_ok/out_newFrad_July2018_${k24}.root /volatile/clas/clase1-6/skorodum/2pi_an_fin_check/two_pi_analysis_code/h10tot21_evt /volatile/clas/clase1-6/skorodum/2pi_an_fin_check/two_pi_analysis_code/norm_nphe_5july18_aft_cc_match_gt50.root /volatile/clas/clase1-6/skorodum/2pi_an_fin_check/two_pi_analysis_code/phel_integr_fract_12July18.txt" >>jsub_new
echo "COMMAND: /volatile/clas/clase1-6/skorodum/2pi_an_fin_check/two_pi_analysis_code/new_jsub/execut_high_evt.sh ${k}" >>jsub_new
echo "OUTPUT_DATA: out.root" >>jsub_new
echo "OUTPUT_TEMPLATE: /mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/aft_2pi_mymain_test9Jul2021_evt/fin_root9jul2021test_evt_${z}.root" >>jsub_new

#if ($size < 220000000) then
if (!(-e  /mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/aft_2pi_mymain_test9Jul2021_evt/fin_root9jul2021test_evt_${z}.root)) then

/site/bin/jsub jsub_new
echo ${z}

#echo $size"   "$z

#jremove /mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/converted_Aug2016/out_newFrad_Aug2016_${z}.root
#rm -f /cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/converted_Aug2016/out_newFrad_Aug2016_${z}.root
endif
#rm jsub_new

end
