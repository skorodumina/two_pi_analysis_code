#!/bin/tcsh -f

setenv k 0
#foreach i (`seq 1 1000`)

#if (!(-e /cache/mss/home/gleb/e1e/sim2014/recsis_sim${i}.bos)) then
#if (!(-e /cache/mss/home/gleb/e1e/sim2014/ceb/out_ceb${i}.hbook)) then

#echo $i
#sed -e "s/test1/sim${i}/g" jsub_test > jsub${i}

#jsub jsub${i}

#rm jsub${i}

#@ j++

#endif

#end

#echo $j

#setenv k 0

foreach i (`seq 1 33000`)

#foreach file ( /cache/mss/home/gleb/e1e/sim2014_w_1675_1825/goa/goa_out*.hbook )

setenv size 0

if ((-e /mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/converted_Aug2016/out_newFrad_Aug2016_${i}.root)) then
setenv size `cat  /mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/converted_Aug2016/out_newFrad_Aug2016_${i}.root  | grep size | sed -e 's/size=//g'`
#echo $size" "$i
endif

if  ($size < 217000000) then

echo $size" "$i"  "
#jremove /mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/converted_Aug2016/out_newFrad_Aug2016_${i}.root
#rm -f /cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/converted_Aug2016/out_newFrad_Aug2016_${i}.root


@ k++

endif


end
echo $k
