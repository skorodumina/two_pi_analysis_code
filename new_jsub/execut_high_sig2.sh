#!/bin/tcsh -f

setenv dir `pwd`
#cd /apps/root/PRO/root/
cd /apps/root/5.34.36/root/
source bin/thisroot.csh
cd $dir

echo "1" > inp
echo "2.039" >> inp
echo "1" >> inp
echo "/cache/mss/home/skorodum/e1e/data_2pi_conv_2July2018/out_conv_2pi_2July2018_9.root" >> inp
echo "1" >> inp
echo "/cache/mss/home/skorodum/e1e/data_2pi_conv_2July2018/out_conv_2pi_2July2018_9.root" >> inp
echo "25" >> inp
@ z = $1 + 24
foreach k (`seq  $1 $z`)
echo "out_newFrad_July2018_${k}.root" >> inp
end
echo "out.root" >> inp

./h10tot21_sig2<inp
