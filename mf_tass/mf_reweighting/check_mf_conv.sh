#!/bin/bash/
for i in `seq 500000 500000 3500000`; do

  echo "working for $i number of MD steps"
  echo '*************************************************'
  
  j=$(bc<<<"$i/100000")
  cp run_reweight_conv.sh run_reweight_modified.sh
  sed -i -e "s/XXXX/${i}/g" run_reweight_modified.sh
  sh run_reweight_modified.sh
  mv free_energy.dat ${j}ns_free_energy.dat
  mv free_energy_1.dat ${j}ns_free_energy_1.dat
done
