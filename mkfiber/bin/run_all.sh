#!/bin/bash
: '
This script generates 8 fibres with different number of motors to start with (80-110), and 9 fibres with 110 motor of different percent rotated (5%-95%).
As the fibres are generated, they are compressed, duplicated, solvated, and equilibriated. This yields a fibre with 220 motors (160-220 with varying starting motors), and a radius 2-3 nm.
There are 2 limiting factors:
    Packmoling the initial structure takes a long time, and sometimes yields unstable systems.
    NPT equillibration is run for 25 ns (12 500 000 steps). One 16 cores, with one GPU, this takes 4-5 hours.
Approximately 20-30% of the fibres generated from packmol are unstable, 
    and another 20-30% either do not equilibriate, or equilibriate to a very small radius as motors leave the fibres, again due to instabilities.
    These will have to be manually re-equillibrated from .tpr files in mkfiber/system/or completely re-run using test_frac_rot.sh and test_num_motors.sh.
'
while [[ $# -gt 0 ]]; do
  case $1 in
    -s|--subtype)
      shift
      SUBTYPE=$1
      shift
      ;;
    -gpu_id)
        shift
        GPU_ID=$1
        shift
        ;;
    -*|--*)
      echo "Unknown option $1"
      return 1
      ;;
    *)
      echo "Unknown value $1"
      return 1
      ;;
  esac
done

for frac in 96 91 75 50 25 10; do
    bash run_frac_rot.sh -s $SUBTYPE -gpu_id $GPU_ID -frac $frac -deffnm frac_${frac} #Run varying percent rotated
done

for nscale in 9 12 15; do
    bash run_frac_rot.sh -s $SUBTYPE -gpu_id $GPU_ID -frac 0 -n_scales $nscale -deffnm nscale_${nscale} #Run varying compression scales with 0% rotated
done
