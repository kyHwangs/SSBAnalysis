#!/bin/tcsh


mkdir -p output
./ssb_analysis sample/UL2016APV_NNLO_inc.list UL2016APV_NNLO_inc_sample.root 0 false UL2016APV UL2016APV_NNLO_inc false false false 10 
./ssb_analysis sample/UL2016_NNLO_inc.list UL2016_NNLO_inc.root 0 false UL2016 UL2016_NNLO_inc false false false 10 
./ssb_analysis sample/UL2017_NNLO_inc.list UL2017_NNLO_inc.root 0 false UL2017 UL2017_NNLO_inc false false false 10 
./ssb_analysis sample/UL2018_NNLO_inc.list UL2018_NNLO_inc.root 0 false UL2018 UL2018_NNLO_inc false false false 10 
