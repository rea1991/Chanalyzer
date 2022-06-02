#!/bin/bash


cd ./example/IonicChannels/mscl/5ligand/frame.002


rm -rf ./ED
mkdir ED ED/ED_CCs ED/ED_tet


#IonicChannels (MSCL)
../../../../../build/chanalyzer -conf conf.prm -find_geo_pockets 1 -num_c2r 10
