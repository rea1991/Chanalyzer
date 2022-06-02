#!/bin/bash


cd ./example/IonicChannels/mscl/5ligand/frame.001


rm -rf ./ED
mkdir ED ED/shape_descriptors ED/ED_CCs ED/ED_tet ED/mouths
mkdir  ED/mouths/0 ED/mouths/1 ED/mouths/2 ED/mouths/3 ED/mouths/4 ED/mouths/5 ED/mouths/6 ED/mouths/7 ED/mouths/8 ED/mouths/9


#IonicChannels (MSCL)
../../../../../build/chanalyzer -conf conf.prm -find_geo_pockets 1 -num_c2r 10
