#!/usr/bin/sh

PHASE=8
#./pipeline.sh 104 1 trialsfile-lup2s008_fmri.csv 28 $PHASE "1 2 3" organize_responses_v0.py "2.8 0.5 0 25 6" "28,29,30,31,32,33"&

PHASE=3
#./pipeline.sh 105 3 trialsfile-lup2s007_fmri.csv 32 $PHASE "1 2 3 4 5" organize_responses.py "2.1 0.5 0 6.0" "26,27,28,29,30,31"&

PHASE=6
#./pipeline.sh 105 3 trialsfile-lup2s007_fmri.csv 32 $PHASE "1 2 3 4 5" organize_responses.py "2.1 0.5 0 6.0" "26,27,28,29,30,31"&

#./pipeline.sh 105 2 trialsfile-lup2s007_fmri.csv 28 $PHASE "1 2 3 4" organize_responses_v0.py "2.8 0.5 0 25 6" "28,29,30,31,32,33"&
#./pipeline.sh 105 1 trialsfile-lup2s007_fmri.csv 28 $PHASE "1 2 3 4" organize_responses_v0.py "2.8 0.5 0 25 6" "28,29,30,31,32,33"&
#./pipeline.sh 106 1 trialsfile-lup2s009_fmri.csv 32 $PHASE "1 2 3 4 5" organize_responses.py "2.1 0.5 0 6.0" "28,29,30,31,32,33"&

#./pipeline.sh 103 1 trialsfile-lup2s002_fmri.csv 28 $PHASE "1 2 3" organize_responses_v0.py "2.8 0.5 0 25 6" "28,29,30,31,32,33"&
PHASE=7
#PHASE=6
#./pipeline.sh 106 1 trialsfile-lup2s009_fmri.csv 32 $PHASE "2 3 4 5" organize_responses.py "2.1 0.5 0 6.0" "26,27,28,29,30,31"&

PHASE=5
#./pipeline.sh 107 1 trialsfile-lup2s001_fmri.csv 32 $PHASE "1 2 3" organize_responses.py "2.1 0.5 0 6.0" "26,27,28,29,30,31"&
#./pipeline.sh 106 1 trialsfile-lup2s009_fmri.csv 32 $PHASE "2 3 4 5" organize_responses.py "2.1 0.5 0 6.0" "26,27,28,29,30,31"&

PHASE=6
#./pipeline.sh 106 1 trialsfile-lup2s009_fmri.csv 32 $PHASE "2 3 4 5" organize_responses.py "2.1 0.5 0 6.0" "26,27,28,29,30,31"&

#./pipeline.sh 103 1 trialsfile-lup2s002_fmri.csv 28 $PHASE "1 2 3" organize_responses_v0.py "2.8 0.5 0 25 6" "28,29,30,31,32,33"
PHASE=8
#./pipeline.sh 103 3 trialsfile-lup2s002_fmri.csv 32 $PHASE "1 2 3 4 5" organize_responses.py "2.1 0.5 0 6.0" "26,27,28,29,30,31"&
#./pipeline.sh 106 1 trialsfile-lup2s009_fmri.csv 32 $PHASE "2 3 4 5" organize_responses.py "2.1 0.5 0 6.0" "26,27,28,29,30,31"
#./pipeline.sh 103 1 trialsfile-lup2s002_fmri.csv 28 $PHASE "1 2 3" organize_responses_v0.py "2.8 0.5 0 25 6" "28,29,30,31,32,33"&
./pipeline.sh 107 1 trialsfile-lup2s001_fmri.csv 32 $PHASE "1 2 3 4 5" organize_responses.py "2.1 0.5 0 6.0" "26,27,28,29,30,31"
./pipeline.sh 107 2 trialsfile-lup2s001_fmri.csv 32 $PHASE "1 2 3 4" organize_responses.py "2.1 0.5 0 6.0" "26,27,28,29,30,31"&
#./pipeline.sh 105 3 trialsfile-lup2s007_fmri.csv 32 $PHASE "1 2 3 4 5" organize_responses.py "2.1 0.5 0 6.0" "26,27,28,29,30,31"&

