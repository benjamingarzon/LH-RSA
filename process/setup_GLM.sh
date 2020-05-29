#!/bin/sh
NEVS_ORIG=$1
NEVS_REAL=$2
echo $NTRIALS
for index in `seq $NEVS_ORIG`; do
echo -e "# EV $index title
set fmri(evtitle$index) \"EV $index\"

# Basic waveform shape (EV $index)
# 0 : Square
# 1 : Sinusoid
# 2 : Custom (1 entry per volume)
# 3 : Custom (3 column format)
# 4 : Interaction
# 10 : Empty (all zeros)
set fmri(shape$index) 3

# Convolution (EV $index)
# 0 : None
# 1 : Gaussian
# 2 : Gamma
# 3 : Double-Gamma HRF
# 4 : Gamma basis functions
# 5 : Sine basis functions
# 6 : FIR basis functions
set fmri(convolve$index) 3

# Convolve phase (EV $index)
set fmri(convolve_phase$index) 0

# Apply temporal filtering (EV $index)
set fmri(tempfilt_yn$index) 1

# Add temporal derivative (EV $index)
set fmri(deriv_yn$index) 1

# Custom EV file (EV $index)
set fmri(custom$index) \"@EV${index}\"
"
for index2 in 0 `seq $NEVS_ORIG`; do
echo -e "# Orthogonalise EV $index wrt EV $index2
set fmri(ortho$index.$index2) 0
"
done
done > fsfchunk

echo -e "# Contrast & F-tests mode
# real : control real EVs
# orig : control original EVs
set fmri(con_mode_old) orig
set fmri(con_mode) orig

# Display images for contrast_real 1
set fmri(conpic_real.1) 0

# Title for contrast_real 1
set fmri(conname_real.1) \"Mean>0\"
" >> fsfchunk

echo -e "# Real contrast_real vector 1 element 1
set fmri(con_real1.1) 1
" >> fsfchunk

for index in `seq 2 $NEVS_REAL`; do
echo -e "# Real contrast_real vector 1 element $index
set fmri(con_real1.$index) 0
"; done >> fsfchunk

echo -e "# Title for contrast_orig 1
set fmri(conname_orig.1) \"Mean>0\"
" >> fsfchunk

echo -e "# Real contrast_orig vector 1 element 1
set fmri(con_orig1.1) 1
" >> fsfchunk

for index in `seq 2 $NEVS_ORIG`; do
echo -e "# Real contrast_orig vector 1 element $index
set fmri(con_orig1.$index) 0
"; done >> fsfchunk



