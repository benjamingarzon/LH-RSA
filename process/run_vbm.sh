#!/bin/sh

NPROCS=10

VBM_DIR="/home/benjamin.garzon/Data/LeftHand/Lund1/vbm"
STRUC_DIR="$VBM_DIR/struc"

mkdir $STRUC_DIR
scan_list_file="$VBM_DIR/scans.txt"
subject_list_file="$VBM_DIR/subjects.txt"

printf '%s\n' sub-*/ses* > $scan_list_file
printf '%s\n' sub-* > $subject_list_file

cd $VBM_DIR

myimage=T1wtotemplate

scanlist=""
scans=`cat $scan_list_file`
#create links
for scan in $scans; do
    cd $VBM_DIR/${scan}
    NAME=`echo $scan | sed 's@/ses-@.@g'`
    echo $NAME    
    if [ -e "${myimage}.nii.gz" ]; then
        imcp $VBM_DIR/${scan}/${myimage}.nii.gz $STRUC_DIR/${NAME}_T1w_struc.nii.gz  
        imcp $VBM_DIR/${scan}/${myimage}_brain.nii.gz $STRUC_DIR/${NAME}_T1w_struc_brain.nii.gz

        scanlist="$scanlist $STRUC_DIR/${NAME}_T1w_struc $STRUC_DIR/${NAME}_T1w_struc_brain"
     fi   
    
done

echo $scanlist
echo $scanlist | wc -w

rm -r slicesdir
slicesdir -o $scanlist &

# create template with average images
subjects=`cat $subject_list_file`

cd $STRUC_DIR
T=${FSLDIR}/data/standard/tissuepriors/avg152T1_gray

### segmentation
rm -f fslvbm2a
for g in `imglob *_struc_brain.*` ; do
    if [ ! -e "${g}_GM.nii.gz" ]; then
        echo ${g}
        echo "fslreorient2std ${g} ${g}; \
        $FSLDIR/bin/fast -R 0.3 -H 0.1 ${g}; \
        $FSLDIR/bin/immv ${g}_pve_1 ${g}_GM >> fslvbm2a
     fi
done
chmod a+x fslvbm2a
cat fslvbm2a | xargs -n 1 -P $NPROCS -i -t sh -c "{}"

rm *pve* *seg* *mixel*
mkdir $VBM_DIR/templates

for subject in $subjects; do
  fslmerge -t $VBM_DIR/templates/${subject}_GM_all $STRUC_DIR/${subject}*_GM.nii.gz 
  fslmaths $VBM_DIR/templates/${subject}_GM_all -Tmean $VBM_DIR/templates/${subject}_GM
done

rm $VBM_DIR/templates/*_GM_all.nii.gz

cd $VBM_DIR/templates

### Estimation of the registration parameters of GM to grey matter standard template
/bin/rm -f fslvbm2b
for g in $subjects ; do
  if [ ! -e "${g}_GM_to_T.nii.gz" ]; then
      echo ${g}
      echo "${FSLDIR}/bin/fsl_reg ${g}_GM $T ${g}_GM_to_T -a" >> fslvbm2b
  fi
done
chmod a+x fslvbm2b
cat fslvbm2b | xargs -n 1 -P $NPROCS -i -t sh -c "{}"

echo Running initial registration 


### Creation of the GM template by averaging all (or following the template_list for) 
cat <<stage_tpl3 > fslvbm2c
#!/bin/sh
mergelist="*_GM_to_T.nii.gz"
\$FSLDIR/bin/fslmerge -t template_4D_GM \$mergelist
\$FSLDIR/bin/fslmaths template_4D_GM -Tmean template_GM
\$FSLDIR/bin/fslswapdim template_GM -x y z template_GM_flipped
\$FSLDIR/bin/fslmaths template_GM -add template_GM_flipped -div 2 template_GM_init
stage_tpl3
chmod +x fslvbm2c
echo Creating first-pass template 
./fslvbm2c

### Estimation of the registration parameters of GM to grey matter standard template
/bin/rm -f fslvbm2d
T=template_GM_init
for g in $subjects; do
  if [ ! -e "${g}_GM_to_T_init.nii.gz" ]; then
  echo ${g}
  echo "${FSLDIR}/bin/fsl_reg ${g}_GM $T ${g}_GM_to_T_init $REG -fnirt \"--config=GM_2_MNI152GM_2mm.cnf\"" >> fslvbm2d
  fi
done
chmod a+x fslvbm2d

echo Running registration to first-pass template
cat fslvbm2d | xargs -n 1 -P $NPROCS -i -t sh -c "{}"

cat <<stage_tpl4 > fslvbm2e
#!/bin/sh
mergelist="*_struc_GM_to_T"
\$FSLDIR/bin/fslmerge -t template_4D_GM \$mergelist
\$FSLDIR/bin/fslmaths template_4D_GM -Tmean template_GM
\$FSLDIR/bin/fslswapdim template_GM -x y z template_GM_flipped
\$FSLDIR/bin/fslmaths template_GM -add template_GM_flipped -div 2 template_GM
stage_tpl4
chmod +x fslvbm2e
echo Creating second-pass template
./fslvbm2e

echo "Study-specific template will be created"
echo "fsleyes " ${FSLDIR}/data/standard/tissuepriors/avg152T1_gray " template_GM"

ln -s $VBM_DIR/templates/template_GM.nii.gz $STRUC_DIR/template_GM.nii.gz

cd $STRUC_DIR

echo "Now running the preprocessing steps and the pre-analyses"
rm -f fslvbm3a
for scan in $scans; do
    g=`echo $scan | sed 's@/ses-@.@g'`
    echo $g    
  if [ ! -e "${g}_GM_to_template_GM_mod.nii.gz" ]; then
  echo $g
  echo "${FSLDIR}/bin/fsl_reg ${g}_T1w_struc_brain_GM template_GM ${g}_GM_to_template_GM -fnirt \"--config=GM_2_MNI152GM_2mm.cnf --jout=${g}_JAC_nl\"; \
        $FSLDIR/bin/fslmaths ${g}_GM_to_template_GM -mul ${g}_JAC_nl ${g}_GM_to_template_GM_mod -odt float" >> fslvbm3a
  fi
done
chmod a+x fslvbm3a
#head -n400 fslvbm3a > fslvbm3a_short 
echo Doing registrations
cat fslvbm3a | xargs -l -P $NPROCS -0 -i -t sh -c '{}'

mkdir $VBM_DIR/stats
cd $VBM_DIR/stats

cat <<stage_preproc2 > fslvbm3b
#!/bin/sh

\$FSLDIR/bin/imcp ../struc/template_GM template_GM

\$FSLDIR/bin/fslmerge -t GM_merg     \`\${FSLDIR}/bin/imglob ../struc/*_GM_to_template_GM.*\`
\$FSLDIR/bin/fslmerge -t GM_mod_merg \`\${FSLDIR}/bin/imglob ../struc/*_GM_to_template_GM_mod.*\`

\$FSLDIR/bin/fslmaths GM_merg -Tmean -thr 0.01 -bin GM_mask -odt char



for i in GM_mod_merg ; do
  for j in 2 3 4 ; do
    \$FSLDIR/bin/fslmaths \$i -s \$j \${i}_s\${j} 
  done
done

stage_preproc2

chmod a+x fslvbm3b

#fslvbm3b_id=`${FSLDIR}/bin/fsl_sub -T 15 -N fslvbm3b -j $fslvbm3a_id ./fslvbm3b`
./fslvbm3b
echo Doing subject concatenation and initial randomise: ID=$fslvbm3b_id

echo "Once this has finished, run randomise with 5000 permutations on the 'best' smoothed 4D GM_mod_merg. We recommend using the -T (TFCE) option. For example:"
echo "randomise -i GM_mod_merg_s3 -o GM_mod_merg_s3 -m GM_mask -d design.mat -t design.con -n 5000 -T -V"


exit 1

    # segment for VBM and get 'myelin' maps
    fslmaths $OUTDIR/T1wtotemplate.nii.gz -div $OUTDIR/T2wtoT1w_template_brain.nii.gz $MASK -mas $OUTDIR/myelin.nii.gz

    fslmaths $OUTDIR/T1wtotemplate.nii.gz -mas $MASK $OUTDIR/T1wtotemplate_brain.nii.gz 
    fslmaths $OUTDIR/T2wtoT1w_template.nii.gz -mas $MASK $OUTDIR/T2wtoT1w_template_brain.nii.gz 
    fast -S 2 -o multipar $OUTDIR/T1wtotemplate.nii.gz $OUTDIR/T2wtoT1w_template_brain.nii.gz 
    fast -o unipar $OUTDIR/T1wtotemplate.nii.gz
    rm $OUTDIR/*pve0* $OUTDIR/*mixel* $OUTDIR/*_seg.nii.gz 


