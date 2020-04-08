#!/bin/sh
# compress folders to save space

WD=/home/share/MotorSkill/dicoms_wave3/
cd $WD

for SUBJECT in lue*; do
  cd $SUBJECT
  FOLDERS=*_??.??.??
  echo $FOLDERS
  for f in $FOLDERS; do
    if [ ! -e ${f}.tar.gz ]; then
      echo "Compressing $SUBJECT/${f}"
      tar -czf ${f}.tar.gz ${f} && rm -R $f  
    fi
  done
  cd ..
done
