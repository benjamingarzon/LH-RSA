cd /home/xgarzb@GU.GU.SE/Data/LeftHand/Lund1/fmriprep/fmriprep

rm rejected_regs.txt
for i in *html; do    
  LINE=`cat $i | grep rejected | cut -d'.' -f4 |cut -d'/' -f4 `;   
  echo $LINE | tr " " "\n"| cut -d'_' -f1,2,4 --output-delimiter=' ' >> rejected_regs.txt;  
done

# average
# compute NMI

