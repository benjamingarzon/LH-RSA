

Pipeline:

clean after phase 5 and 6

BAD RECONS lue1101 long 2 

# use residuals from GLM (step 5) for connectivity analysis

check that strutctural images are complete
NOT USING HIGH PASS, SHOULD IT BE USED foR ERP?? REPASSAR parametres

repassar motion parameters from fmriprep

think about using fmriprerp 1.4 whe it is ready
improve run_fs pipeline???

xxconversions.sh: should read from file which are the correct sequences
download data from 
should read how many runs to put in fieldmap


fix in line 168
		try:
                    fp.write('\t'.join([str(val) for val in seq]) + '\n')
		except UnicodeEncodeError:
		    print(seq)
		    print("Recoding")
                    fp.write('\t'.join([str(val.encode('utf-8')) if not isinstance(val, (int, float,tuple)) else str(val) for val in seq]) + '\n')

look at timecourses, return to baseline
remove derivatives from fixation ???
organize responses
output within and between values for all subjects and show bars with errors .. extract surface
collapse for different sequences

correct LSS
run with and without PCA, keep it simple, with and without derivatives
fMRI = shorter runs, more trials, random ITI, no locking
fix problem with sessions of different subject
for i in */sub-105*; do mv $i `echo $i | sed 's/sub-105/sub-105.1/g'`; done


xx1. test folder
xx2. use tstat
xx3. add pca

4. larger radius : 20
add midthickness

ridge with multilabel
amplify correlation



xxuse cortex labels
mixture of gaussians and compare IC
do surface based fMRI analysis and use activation sites as mask
xxdo 'prewhitening'
permutation and rsa

are some maps 0?? no, open again with freesurfer




Find centers without spatiotemporal and then do spatiotemporal 

Fix volume searchlight nicely
Run fMRI
Revisar cada pas en les dades noves
Using euclidean and spatiotemporal: spread is high in somatosensory and motor cortex!


Do I need to normalize the features (e. g. to 1 )

look at tutorial

shuffle (within each chunk)? and fit again and see if statistically significant, use parallel toolbox 
xx Plot MDS
Check consistency scores again

Remove outliers

Do spatiotemporal heatmap fig
put in fsaverage to plot better
Try Mahalanobis distance

try accuracy on extract

Run chunks separately and plot

Run the two hemispheres
Add and plot the decoding results
Check the multiclass
Run also volume searchlights
Look at trained vs untrained

xxPCA
xx add derivatives

xxxRun with 5 mm
xxPlot matrices
xxxCan use rest scans also when accuracy was bad
xxxDo maps of within, between and score
xxxGo back to Searchlight and try what you learned, also classification
xxxCheck decoding for spatio-temporal patterns ---> but think about timing issues!
xxRedo surface_searchlight with new thinking about zsxores
xxDo a loop with different vertices and their plots
xxTry multiclass and ridge
xxhow to use multivariate searchlight?

Do I improve sensitivity if I include more data?

xxCompute within vs between measure

At the end run also with GLM

xx Plot time courses

Check decoding as we move far away in time


Talk to Peter: unwarping, dummies, Sync the trigger
Karin: dummies


Backup !!