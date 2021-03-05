#!/usr/bin/env python

# a flexible version of FSL's overlay
import sys
from nilearn import plotting, datasets, image
import numpy as np
import matplotlib.pyplot as plt
import os
from nilearn.datasets import MNI152_FILE_PATH
from imageio import imwrite
from PIL import Image, ImageFont, ImageDraw

FIGS_DIR='/home/benjamin.garzon/Data/LeftHand/Lund1/figs/structure'
bg_img = MNI152_FILE_PATH

if False:
    bg_img = sys.argv[1] 
    img = sys.argv[2]
    sl = int(sys.argv[3])
    outpng = sys.argv[4]
    lthresh = float(sys.argv[5])
    useatlas = sys.argv[6]
    annotate = sys.argv[7]
    alpha = sys.argv[8]
    alpha = 0.9 if alpha == '' else float(alpha)

else:
    os.chdir('/home/benjamin.garzon/Data/LeftHand/Lund1/cat12/tests/reliability')
    img = "ICC.nii.gz"
    sl = int(7)
    outpng = os.path.join(FIGS_DIR, "reliability-VBM")
    lthresh = float(4)
    mask = "/home/benjamin.garzon/Data/LeftHand/Lund1/cat12/mask.nii.gz"
    annotate = "0"
    alpha = 0.5
    
img = image.load_img(img)
(x, y, z) = image.coord_transform(0, 0, sl, img.affine)
print(sl, x, y, z)

mask_img = image.load_img(mask)

imgs = []
for mode in ['x', 'y', 'z']:
  fig = plt.figure()

  display = plotting.plot_stat_map(img,
    bg_img=bg_img,
    colorbar=False,
    black_bg=True,
    display_mode = mode,
    cut_coords = sl,
    annotate = False,
    draw_cross = False,
    alpha = alpha, 
    figure = fig, 
    symmetric_cbar = True,
    vmax = 1,
    threshold = 0.5,
    cmap = 'cold_hot')
    
  display.add_contours(mask_img,
    filled = False,
    cmap = 'gray', linewidths = 0.5)

  display.savefig(outpng + '-' + mode, dpi = 1000)
  img2 = Image.open(outpng +  '-' + mode + '.png')
  h = img2.height 
  w = img2.width
  rh = 0.35
  imgs.append(img2.crop((0, int(h*rh), w, int(h*(1-rh)))))
  
display.close()
full_img = np.concatenate(imgs, axis = 0)
imwrite(outpng + '.png', full_img)

