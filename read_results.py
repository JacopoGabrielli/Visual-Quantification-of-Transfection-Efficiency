import os
import pickle
import matplotlib.pyplot as plt
import numpy as np
import skimage.io
import skimage.color

def combine_masks(msks):
    sz = msks[0].mask.shape
    out_bin = np.zeros(sz)
    out_label = np.zeros(sz)
    id = 1
    for m in msks:
        out_bin[m.mask > 0] = 1
        out_label[m.mask > 0] = id
        id += 1
    return out_bin, out_label


folder_images = '/PATH/TO/SCALED/IMAGES/'
folder_masks = '/PATH/TO/MASK/FILES/'
folder_out = '/PATH/WHERE/TO/SAVE/RESULTS/'

files = os.listdir(folder_images)
for f in files:
    if not f.startswith('.'):
        image = skimage.io.imread(folder_images + f)
        try:
            with open(folder_masks + f+'.npy', 'rb') as F:
                masks = pickle.load(F)
            mask, label = combine_masks(masks)
            label_overlay = skimage.color.label2rgb(mask, image, bg_label=0)
            with open(folder_out + f + '.npy', 'wb') as F:
                pickle.dump(mask, F)

        except FileNotFoundError:
            print(f)


