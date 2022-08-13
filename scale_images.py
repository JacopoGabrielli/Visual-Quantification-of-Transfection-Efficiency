import os
import numpy as np
import skimage.io
import skimage.exposure
import skimage.transform
import skimage.morphology
import matplotlib.pyplot as plt

'''
#scale images

folder_in = ''
folder_out = ''

files = os.listdir(folder_in)
for f in files:
    if not f.startswith('.'):
        image = skimage.io.imread(folder_in + f)
        scaled_image = 255*skimage.transform.rescale(image, 0.1)
        skimage.io.imsave(folder_out + f, scaled_image.astype(np.uint8))

'''
'''
#update .json

file_conf = ''
file_out = ''
new_height = '313'
new_width = '471'
data ={}
with open(file_conf, 'r') as F:
    for rr, r in enumerate(F.readlines()):
        if 'height' in r:
            before = r.split(':')[0]
            data[rr] = before + ': ' + new_height + ',\n'
        elif 'width' in r:
            before = r.split(':')[0]
            data[rr] = before + ': ' + new_width + ',\n'
        else:
            data[rr] = r

K = sorted(list(data.keys()))
with open(file_out, 'w') as F:
    for k in K:
        F.write(data[k])

'''
# increase the contrast
folder_in = '/PATH/TO/SCALED/IMAGES/'
folder_nuclei = '/PATH/TO/NUCLEI/STAINED/IMAGES/'
folder_out = '/PATH/TO/OUTPUT/FOLDER/'
files = os.listdir(folder_in)
for f in files:
    if not f.startswith('.'):
        image = skimage.io.imread(folder_in + f)
        nuclei = skimage.io.imread(folder_nuclei + f.split('C2')[0] +'C3'+f.split('C2')[1])
        nuclei = skimage.transform.rescale(nuclei, 0.1)
        scaled = (image - np.amin(image))/np.amax(image)
        h, bin = skimage.exposure.cumulative_distribution(nuclei)
        for ihh, hh in enumerate(h):
            if hh > 0.95:
                th_bkg = bin[ihh-1]
                break
        corrected = np.copy(scaled)
        #corrected[nuclei > th_bkg] = np.mean(image)
        #corrected = skimage.morphology.dilation(corrected, skimage.morphology.disk(7))
        #percentile = np.percentile(corrected, (0.05, 99.5))
        scaled = (corrected - np.amin(corrected))/np.amax(corrected) #skimage.exposure.rescale_intensity(corrected, in_range=tuple(percentile))
        skimage.io.imsave(folder_out + f, corrected)

