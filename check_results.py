import os
import pickle

import matplotlib.pyplot as plt
import numpy as np
import skimage

def filter_area(out, label, th_area_min, th_area_max):
	for i in range(np.amax(label)):
		pixels = np.where(label == i)
		area = len(pixels[0])
		if area < th_area_min or area > th_area_max:
			out[pixels] = 0
	return out

def segment_nuclei(img):
	out = np.zeros_like(img)
	scaled = np.copy(img)
	h, bin = skimage.exposure.cumulative_distribution(img)
	for ihh, hh in enumerate(h):
		if hh > 0.95:
			th_bkg = bin[ihh-1]
			break
	scaled[scaled < th_bkg] = np.mean(img)
	#f, ax = skimage.filters.try_all_threshold(img)
	th = skimage.filters.threshold_triangle(scaled)
	out[scaled > th] = 1
	label = skimage.measure.label(out)
	out = filter_area(out, label, 100, 1000)
	return out

def segment_plasmid(img):
	#f, ax = skimage.filters.try_all_threshold(img)
	#plt.show()
	th = skimage.filters.threshold_yen(img)
	out = np.zeros_like(img)
	out[img > th] = 1
	return out

folder = '/PATH/TO/BRIGHTFIELD/IMAGES/'
folder_nuclei = '/PATH/TO/NUCLEI/STAINED/IMAGES/'
folder_plasmid = '/PATH/TO/PLASMID/STAINED/IMAGES/'
folder_segment = '/PATH/TO/SEGMENTED/IMAGES/''
folder_out = 'PATH/WHERE/TO/SAVE/RESULTS/'
images = os.listdir(folder)
for i in images:
	if not i.startswith('.'):
		image_raw = skimage.io.imread(folder + i)
		image_nuclei = skimage.io.imread(folder_nuclei + i.split('C2')[0] + 'C3' + i.split('C2')[1])
		image_plasmid = skimage.io.imread(folder_plasmid + i.split('C2')[0] + 'C1' + i.split('C2')[1])

		image_scaled = skimage.transform.rescale(image_raw, 0.1)
		image_nuclei_scaled = skimage.transform.rescale(image_nuclei, 0.1)

		image_plasmid_scaled = skimage.transform.rescale(image_plasmid, 0.1)
		with open(folder_segment + i + '.npy', 'rb') as F:
			segmentation_cell = pickle.load(F).astype(bool)
		segmentation_nuclei = segment_nuclei(image_nuclei_scaled)
		segmentation_plasmid = segment_plasmid(image_plasmid_scaled)
		label_plasmid = skimage.measure.label(segmentation_plasmid)
		segmentation_plasmid.astype(bool)

		plasmid_in_nuclei = np.logical_and(segmentation_plasmid, segmentation_nuclei)
		label_pin = skimage.measure.label(plasmid_in_nuclei)
		plasmid_in_cell = np.logical_and(segmentation_plasmid, segmentation_cell)
		label_pic = skimage.measure.label(plasmid_in_cell)
		plasmids_out = np.amax(label_plasmid) - np.amax(plasmid_in_nuclei) - np.amax(plasmid_in_cell)
		to_save = {'nuclei': segmentation_nuclei, 'cell': segmentation_cell, 'plasmid': segmentation_plasmid}
		with open(folder_out+i+'.npy', 'wb') as F:
			pickle.dump(to_save, F)


