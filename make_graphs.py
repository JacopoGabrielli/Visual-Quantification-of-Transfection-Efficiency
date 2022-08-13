import os
import pickle
import skimage
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def get_results(data, intensity):
	out = {}
	out_red = {}
	out_area = {}
	for d in data:
		label_plasmid = skimage.measure.label(data[d]['plasmid'])
		out[d] = []
		out_red[d] = []
		out_area[d] = []
		for i in range(1, np.amax(label_plasmid)):
			pixels = np.where(label_plasmid == i)
			out_red[d].append(np.mean(intensity[d][pixels]))
			out_area[d].append(len(pixels[0]))
			in_cell = data[d]['cell'][pixels]
			fraction_true_cell = np.sum(in_cell)/len(in_cell)
			in_nucleus = data[d]['nuclei'][pixels]
			fraction_true_nucleus = np.sum(in_nucleus)/len(in_nucleus)
			if fraction_true_nucleus > 0.0:
				out[d].append(2)
			else:
				if fraction_true_cell > 0.0:
					out[d].append(1)
				else:
					out[d].append(0)
	return out, out_red, out_area

def sort_results(results, intensity, area):
	out ={1: [0, 0, 0], 2: [0, 0, 0], 4: [0, 0, 0]}
	out_intensity = {1:{'O':[], 'C': [], 'N':[]}, 2: {'O':[], 'C':[], 'N':[]}, 4:{'O':[], 'C':[], 'N':[]}}
	out_area = {1:{'O':[], 'C': [], 'N':[]}, 2: {'O':[], 'C':[], 'N':[]}, 4:{'O':[], 'C':[], 'N':[]}}

	for r in results:
		if r.split('HC')[0][-1] =='1':
			for irr, rr in enumerate(results[r]):
				out[1][rr] += 1
				if rr == 0:
					out_intensity[1]['O'].append(intensity[r][irr])
					out_area[1]['O'].append(area[r][irr])
				elif rr == 1:
					out_intensity[1]['C'].append(intensity[r][irr])
					out_area[1]['C'].append(area[r][irr])
				else:
					out_intensity[1]['N'].append(intensity[r][irr])
					out_area[1]['N'].append(area[r][irr])
		elif r.split('HC')[0][-1] == '2':
			for irr, rr in enumerate(results[r]):
				out[2][rr] += 1
				if rr == 0:
					out_intensity[2]['O'].append(intensity[r][irr])
					out_area[2]['O'].append(area[r][irr])
				elif rr == 1:
					out_intensity[2]['C'].append(intensity[r][irr])
					out_area[2]['C'].append(area[r][irr])
				else:
					out_intensity[2]['N'].append(intensity[r][irr])
					out_area[2]['N'].append(area[r][irr])
		elif r.split('HC')[0][-1] == '4':
			for irr, rr in enumerate(results[r]):
				out[4][rr] += 1
				if rr == 0:
					out_intensity[4]['O'].append(intensity[r][irr])
					out_area[4]['O'].append(area[r][irr])
				elif rr == 1:
					out_intensity[4]['C'].append(intensity[r][irr])
					out_area[4]['C'].append(area[r][irr])
				else:
					out_intensity[4]['N'].append(intensity[r][irr])
					out_area[4]['N'].append(area[r][irr])
		else:
			raise ValueError('unrecognized time point')
	return out, out_intensity, out_area

def print_results(results, data, folder_im, PDF):
	images = os.listdir(folder_im)
	for i in images:
		if not i.startswith('.'):
			bf_image = skimage.io.imread(folder_im + i)
			image_scaled = skimage.transform.rescale(bf_image, 0.1)
			combined_image = make_combined_image(data[i.split('.jpg')[0]])
			image_results = get_image_results(results[i.split('.jpg')[0]])
			to_print = skimage.color.label2rgb(combined_image, image=image_scaled, colors=[(0,1,0), (0,0,1), (1,0,0)])
			f,ax = plt.subplots()
			ax.imshow(to_print)
			ax.set_title(i)
			PDF.savefig()

def get_image_results(res):
	out = [0, 0, 0]
	for r in res:
		out[r] += 1
	s = sum(out)
	if s == 0:
		s = 0.1
	return [x/s for x in out]

def make_combined_image(dt):
	mask = np.zeros_like(dt['plasmid'])
	mask[dt['cell']] = 1
	mask[dt['nuclei'].astype(bool)] = 2
	mask[dt['plasmid'].astype(bool)] = 3
	return mask

def print_sorted(sorted, sorted_intensity, sorted_area, PDF):
	f, ax = plt.subplots()
	out = [sorted[1][0], sorted[2][0], sorted[4][0]]
	cell = [sorted[1][1], sorted[2][1], sorted[4][1]]
	nuclei = [sorted[1][2], sorted[2][2], sorted[4][2]]
	time = [1, 2, 4]
	ax.plot(time, out, 'ok', label='out of cell')
	ax.plot(time, cell, 'og', label='in cell')
	ax.plot(time, nuclei, 'ob', label='in nucleus')
	ax.set_title('number of plasmid "clumps"')
	ax.legend()
	PDF.savefig()
	f, ax = plt.subplots()
	ax.plot(time, cell, 'og', label='in cell')
	ax.plot(time, nuclei, 'ob', label='in nucleus')
	ax.set_title('number of plasmid "clumps"')
	ax.legend()
	PDF.savefig()
	f, ax = plt.subplots(3, 1)
	f2, ax2 = plt.subplots(3, 1)
	for tt, t in enumerate(sorted_intensity):
		data_intensity = [sorted_intensity[t]['O'], sorted_intensity[t]['C'], sorted_intensity[t]['N']]
		data_area = [sorted_area[t]['O'], sorted_area[t]['C'], sorted_area[t]['N']]
		ax[tt].boxplot(data_intensity)
		ax[tt].set_xticklabels(['out', 'cell', 'nucleus'])
		ax[tt].set_ylabel(str(t))
		ax2[tt].boxplot(data_area)
		ax2[tt].set_xticklabels(['out', 'cell', 'nucleus'])
		ax2[tt].set_ylabel(str(t))
	ax[0].set_title('plasmid signal intensity')
	ax2[0].set_title('plasmid "clumps" area')
	PDF.savefig(f)
	PDF.savefig(f2)


folder = '/PATH/TO/SEGMENTED/IMAGES/'
folder_imgs = '/PATH/TO/CUT/BRIGHTFIELD/IMAGES/'
folder_plasmid = '/PATH/TO/CUT/PLASMID/STAINED/IMAGES/'
PDF = PdfPages('/PATH/TO/OUTPUT/FILE.pdf')
files = os.listdir(folder)
data = {}
red_intensity = {}
for f in files:
	if not f.startswith('.'):
		with open(folder + f, 'rb') as F:
			data[f.split('.jpg')[0]] = pickle.load(F)
		name_red = f.split('C2')[0] + 'C1' + f.split('C2')[1].split('.npy')[0]
		red_intensity[f.split('.jpg')[0]] = skimage.io.imread(folder_plasmid + name_red)
results, intensities, areas = get_results(data, red_intensity)
sorted, sorted_intensity, sorted_area = sort_results(results, intensities, areas)
print_sorted(sorted, sorted_intensity, sorted_area, PDF)
print_results(results, data, folder_imgs, PDF)
data_out = {'masks': data, 'sorted_nclumps': sorted, 'sorted_intensities': sorted_intensity, 'sorted_area': sorted_area}
with open('/Volumes/Surtsey/SBImperial/Jacopo_images/summary_results.pkl', 'wb') as F:
	pickle.dump(data_out, F)
PDF.close()
