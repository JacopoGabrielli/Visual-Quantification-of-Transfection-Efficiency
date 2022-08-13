import os
import skimage

def print_file(data, fileout):
	with open(fileout, 'w') as F:
		F.write('{\n')
		F.write('\t"images": [\n')
		for d in data:
			F.write('\t\t{\n')
			F.write('\t\t\t"id": '+ str(d) +',\n')
			F.write('\t\t\t"height": ' + str(data[d]['shape'][0]) + ',\n')
			F.write('\t\t\t"width": ' + str(data[d]['shape'][1]) + ',\n')
			F.write('\t\t\t"file_name": "' + data[d]['name'] + '"\n')
			if d<len(data)-1:
				F.write('\t\t},\n')
			else:
				F.write('\t\t}\n')
		F.write('\t],\n')
		F.write('\t"annotations": [],\n')
		F.write('\t"categories": [\n')
		F.write('\t\t{\n')
		F.write('\t\t\t"id": '+str(1)+',\n')
		F.write('\t\t\t"name": "Cell",\n')
		F.write('\t\t\t"supercategory": "Cell"\n')
		F.write('\t\t}\n')
		F.write('\t]\n')
		F.write('}\n')


folder = '/PATH/TO/CUT/IMAGES/'
file_out = '/PATH/TO/OUTPUT/ANNOTATION/FILE.json'

files = os.listdir(folder)
idx = 0
annotations = {}
for f in files:
	if not f.startswith('.'):
		img = skimage.io.imread(folder+f)
		annotations[idx] = {}
		annotations[idx]['name'] = f
		annotations[idx]['shape'] = img.shape
		idx+=1

print_file(annotations, file_out)
