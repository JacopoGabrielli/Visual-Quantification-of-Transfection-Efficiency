import skimage
import os
import shutil
from PIL import Image
Image.MAX_IMAGE_PIXELS = None

def get_idxs(idx, r, c, subimages):
	r_sub = int(r / subimages)
	c_sub = int(c / subimages)
	r_idx = int(idx/subimages)
	c_idx = idx%subimages
	return [r_idx*r_sub, (r_idx+1)*r_sub], [c_idx*c_sub, (c_idx+1)* c_sub]


folder_in = '/PATH/TO/IMAGES/TO/CUT/'
folder_out = '/PATH/WHERE/TO/SAVE/CUT/IMAGES/'
sub_images_per_side = 5
'''
files = os.listdir(folder_in)
for f in files:
	if not f.startswith('.') and 'tif' in f:
		image = skimage.io.imread(folder_in + f)
		[r, c] = image.shape
		idx = 0
		while idx < sub_images_per_side**2:
			r_idxs, c_idxs = get_idxs(idx, r, c, sub_images_per_side)
			new_image = skimage.util.img_as_ubyte(image[r_idxs[0]:r_idxs[1], c_idxs[0]:c_idxs[1]])
			skimage.io.imsave(folder_out + f.split('.tif')[0] + '_id_' + str(idx) + '.jpg', new_image)
			idx+=1

'''


jpg_imgs = os.listdir(folder_out)
for j in jpg_imgs:
	if os.path.isfile(folder_out +j):
		if 'empty' in j:
			shutil.move(folder_out + j, folder_out+'empty/'+ j)
		else:
			if 'C1_' in j:
				shutil.move(folder_out + j, folder_out + 'C1/' + j)
			elif 'C2_' in j:
				shutil.move(folder_out + j, folder_out + 'C2/' + j)
			elif 'C3_' in j:
				shutil.move(folder_out + j, folder_out + 'C3/' + j)
			else:
				print(j)

