import cv2
import numpy as np
import matplotlib.pyplot as plt
import dill
from skimage.segmentation import clear_border

def remove_small_objects(image, areas, th1, th2):
    out = np.zeros_like(image)
    for ia, a in enumerate(areas):
        if ia > 0: # the first element is the background
            if th1 < a < th2:
                pixels = np.where(image == ia)
                out[pixels] = 1
    return out

def mask_image(image, mask):
    out = np.copy(image)
    [zeros_x, zeros_y] = np.where(mask == 0)
    out[zeros_x, zeros_y] = 0
    return out

def separate_images(image_bright, n_ims = 5, border = 70):
    [w, h] = image_bright.shape
    h_width = int(h/n_ims)
    w_width = int(w/n_ims)
    out = {}
    out_coordinates = {}
    for i in range(n_ims**2):
        c = i % n_ims
        r = int(i/n_ims)
        x_begin = w_width*c + border
        y_begin = h_width*r + border
        x_range =[x_begin, x_begin+w_width]
        y_range = [y_begin, y_begin+h_width]
        out[i] = image_bright[x_range[0]:x_range[1], y_range[0]:y_range[1]]
        out_coordinates[i] = {'x': x_range, 'y': y_range}
    return out, out_coordinates

def combine_images(input, dims, coordinates, n_ims = 5, border = 70):
    out = np.zeros(dims)
    for i in input:
        out[coordinates[i]['x'][0]:coordinates[i]['x'][1], coordinates[i]['y'][0]:coordinates[i]['y'][1]] = input[i]
    return out


def remove_outside_of_cells(nuclei, bright):
    out = np.zeros_like(nuclei)
    ones_nuclei = np.where(nuclei > 0)
    for o in range(len(ones_nuclei[0])):
        if bright[ones_nuclei[0][o], ones_nuclei[1][o]] == 1:
            out[ones_nuclei[0][o], ones_nuclei[1][o]] = 1
    return out



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    file_name = '/Volumes/Surtsey/SBImperial/Jacopo_images/exported/T3_1_Hoescht_Cy3_15_1H'
    file_name_red = file_name + 'C1.tif'
    file_name_blue = file_name + 'C3.tif'
    file_name_bright = file_name + 'C2.tif'
    file_out = file_name+'_output.pickle'

    image_bright = cv2.imread(file_name_bright)
    image_bright = image_bright[:, :, 0]

    images, coordinates = separate_images(image_bright)
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (35, 35))
    kernel2 = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (51, 51))
    bkg = {}
    blobs = {}
    roughly_segmented = {}
    binary = {}
    post_process = {}
    total = {}
    label = {}
    stats = {}
    centroid = {}

    for i in images:
        print(i)
        bkg[i] = cv2.morphologyEx(images[i], cv2.MORPH_BLACKHAT, kernel)
        blobs[i] = cv2.adaptiveThreshold(bkg[i], 1, cv2.ADAPTIVE_THRESH_MEAN_C, cv2.THRESH_BINARY_INV, 261, 4)
        roughly_segmented[i] = mask_image(images[i], blobs[i])
        [total_temp, label_temp, stats_temp, centroid_temp] = cv2.connectedComponentsWithStats(roughly_segmented[i], connectivity=8)
        binary[i] = remove_small_objects(label_temp, stats_temp[:,-1], 50, 10000)
        post_process[i] = cv2.morphologyEx(np.uint8(binary[i]), cv2.MORPH_CLOSE, kernel2)
        [total[i], label[i], stats[i], centroid[i]] = cv2.connectedComponentsWithStats(roughly_segmented[i], connectivity=8)
    final_image_bright = combine_images(post_process, image_bright.shape, coordinates)
    f, ax = plt.subplots(1, 2, sharex=True, sharey=True)
    ax[0].imshow(image_bright, cmap='gray')
    ax[0].set_title('original image: brightfield')
    ax[1].imshow(final_image_bright, cmap='gray')
    ax[1].set_title('segmented image: brightfield')
    image_red = cv2.imread(file_name_red)
    image_blue = cv2.imread(file_name_blue)

    s_cy3_dna = image_red[:, :, 0]
    s_nuclei = image_blue[:, :, 0]

    thresh_nuclei = cv2.adaptiveThreshold(s_nuclei, 1, cv2.ADAPTIVE_THRESH_MEAN_C, cv2.THRESH_BINARY_INV, 31, 2)#cv2.threshold(s_nuclei, cv2.THRESH_OTSU, 255, cv2.THRESH_BINARY_INV)
    thresh_cy3 = cv2.adaptiveThreshold(s_cy3_dna, 1, cv2.ADAPTIVE_THRESH_MEAN_C, cv2.THRESH_BINARY_INV, 31, 4)
    f2, ax2 = plt.subplots(2,1)
    kernel_blue_red = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (15, 15))
    closing_nuclei = cv2.morphologyEx(thresh_nuclei, cv2.MORPH_CLOSE, kernel_blue_red, iterations=2)
    closing_nuclei = clear_border(closing_nuclei)  # Remove edge touching grains
    mask_nuclei = remove_outside_of_cells(closing_nuclei, final_image_bright)
    ax2[0].imshow(mask_nuclei, cmap='gray')  # This is our image to be segmented further using watershed
    ax2[1].imshow(thresh_cy3, cmap='gray')

    reds = np.where(thresh_cy3 > 0)
    in_cell = []
    in_nucleus = []
    out_cell = []
    for r in range(len(reds[0])):
        if mask_nuclei[reds[0][r], reds[1][r]] > 0:
            in_nucleus.append(image_red[reds[0][r], reds[1][r]])
        else:
            if final_image_bright[reds[0][r], reds[1][r]] > 0:
                in_cell.append(image_red[reds[0][r], reds[1][r]])
            else:
                out_cell.append(image_red[reds[0][r], reds[1][r]])
    print('plasmids in nucleus:' + str(len(in_nucleus)/(len(in_nucleus) + len(in_cell) + len(out_cell))))
    print('plasmids in cells, out of nucleus:' + str(len(in_cell)/(len(in_nucleus) + len(in_cell) + len(out_cell))))
    print('plasmids out of cells:' + str(len(out_cell)/(len(in_nucleus) + len(in_cell) + len(out_cell))))

    dill.dump_session(file_out)
    plt.show()