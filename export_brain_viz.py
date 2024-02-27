# Joshua Gobé 2023: generate .txt out of .raw files
# Alex Andrix 2023: adaptation for sparser sampling
# Alex Andrix 2024: generate .txt out of .nii mask files with nibabel
# Alex Andrix 2024: generate .txt out of .trk files
# Alex Andrix 2024: generate .txt out of .tif files
#
# Note 1: our custom filetypes are prefixed by 'andrix-' solely to remind that it is not a standard format from a library but some attempt to classify source data formats. Can be improved.
#
# Note 2: global variable is_file_blank is not reinitialized after a file is generated, this means you should only generate 1 file at a time by commenting all other calls. There probably is a better way to do it.

import numpy as np
import nibabel as nib
import math
import random
from PIL import Image


is_file_blank = True

# Open method for .raw files, outputs the volume as a Numpy array
# DTYPE is BYTE DEPTH + ORDER
def open_raw (path, dimx, dimy, dimz, dt):
#    dim_img = (dimz,dimy,dimx)
    dim_img = (dimx,dimy,dimz)
    data_type = np.dtype(dt)
    image = np.fromfile(path, dtype = data_type) #dtype 16bits unsigned
    image = np.reshape(image, dim_img)
    return image

# Open method for .nii mask files
def open_mask (path):
    mask = nib.load(path)
    mask_data = mask.get_fdata()
    print('NII - Mask array shape is ' + str(mask_data.shape))
    return mask_data
    
# Open method for .trk tractogram files
def open_trk (path):
    print('TRK - File format is ' + str(nib.streamlines.detect_format(path)))
    print('TRK - Loading full .trk file...')
    tractogram_file = nib.streamlines.load(path)
    streamlines = tractogram_file.streamlines
    print('TRK - Tractogram has ' + str(streamlines.total_nb_rows) + ' rows.')
    print('TRK - Tractogram common shape is ' + str(streamlines.common_shape) + '.')
    return streamlines
    
# Write a line with (x, y, z) cartesian information and extra 'value' free field
# is_file_blank permet de s'affranchir automatiquement du saut de ligne généré en fin de fichier
def write_line (cx, cy, cz, value, file):
    global is_file_blank
    
    if not(is_file_blank):
        file.write("\n")
    else:
        is_file_blank = False
    
    end = (" " + str(value)) if (value != "") else ""
    line = str(cx) + " " + str(cy) + " " + str(cz) + end
    file.write(line)
    
    return
    
# Write a line with segment information, i.e. start and end point + some extra values
def write_segment (x1, y1, z1, x2, y2, z2, values, file):
    global is_file_blank
    
    if not(is_file_blank):
        file.write("\n")
    else:
        is_file_blank = False
    
    line = str(x1) + " " + str(y1) + " " + str(z1) + " " + str(x2) + " " + str(y2) + " " + str(z2) + " " + str(values)
    file.write(line)
    
    return

# Main function for .raw and .nii file generation
def main (path, dimx, dimy, dimz, dn, datatype, filetype, filename):
    # path: source file path
    # dimx: data X dimension
    # dimy: data Y dimension
    # dimz: data Z dimension
    # dx: sample every dn on X, Y and Z (takes 1 point in dn³)
    # datatype: data type for raw file reading (https://numpy.org/doc/stable/reference/generated/numpy.fromfile.html) (https://numpy.org/doc/stable/reference/arrays.dtypes.html)
    # filetype: can be 'andrix-raw', 'andrix-mask' or 'andrix-mask-dir' (I add andrix just to remember it's not a native param)
    # filename: output file path
    
    if filetype == 'andrix-raw' or filetype == 'andrix-only-coords':
        vol = open_raw(path, dimx, dimy, dimz, datatype)
    elif filetype == 'andrix-mask' or filetype == 'andrix-mask-dir':
        vol = open_mask(path)
    
    f = open(filename, "w")
    x = 0
    dx = dn
    dy = dn
    dz = dn
    while x <= dimx - dx:
        print("Whiling x! %d/max" % x)
        y = 0
        while y <= dimy - dy:
            z = 0
            while z <= dimz - dz:
                v = vol[x, y, z]
                
                if filetype == 'andrix-raw':
                    if v >= threshold:
                        write_line(x, y, z, v, f)
                        
                elif filetype == 'andrix-only-coords':
                    if v >= threshold:
                        write_line(x, y, z, "", f);
                        
                elif filetype == 'andrix-mask':
                    vv = math.trunc(v[0])
                    if vv >= threshold:
                        write_line(x, y, z, vv, f)
                        
                elif filetype == 'andrix-mask-dir':
                    dir0 = math.trunc(v[0])
                    dir1 = math.trunc(v[1])
                    dir2 = math.trunc(v[2])
                    if dir0 != 0 and dir1 != 0 and dir2 != 0:
                        vv = str(dir0) + ' ' + str(dir1) + ' ' + str(dir2)
                        write_line(x, y, z, vv, f)
                        
                z += dz
            y += dy
        x += dx
    return

# Main function for .trk files
def parse_trk (path, file_name_token, dn, threshold):
    trk = open_trk(path)
    rows = trk.total_nb_rows
    scanned_fiber_count = 0
    retained_fiber_count = 0
    node_count = 0
    file_to_write = open('output/aa-fibers-' + file_name_token + '-522x448x400-1-in-' + str(dn) + '-thresh-' + str(threshold) + '.txt', "w")
    
    print('Writing fibers longer than ' + str(threshold) + ' nodes.')
    for fiber in trk[::dn]:

        if scanned_fiber_count % 10000 == 0:
            print('Scanning fiber ' + str(scanned_fiber_count) + '...')
        
        l = len(fiber)
        
        if l >= threshold:
            for n in range(0, l - 1):
                node = fiber[n]
                next_node = fiber[n + 1]
                x1 = np.trunc(10 * node[0]) / 10
                y1 = np.trunc(10 * node[1]) / 10
                z1 = np.trunc(10 * node[2]) / 10
                x2 = np.trunc(10 * next_node[0]) / 10
                y2 = np.trunc(10 * next_node[1]) / 10
                z2 = np.trunc(10 * next_node[2]) / 10
                
                values = str(scanned_fiber_count) + " " + str(l)
                write_segment(x1, y1, z1, x2, y2, z2, values, file_to_write)
                node_count += 1
                
            retained_fiber_count += 1
            
        scanned_fiber_count += 1

    print(str(retained_fiber_count) + " fibers longer than " + str(threshold) + " nodes were written. The total number of written nodes is " + str(node_count) + ".")

    return

# Main function for .tiff files
def parse_tiff (path, file_name_token, dn):

    image = Image.open(path)
    first_image_array = np.array(image)
    
    print('TIFF - Multi-frame TIFF Image has ' + str(image.n_frames) + ' frames with dimensions ' + str(first_image_array.shape) + '.')

    w, h, d = first_image_array.shape[0], first_image_array.shape[1], image.n_frames

    new_arr = np.zeros((d // dn, w // dn, h // dn), np.int8)
    print('TIFF - Destination cuboid shape is ' + str(new_arr.shape))

    file_to_write = open('output/aa-' + file_name_token + '-' + str(w) + 'x' + str(h) + 'x' + str(d) + '-1-in-' + str(dn) + '.txt', "w")

    for i in range(0, d - 1, dn):

        z = i // dn
        
        if z % 10 == 0:
            print('Scanning image layer ' + str(i) + '...')
        
        image.seek(i)
        img = np.array(image)
        
        for j in range(0, w - 1, dn):
            x = j // dn
            for k in range(0, h - 1, dn):
                y = k // dn
                
                if img[j][k] != 0:
                    new_arr[z][x][y] = img[j][k]
                    write_line(x, y, z, "", file_to_write)
    
    return


# Part of mouse brain // .raw
# dn = 3
# threshold = 0
# datatype = '>u2'
# main('./data/xpct_tracto_data/appps1-mouse1_522x448x400.raw', 400, 448, 522, dn, datatype, 'andrix-raw', 'output/aa-brain-part-cuboid-522x448x400-1-in-' + str(dn) + '-thresh-' + str(threshold) + '.txt')

# Mouse brain // .raw
# dn = 10
# threshold = 1
# datatype = '<u2'
# main('./data/mouse_brain_LPC1_240915_ethanol_pag--deringed--rotated__1179x853x1513.raw', 1513, 853, 1179, dn, datatype, 'andrix-raw', 'output/aa-brain-full-cuboid-1179x853x1513-1-in-' + str(dn) + '-thresh-' + str(threshold) + '.txt')

# Mask with mouse brain envelope (edges) // .raw
# dn = 1
# threshold = 2
# datatype = 'u1'
# filetype = 'andrix-raw'
# main('./data/mouse_brain_ctrl_FA_eth_9_pag_8bit_halfdim_edges_760x755x545.raw', 545, 755, 760, dn, datatype, filetype, 'output/aa-brain-edges-760x755x545-1-in-' + str(dn) + '-thresh-' + str(threshold) + '.txt')
# filetype = 'andrix-only-coords'
# main('./data/mouse_brain_ctrl_FA_eth_9_pag_8bit_halfdim_edges_cleaned_760x755x545.raw', 545, 755, 760, dn, datatype, filetype, 'output/aa-brain-edges-760x755x545-1-in-' + str(dn) + '-thresh-' + str(threshold) + '.txt')

# Mask with white matter // .raw
# dn = 3
# threshold = 1
# datatype = 'u1'
# main('./data/mouse_brain_ctrl_FA_eth_9_pag_8bit_halfdim_wm_erroded_760x755x545.raw', 545, 755, 760, dn, datatype, 'andrix-raw', 'output/aa-white-matter-760x755x545-1-in-' + str(dn) + '-thresh-' + str(threshold) + '.txt')

# Mask with antcom (commissure antérieure) // .raw
# dn = 1
# threshold = 1
# datatype = 'u1'
# filetype = 'andrix-only-coords'
# main('./data/mouse_brain_ctrl_FA_eth_9_pag_8bit_halfdim_antcommask_760x755x545.raw', 545, 755, 760, dn, datatype, filetype, 'output/aa-antcom-760x755x545-1-in-' + str(dn) + '-thresh-' + str(threshold) + '.txt')

# Mask with corpus callosum (corps calleux) // .raw
dn = 2
threshold = 1
datatype = 'u1'
filetype = 'andrix-only-coords'
main('./data/mouse_brain_ctrl_FA_eth_9_pag_8bit_halfdim_corpuscallosummask_760x755x545.raw', 545, 755, 760, dn, datatype, filetype, 'output/aa-corpus-760x755x545-1-in-' + str(dn) + '-thresh-' + str(threshold) + '.txt')

# Mask with intensities // .nii
# dn = 2
# threshold = 1 # Mask has 255s or 0s
# main('./data/xpct_tracto_data/appps1-mouse1_mask.nii', 522, 448, 400, dn, '', 'andrix-mask', 'output/aa-mask-522x448x400-1-in-' + str(dn) + '-thresh-' + str(threshold) + '.txt')

# Mask with directions // .nii
# main('./data/xpct_tracto_data/masked_tensor.nii', 522, 448, 400, dn, '', 'andrix-mask-dir', 'output/aa-mask-dir-522x448x400-1-in-' + str(dn) + '-thresh-' + str(threshold) + '.txt')

# Tractography with fibers // .trk
# dn = 5
# threshold = 0
# parse_trk('./data/sample.trk', 'sample', dn, 0)

# threshold = 7
# dn = 30
# parse_trk('./data/LPC2_demyel.trk', 'demyel', dn, threshold)
# dn = 15
# parse_trk('./data/LPC2_myel.trk', 'myel', dn, threshold)

# Blood vessels // .tif
# dn = 2
# parse_tiff('./data/mouse_brain_ctrl_vaisseaux.tif', 'vessels', dn)

# Blood vessels skeleton // .tif
# dn = 2
# parse_tiff('./data/mouse_brain_ctrl_squelette.tif', 'vessels-skel', dn)
