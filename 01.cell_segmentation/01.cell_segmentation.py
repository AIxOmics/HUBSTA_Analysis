# environment: cellbin
from cellbin.modules.cell_segmentation import CellSegmentation
from cellbin.image.augmentation import f_ij_16_to_8
from cellbin.image.augmentation import f_ij_16_to_8_v2
import tifffile
import matplotlib.pyplot as plt
import cv2
import numpy as np
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--img_file', type=str, help='Path to the ssDNA file path')
parser.add_argument('--mask_file', type=str, help='Path to the cell mask file path')

args = parser.parse_args()
img_file = args.img_file
mask_file = args.mask_file

gpu = 0
num_threads = 0
cell_bcdu = CellSegmentation(
    model_path="/data/work/Cellbin/model/cell_segmetation_v3.0.onnx",
    gpu=gpu,
    num_threads=num_threads
)

file_type = "SSDNA"
img = tifffile.imread(img_file)
img = f_ij_16_to_8_v2(img)
# img = cell_bcdu.get_trace(img)

if file_type != "HE" and img.ndim == 3:
    img = cv2.cvtColor(img, cv2.COLOR_RGB2GRAY)

mask = cell_bcdu.run(img)
tifffile.imwrite(mask_file, mask, compression='zlib')