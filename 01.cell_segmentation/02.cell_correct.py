from cellbin.contrib.fast_correct import Fast
from cellbin.modules.cell_labelling import CellLabelling
import tifffile
import cv2
parser = argparse.ArgumentParser()

parser.add_argument('--gem_file', type=str, help='Path to the gem file path')
parser.add_argument('--mask_file', type=str, help='Path to the cell mask file path')
parser.add_argument('--out_path', type=str, help='Path to the output file path')

args = parser.parse_args()
gem_file = args.gem_file
mask_file = args.mask_file
out_path = args.out_path


mask = tifffile.imread(mask_file)
t = CellLabelling(mask = mask, gene_file = gem_file)
fast = Fast(t.mask)
fast.process()
fast_mask = fast.get_mask_fast()
cv2.imwrite(out_path, fast_mask)


'''
from cellbin.contrib.fast_correct import Fast
from cellbin.modules.cell_labelling import CellLabelling
import tifffile
import cv2
import argparse
def exp(gem_file, mask_file, out_path):
    mask = tifffile.imread(mask_file)
    t = CellLabelling(mask = mask, gene_file = gem_file)
    fast = Fast(t.mask, process = 96)
    fast.process()
    fast_mask = fast.get_mask_fast()
    cv2.imwrite(out_path, fast_mask)
'''