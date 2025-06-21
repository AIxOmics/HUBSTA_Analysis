import copy
import os
import shutil

import tifffile
import cv2
import json
from wsi_split import SplitWSI
import numpy as np
import sys
import glob



def f_ij_16_to_8(img, chunk_size=1000):
    """
    16 bits img to 8 bits

    :param img: (CHANGE) np.array
    :param chunk_size: chunk size (bit)
    :return: np.array
    """

    if img.dtype == 'uint8':
        return img
    dst = np.zeros(img.shape, np.uint8)
    p_max = np.max(img)
    p_min = np.min(img)
    scale = 256.0 / (p_max - p_min + 1)
    for idx in range(img.shape[0] // chunk_size + 1):
        sl = slice(idx * chunk_size, (idx + 1) * chunk_size)
        win_img = copy.deepcopy(img[sl])
        win_img = np.int16(win_img)
        win_img = (win_img & 0xffff)
        win_img = win_img - p_min
        win_img[win_img < 0] = 0
        win_img = win_img * scale + 0.5
        win_img[win_img > 255] = 255
        dst[sl] = np.array(win_img).astype(np.uint8)
    return dst


def feature_dct(coords, id_no):
    g_dct = {
        'type': 'Polygon',
        'coordinates': coords
    }

    p_dct = {
        'objectType': 'annotation',
        'classification': {'name': 'cells', 'color': [255, 0, 0]}
    }

    f_dct = {
        'type': 'Feature',
        'id': id_no,
        'geometry': g_dct,
        'properties': p_dct
    }

    return f_dct


def write_qupath_object(out_file, mask):
    j_dct = {
        'type': 'FeatureCollection',
        'features': []
    }

    cnts, _ = cv2.findContours(mask, cv2.RETR_LIST, cv2.CHAIN_APPROX_TC89_L1)
    for i in range(len(cnts)):
        cnt = cnts[i].tolist()
        cnt = [[val[0] for val in cnt]]

        if len(cnt[0]) < 3:
            continue

        cnt[0].append(cnt[0][0])
        f_dct = feature_dct(cnt, i)
        j_dct['features'].append(f_dct)

    with open("{}".format(out_file), 'w') as jf:
        json.dump(j_dct, jf)


def del_files(filepath):
    """
    删除某一目录下的所有文件或文件夹
    :param filepath: 路径
    :return:
    """
    if os.path.isfile(filepath):
        os.remove(filepath)
    else:
        del_list = os.listdir(filepath)
        for f in del_list:
            file_path = os.path.join(filepath, f)
            if os.path.isfile(file_path):
                os.remove(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        os.rmdir(filepath)
    return


def split_and_save(img_file, mask_file, tar_path, win_size=256):
    name = os.path.split(img_file)[-1]
    name = os.path.splitext(name)[0]

    # img_dir = os.path.join(tar_path, name, "img")
    # json_dir = os.path.join(tar_path, name, "json")

    # if not os.path.exists(json_dir):
    #     os.makedirs(json_dir)
    # if not os.path.exists(img_dir):
    #     os.makedirs(img_dir)

    img = tifffile.imread(img_file)
    img = np.squeeze(img)
    img = f_ij_16_to_8(img)
    mask = tifffile.imread(mask_file)
    mask = np.squeeze(mask)
    mask = np.clip(mask, 0, 1)
    mask = np.uint8(mask * 255)
    
    sp_run = SplitWSI(img, img.shape, 0, 0, False, False, False, np.uint8)
    box_lst, _, _ = sp_run.f_split2run()

    for i in range(len(box_lst)):
        y_begin, y_end, x_begin, x_end = box_lst[i]
        win_mask = mask[y_begin:y_end, x_begin:x_end]
        if not np.sum(win_mask) > 0:
            continue
        # win_img = img[y_begin:y_end, x_begin:x_end]
        write_qupath_object(os.path.join(tar_path, f"{name}.json"), win_mask)
        # tifffile.imwrite(os.path.join(img_dir, f"{name}.tif"), win_img,compression='zlib')



import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--minus_file', type=str, help='Path to the ssDNA file path')
parser.add_argument('--watershed_file', type=str, help='Path to the cell mask file path')
parser.add_argument('--tar_path', type=str, help='Path to the tar file path')

args = parser.parse_args()
minus_file = args.minus_file
watershed_file = args.watershed_file
tar_path = args.tar_path
split_and_save(minus_file, watershed_file, tar_path =tar_path, )