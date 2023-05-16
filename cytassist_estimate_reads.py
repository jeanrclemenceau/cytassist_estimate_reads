"""
cytassist_estimate_seq_depth.py

Estimate number of reads for a Visium FFPE for CytAssist sample according to its image.
Requires  matplotlib, scikit-image, and pandas.

For help, execute: python cytassist_estimate_seq_depth.py -h

Author: Jean R Clemenceau
Date Created: 02/16/2023
"""

import re
import argparse
import numpy as np
import pandas as pd
import logging as log
from pathlib import Path
import matplotlib.pyplot as plt
from skimage.color import rgb2gray
from skimage.filters import threshold_otsu
from skimage.morphology import remove_small_objects, opening,closing, square

# Log Levels
VERBOSE_LEVEL = [log.ERROR, log.WARNING, log.INFO, log.DEBUG]

#Default trimming parameters
VISIUM_PARAMS = {
    "6.5mm": {
        'TOTAL_SPOTS': 4992,
        'ROI_X':230,
        'ROI_Y':768,
        'SIZE_X':1421,
        'SIZE_Y':1496    
    },
    "11mm": {
        'TOTAL_SPOTS': 14336,
        'ROI_X':280,
        'ROI_Y':300,
        'SIZE_X':2492,
        'SIZE_Y':2454    
    }
}
VISIUM_PARAM_OPTIONS = list(VISIUM_PARAMS.keys())

#Visium read depth default parameters
READS_PER_SPOT = 25000

if __name__ == '__main__':
    # Define command line arguments
    ap = argparse.ArgumentParser()
    ap.add_argument('-i', '--input', default=None, help='Input CytAssist image. Required.')
    ap.add_argument('-o', '--output', default="./", help='Output directory Default: [./]')
    ap.add_argument('-r', '--reads_per_spot', default=READS_PER_SPOT, type=int, help='Desired number of reads per spot. Default: [%d]' % READS_PER_SPOT)
    ap.add_argument('-s', '--slide_type', default=VISIUM_PARAM_OPTIONS[0], type=str, choices=VISIUM_PARAMS, help='Total number of spots available. Default: [%s], Options: [%s]' % (VISIUM_PARAM_OPTIONS[0], VISIUM_PARAM_OPTIONS))
    ap.add_argument('-t', '--threshold', default=None, type=float, help='Thresholding value for foreground mask. Default: [Otsu]' )
    ap.add_argument('-v', '--verbose', action='count', help='Print updates and reports as program executes. Provide the following number of "v" for the corresponding settings: [Default: Error, v: Warning, vv: Info, vvv: Debug]')
    args = vars(ap.parse_args())

    # Validate arguments
    if args['verbose'] is not None:
        if args['verbose'] > len(VERBOSE_LEVEL)-1:
            log.getLogger().setLevel(VERBOSE_LEVEL[-1])
        else:
            log.getLogger().setLevel(VERBOSE_LEVEL[args['verbose']])
    else:
        log.getLogger().setLevel(log.ERROR)

    # Validate Visium parameters
    if args['slide_type'] not in VISIUM_PARAM_OPTIONS:
        raise ValueError("Slide type given is not a valid option.")
    slide_params = VISIUM_PARAMS[args['slide_type']]
    
    # Validate and find input image
    if args['input'] is None:
        raise ValueError("No input path was given.")
    else:
        input_path = Path(args['input'])

        if not input_path.is_file():
            raise FileNotFoundError("Mask image file was not found.")
        elif input_path.suffix not in (".tif",".tiff",".png"):
            raise ValueError("File type is not compatible. Must be TIF or PNG.")
        else:
            # Import image
            log.info("Input image found: %s" % input_path)
            tissue_img = np.array(plt.imread(input_path))

    #Create grayscale image
    tissue_gray = (rgb2gray(tissue_img) * 255).astype(np.uint8)

    #Trim image for capture region only
    log.debug("Trimming image to capture area: [%d : %d, %d : %d]" % 
              (slide_params['ROI_Y'],slide_params['ROI_Y']+slide_params['SIZE_Y'],
               slide_params['ROI_X'],slide_params['ROI_X']+slide_params['SIZE_X']))
    tissue_trimmed = tissue_gray[
        slide_params['ROI_Y']:slide_params['ROI_Y']+slide_params['SIZE_Y'],
        slide_params['ROI_X']:slide_params['ROI_X']+slide_params['SIZE_X']
        ]

    #Generate tissue foreground mask
    log.debug("Generating tissue mask")
    
    if args['threshold'] is None:
        threshold = threshold_otsu(tissue_trimmed)*1.
    else:
        threshold = float(args['threshold'])

    tissue_mask = (tissue_trimmed[:, :] < threshold)
    tissue_mask = remove_small_objects(tissue_mask, 50)
    tissue_mask = closing(tissue_mask, square(5))
    tissue_mask = opening(tissue_mask, square(5))

    #Calculate Sequencing Depth
    foreground_ratio = np.sum(tissue_mask)/tissue_mask.size
    dots_used = np.ceil(foreground_ratio * slide_params['TOTAL_SPOTS'])
    seq_depth = np.ceil(dots_used * args['reads_per_spot'])
    seq_data = {
        'specimen': re.sub('^[A|D]1_','',input_path.stem),
        'filename': input_path.name,
        'slide_type': args['slide_type'],
        'mask_threshold': threshold,
        'foreground_ratio': foreground_ratio,
        'spots_used': dots_used,
        'reads_per_spot': args['reads_per_spot'],
        'total_reads': seq_depth,
    }

    #Report results
    for key, value in seq_data.items():
        key_str = key.replace("_"," ")
        if type(value) == str:
            log.info("  %s:  %s" % (key_str,value))
        else:
            log.info(f"  {key_str}:  {value:,.2f}")
        

    #Prepare output directory
    outpath = Path(args["output"])
    outfile_path = None

    if outpath.is_file():
        if outpath.name.endswith("seq_params.tsv"):
            outfile_path = outpath
        outpath = outpath.parent
    if not outpath.is_dir():
        outpath.mkdir(parents=True)

    #Set filename path
    if outfile_path is None:
        seqparam_files = [i for i in outpath.glob("*seq_params.tsv")]
        if len(seqparam_files) < 1:
            outfile_path = outpath / "seq_params.tsv"
        else:
            outfile_path = seqparam_files[0]

    log.debug("Output directory set to: %s" % outpath)
    log.info("Output file set to: %s" % outfile_path)
    

    #Save results to TSV file via pd.Dataframe
    if outfile_path.is_file():
        pd_seqdata = pd.read_csv(outfile_path, sep="\t")
        pd_seqdata = pd_seqdata.append(seq_data,ignore_index=True)
    else:
        for key, value in seq_data.items():
            seq_data[key] = [value]
        pd_seqdata = pd.DataFrame.from_dict(seq_data, orient="columns")

    log.debug("Exporting sequencing parameters to TSV file: %s" % outfile_path)
    pd_seqdata.to_csv(outfile_path,sep="\t", index=False)

    #Save image mask
    mask_path = outpath / ("%s__tissue_mask.png" % input_path.stem)
    log.debug("Exporting tissue mask to PNG: %s" % mask_path)
    plt.imsave(mask_path, tissue_mask, cmap="Greys_r")

    log.info("Processing Complete.")
