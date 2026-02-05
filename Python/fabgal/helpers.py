import numpy as np
from typing import List
from numpy.typing import NDArray
from skimage.restoration import rolling_ball
from skimage.transform import rescale,resize
from skimage.util import img_as_ubyte, img_as_float
from scipy.ndimage import uniform_filter
from pathlib import Path
import logging
from bioio.writers import OmeTiffWriter

#### Log options ####

import logging

logger = logging.getLogger(__name__)


def calculate_bgal(img, bgal_ch: int, bgal_thmin: int, bgal_thmax: int = 254, pxarea: float = None) -> List[float]:
    """
    Calculates the area and intensity of B-gal positive regions in an image.

    This function masks pixels falling within the specific threshold range 
    and calculates the total pixel count, physical area, and intensity sum.

    Args:
        img: A BioImage object containing pixel size metadata.
        parea (float): Pixel area of the image. Defaults to None, where the info will be extracted from the image or, if no info is found, defaults to 1.
        bgal_ch (int): The channel index for the B-gal signal.
        bgal_thmin (int): The minimum pixel intensity threshold.
        bgal_thmax (int, optional): The maximum pixel intensity threshold. 
            Defaults to 254 to exclude burned (saturated) pixels at 255 
            where real values are unknown.

    Returns:
        list: A list containing three values:
            1. Number of positive pixels (NpxPos)
            2. Total number of pixels (NpxTot)
            3. Physical area of positive region (AreaPos)
            4. Physical area of the image (AreaTot)
            5. Pixel area of the image (pxarea)
            6. Raw integrated density of positive pixels (Bgal_RawIntDen)
    """
    ##### Get pixel area info #####

    pxInfo = True

    # If user does not input a pixel area
    if pxarea is None:

        # If there is information in the image and there is no user input, use embedded information.
        if img.physical_pixel_sizes[1] is not None and img.physical_pixel_sizes[2] is not None:
            pxarea = img.physical_pixel_sizes[1] * img.physical_pixel_sizes[2]
            
        # If there is no info about physical pixel size, defaults to 1 to perform calculations.
        else:
            logger.info("\nWARNING: Image does not have pixel physical sizes. Using default of 1 for B-Gal calculations.", end="\r", flush=True)
            pxInfo = False
            pxarea = 1
    
    ##### Extract image data for the specific channel #####

    imgdata = img.get_image_data("YX", C=bgal_ch)
    
    
    ##### Create a boolean mask for pixels within the threshold #####

    mask = (imgdata >= bgal_thmin) & (imgdata <= bgal_thmax)
    
    ##### Calculate stats #####

    NpxPos = np.sum(mask).item()
    NpxTot = imgdata.size

    if pxInfo is True:
        AreaPos = NpxPos * pxarea
        AreaTot = NpxTot * pxarea
    else:
        AreaPos = None
        AreaTot = None
    
    RawIntDen = np.sum(imgdata[mask]).item()
    
    return [NpxPos, NpxTot, AreaPos, AreaTot, pxarea, RawIntDen]


def subtract_background(img, radius : int, scale : int = 8) -> NDArray[np.uint8]:
    """
    Subtracts the background signal of an input image using the rolling ball algorithm.

    This function smooths the image using an uniform filter, estimates background signal via a downscaled 
    rolling ball, and returns the background-subtracted image as a numpy array of uint8 format in order to be saved.

    Args:
        img (np.ndarray): The input image array (pixel data).
        radius (int): Radius to be applied in the rolling_ball algorithm.
        scale (int, optional): Factor by which the image is downscaled to speed up 
            background estimation. Defaults to 8.

    Returns:
        np.ndarray: The background-subtracted image as ubyte (uint8).

    """
    img = uniform_filter(img, size=3)
    img = img_as_float(img)

    #Downscale
    bg = rescale(img,
                 1/scale,
                 order=1,
                 preserve_range=True)
    bg = rolling_ball(bg, radius=radius, nansafe=False)
    bg = resize(bg,
                img.shape,
                order=1,
                preserve_range=True)
    return img_as_ubyte(img-bg)

def generate_biapy_input(img, nuclei_ch: int, apply_sbg: bool, sbg_rad: int, out_path: str):
    """
    Docstring para process_biapy_input
    """
    # Load nuclei channel
    try:
        imgdata = img.get_image_data("YX",C=nuclei_ch)
    except IndexError as e:
        raise IndexError(f"\nERROR: Image does not have nuclei input channel. Please check that nuclei channel parameters are correct.") from e

    # Apply subtract background if needed
    if apply_sbg:
        OmeTiffWriter.save(
            subtract_background(imgdata,
                                radius=sbg_rad), 
                    str(out_path),
                    dim_order="YX")
    else:
        OmeTiffWriter.save(imgdata,
                    str(out_path),
                    dim_order="YX")
        

def load_input(input_folder: str) -> List:
    """
    Docstring para load_input
    """
    # Check input folder
    try:
        input_folder = Path(input_folder)

        if not input_folder.exists():
            raise FileNotFoundError(f"Folder does not exist: {input_folder}")

        if not input_folder.is_dir():
            raise NotADirectoryError(f"Not a directory: {input_folder}")

    except Exception as e:
        logger.info(f"ERROR: Folder validation failed: {e}")

    # Check file extension
    allowed_ext = {".tif", ".tiff"}
    myfiles = [f for f in input_folder.iterdir() if not f.name.startswith(".")] # Read visible files

    for inf in myfiles:
        if inf.is_file() and inf.suffix.lower() not in allowed_ext:
            raise ValueError(f"Invalid file found: {inf.name}. Please ensure all items in input folder are .tif or .tiff.")
    return myfiles