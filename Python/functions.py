# Load libraries

from skimage.restoration import rolling_ball
from skimage.transform import rescale,resize
from skimage.util import img_as_ubyte, img_as_float
from scipy.ndimage import uniform_filter
import numpy as np

from numpy.typing import NDArray
from typing import List

import functions

# Functions

def calculate_bgal(img, bgal_ch: int, bgal_thmin: int, bgal_thmax: int = 254, psize: float = None) -> List[float]:
    """
    Calculates the area and intensity of B-gal positive regions in an image.

    This function masks pixels falling within the specific threshold range 
    and calculates the total pixel count, physical area, and intensity sum.

    Args:
        img: A BioImage object containing pixel size metadata.
        psize (float): Pixel size of the image. Defaults to None, where the info will be extracted from the image or, if no info is found, defaults to 1.
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
            5. Raw integrated density of positive pixels (Bgal_RawIntDen)
    """
    # Calculate area of a single pixel (Y size * X size).
    pixInfo = True
    if psize is not None:
        print(f"Using user input {psize} as physical size.")
        pxarea = psize * psize #  If pixel size information is manually entered, use user input.

    elif img.physical_pixel_sizes[1] is not None and img.physical_pixel_sizes[2] is not None:
         # If there is information in the image and there is no user input, use embedded information.
        print(f"Image pixel sizes are {img.physical_pixel_sizes}")
        pxarea = img.physical_pixel_sizes[1] * img.physical_pixel_sizes[2]

    else:
        print("\nWARNING: Image does not have pixel physical sizes. Using default of 1 for B-Gal calculations.")
        pixInfo = False
        pxarea = 1 #  If there is no info about physical pixel size, defaults to 1 to perform calculations.
    
    # Extract image data for the specific channel
    imgdata = img.get_image_data("YX", C=bgal_ch)
    
    
    # Create a boolean mask for pixels within the threshold
    mask = (imgdata >= bgal_thmin) & (imgdata <= bgal_thmax)
    
    # Calculate stats
    NpxPos = np.sum(mask).item()
    NpxTot = imgdata.size

    if pixInfo is True:
        AreaPos = NpxPos * pxarea
        AreaTot = NpxTot * pxarea
    else:
        AreaPos = None
        AreaTot = None
    
    RawIntDen = np.sum(imgdata[mask]).item()
    
    return [NpxPos, NpxTot, AreaPos, AreaTot, RawIntDen]

## Function to substract background
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