import numpy as np
from typing import List
from numpy.typing import NDArray
from skimage.restoration import rolling_ball
from skimage.transform import rescale,resize

from scipy.ndimage import uniform_filter
from pathlib import Path
import logging
from bioio.writers import OmeTiffWriter

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Safely check if we are in Jupyter to handle screen clearing
try:
    from IPython.display import clear_output
    is_jupyter = True
except ImportError:
    is_jupyter = False

#### Log options ####

import logging

logger = logging.getLogger(__name__)


def calculate_bgal(img, bgal_ch: int, bgal_thmin: int, bgal_thmax: int = 254, pxarea: float = None, pxunit: str = None) -> List[float]:
    """
    Calculates the area and intensity of B-gal positive regions in an image.

    This function masks pixels falling within the specific threshold range 
    and calculates the total pixel count, physical area, and intensity sum.

    Args:
        img: A BioImage object containing pixel size metadata.
        bgal_ch (int): The channel index for the B-gal signal.
        bgal_thmin (int): The minimum pixel intensity threshold.
        bgal_thmax (int, optional): The maximum pixel intensity threshold. 
            Defaults to 254 to exclude burned (saturated) pixels at 255 
            where real values are unknown.
        pxarea (float): Pixel area of the image. Defaults to None, where the info will be extracted from the image or, if no info is found, defaults to 1.
        pxunit (str): Units of pixel length/width.


    Returns:
        list: A list containing three values:
            1. Number of positive pixels (NpxPos)
            2. Total number of pixels (NpxTot)
            3. Physical area of positive region (AreaPos)
            4. Physical area of the image (AreaTot)
            5. Pixel area of the image (pxarea)
            6. Pixel area units (pxunit)
            7. Raw integrated density of positive pixels (Bgal_RawIntDen)
            8. Mean intensity of signal in the image (Mean_Intens)
    """
    ##### Get pixel area info #####

    pxInfo = True

    # If user does not input a pixel area
    if pxarea is None:

        # If there is information in the image and there is no user input, use embedded information.
        if img.physical_pixel_sizes[1] is not None and img.physical_pixel_sizes[2] is not None:
            pxarea = img.physical_pixel_sizes[1] * img.physical_pixel_sizes[2]
            
        # If there is no info about physical pixel size, create a warning.
        else:
            print("\nWARNING: Image does not have pixel physical sizes. No calculations per area will be performed.", end="\r", flush=True)
            pxInfo = False
    
    ##### Pixel area conversion if not um2 #####

    pxareaunit = "µm^2"
    area_conv = 1 # Conversion factor for pixel area (Defaults to 1 for µm^2)

    if pxInfo and pxunit != "µm":
        if pxunit == "mm":
            pxarea = pxarea * 10^6
        elif pxunit == "cm":
            pxarea = pxarea * 10^8
        else:
            pxareaunit = None
    
    ##### Extract image data for the specific channel #####

    try:
        imgdata = img.get_image_data("YX",C=bgal_ch)
    except IndexError as e:
        raise IndexError(f"\nERROR: Image does not have B-gal input channel. Please check that B-gal channel parameters are correct.") from e

    ##### Create a boolean mask for pixels within the threshold #####

    mask = (imgdata >= bgal_thmin) & (imgdata <= bgal_thmax)
    
    ##### Calculate stats #####

    NpxPos = np.sum(mask).item()
    NpxTot = imgdata.size

    # Area calculations based on pixel info
    if pxInfo is True:
        AreaPos = NpxPos * pxarea * area_conv
        AreaTot = NpxTot * pxarea * area_conv
    else:
        AreaPos = None
        AreaTot = None
   
    
    RawIntDen = np.sum(imgdata[mask]).item()
    MeanIntens = np.mean(imgdata).item()
    
    return [NpxPos, NpxTot, AreaPos, AreaTot, pxarea, pxareaunit, RawIntDen,MeanIntens]


def subtract_background(image, radius : int) -> NDArray[np.generic]:
    """
    Subtracts the background signal of an input image using the rolling-ball
    algorithm with a previous downscaling to improve speed.

    This function applies mild uniform smoothing to the input image, estimates
    the background using a downscaled rolling-ball algorithm, and returns the background-subtracted image.

    This function follows the downscaling strategy proposed by scikit-image developers (see PR #7954 by jouyun):
    https://github.com/scikit-image/scikit-image/pull/7954

    The output image has the same dtype as the input image.

    Parameters
    ----------
        image : np.ndarray
            Input image array.
        radius : int
            Radius of the rolling ball, in pixels.

    Returns
    -------
        np.ndarray
            Background-subtracted image with the same dtype as the input.
    """
    # Pre-smoothing
    image = uniform_filter(image, size=3)

    # Convert image info
    img = np.asarray(image)
    img = img.astype(np.float64, copy = False)

    # Calculate downscale factor
    down_scale_factor = max(int(round(0.5 * (np.sqrt(radius) - 1))), 1)

    # Downscale
    img_small = rescale(
        img,
        1 / down_scale_factor,
        order=1,
        preserve_range=True,
        anti_aliasing=False,
    )

    # Apply rolling ball with scaled radius
    bg_small = rolling_ball(
        img_small,
        radius=radius / down_scale_factor,
        nansafe=False
    )

    # Upscale
    bg = resize(
        bg_small,
        img.shape,
        order=1,
        preserve_range=True,
        anti_aliasing=False,
    )

    bg = np.minimum(img, bg)

    # Subtract background
    result = img - bg

    return result.astype(image.dtype, copy=False)

def generate_biapy_input(img, nuclei_ch: int, apply_sbg: bool, sbg_rad: int, out_path: str):
    """
    Generates biapy input images.
    
    Parameters
    ----------
        img : BioImage
            Input image loaded as BioImage object.
        nuclei_ch : int
            Number of the channel containing nuclei staining signal.
        apply_sbg : bool
            If `True`, applies subtract background algorithm before saving BiaPy input.
        sbg_rad : int
            Radius of the rolling ball for the subtract background algorithm, in pixels.
        out_path : str
            Path indicating the folder where BiaPy input images will be stored.
    """
    # Load nuclei channel
    try:
        imgdata = img.get_image_data("YX",C=nuclei_ch)
    except IndexError as e:
        raise IndexError(f"\nERROR: Image does not have nuclei input channel. Please check that nuclei channel parameters are correct.") from e

    # Apply subtract background if needed
    if apply_sbg:
        OmeTiffWriter.save(
            subtract_background(imgdata,radius=sbg_rad),
            str(out_path),
            dim_order="YX")
    else:
        OmeTiffWriter.save(imgdata,
                           str(out_path),
                           dim_order="YX")
        

def load_input(input_folder: str) -> List:
    """
    This function accepts as parameter an existing path to a folder containing TIFF images and returns a list with the path to every image. It checks for folder existence, TIFF extension of files, and also skips hidden files that might be created by the OS in the folder.

    Parameters
    ----------
        input_folder: str
            Path to input folder containing TIFF images.

    Returns
    -------
        List:
            List of all (TIFF) files present in such folder.
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


def choose_threshold(df: pd.DataFrame) -> float:
    """
    This function helps the user to choose an specific nuclei threshold based on the density plot of nuclei area detected by BiaPy. First, it prints a kdeplot and asks the user to input a nuclei threshold (in um). Then it re-prints the plot showing the input threshold and asks the user for confirmation. If the user wants to change, by entering `n` it repeats the input procedure until positive (`y`) confirmation. Once given, it returns the nuclei threshold (in um) to be used in the `calculate_CTF` module.

    Parameters
    ----------
        df: pd.DataFrame
            DataFrame containing BiaPy nuclei quantification results.

    Returns
    -------
        float:
            Chosen nuclei threshold (in um).
    """

    df['area_um'] = df.area * df.PxArea
    sns.kdeplot(data=df, x="area_um", fill=True)
    plt.show(block=False)

    while True:
        user_input = input(f"\nEnter nuclei threshold (in um): ")
                
        try:
            thr = float(user_input)
            
            # If Jupyter, clear output
            if is_jupyter:
                clear_output(wait=True)
            
            # Close prior plots
            plt.close('all')
            
            # 3. Create the new plot
            sns.kdeplot(data=df, x="area_um", fill=True)
            plt.axvline(thr, color='red', linestyle='--', label=f'Threshold: {thr}')
            plt.title(f"Evaluating Threshold: {thr}")
            plt.legend()
            
            # Show plot
            plt.show(block=False)
            
            # 5. Get confirmation
            confirm = input("Keep this threshold? (y/n): ").lower()
            if confirm == 'y':
                plt.close('all')
                return thr
                
        except ValueError:
            print("Error: Please enter a valid number.")