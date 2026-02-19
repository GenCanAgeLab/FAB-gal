from dataclasses import dataclass

###### Create pipeline config class ######

@dataclass
class FABgalConfig:
    """
    Dataclass for storing all parameters required for running FAB-gal
    
    ## Variables:

    ### Run parameters:
        **apply_subtract_background** (`bool`): If `True`, subtract background algorythm will be applied to input images before BiaPy runs.

        **run_name** (`str`): Name of run.

    ### Run parameters:
        **is_rerun** (`bool`): Set as `True` if the current run is a re-run of prior(s) run(s). If `True`, then the user needs to indicate which info (B-gal, BiaPy, or both) the program is reusing in `Bgal_run` and `Biapy_run`, respectively.

        **Bgal_run** (`str`): Name of the run from which to get B-gal results. If `None`, then B-gal quantification will be run again using current user configuration.

        **Biapy_run** (`str`): Name of the run from which to get BiaPy results. If `None`, then BiaPy nuclei count will be run again using current user configuration.

    ### File parameters:
        **input_folder** (`str`): Path to folder containing input images in `TIFF` format.

        **experiment_name** (`str`): Name of the experiment.

        **img_to_ind** (`str`): Path to `.tsv` file with two columns: `File` and `Individual`, each containing filenames (eg. `image1.tif`) and the corresponding individual (eg. `WT 1`). If `img_to_ind` is not `None`, then FAB-gal will compute CTF per individual as well as per image, and return them in a separate file.

        **out_path** (`str`): Output folder where the results folder will be generated.

    ### Image parameters:
        **nuclei_ch** (`int`): Number of the channel containing nuclei staining signal. If `None` (eg. tissue images with no nuclei counterstain), then FAB-gal will not run the nuclei count and will not compute CTF per nuclei.

        **bgal_ch** (`int`): Number of the channel containing B-gal staining signal.

        **pixel_area** (`float`): Pixel size of input images. For metadata-enriched images, it can be left as `None`, and FAB-gal will get it from metadata information. If there is no `pixel_area`, either from user input or metadata, then it will default to 1, and a warning will be printed.

    ### β-Gal parameters:
        **backgr_val** (`int`): β-Gal mean fluorescence intensity of a background image. If `None`, FAB-gal will try to get it from `backgr_img`. If `backgr_img` is also `None`, then it will print a warning and not compute CTF.

        **backgr_img** (`str`): Name of the image to be used to calculate the β-Gal background mean fluorescence intensity.

        **bgal_th** (`int`): Intensity value to use as threshold to determine positive β-Gal pixels.

    ### Subtract background parameters (for nuclei channel processing)
        **sbg_rad** (`int`): Radius of subtract background algorithm to be used if subtract_background = `True`. Should be adjusted depending on nuclei size.

    ### BiaPy parameters
        **config_file** (`str`): Path to BiaPy YAML configuration file.

        **gpu** (`str`): GPU to use. Obtained from `nvidia-smi` command.

        **keep_masks** (`bool`): If `True`, masks from BiaPy nuclei quantification will be kept in the results folder.

        **nuclei_thr** (`float`): Area value for a predicted nuclei to be considered as such. Has to be in the same units as the `pixel_area` the user has input or the image metadata contains. If none, it executes a function to explore nuclei area distribution and inputing a threshold value based on it.
    """

    # Run parameters
    apply_subtract_background: bool
    run_name: str

    # Re-run options
    is_rerun: bool
    Bgal_run: str | None
    Biapy_run: str | None

    #File parameters
    input_folder: str
    experiment_name: str
    img_to_ind: str
    out_path: str

    # Image parameters
    nuclei_ch: int | None
    bgal_ch: int 
    pixel_area: float | None

    # β-Gal parameters
    backgr_val: float | None 
    backgr_img: str | None
    bgal_th: int

    # Subtract background parameters (for nuclei channel processing)
    sbg_rad: int

    # BiaPy parameters
    config_file: str
    gpu: str
    keep_masks: bool
    nuclei_thr: float | None