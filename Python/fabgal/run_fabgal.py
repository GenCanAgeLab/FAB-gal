from pathlib import Path
from re import sub
from bioio import BioImage
from bioio_base.exceptions import UnsupportedFileFormatError
import logging
from dataclasses import asdict
import yaml

from .helpers import calculate_bgal,generate_biapy_input,load_input
from .run_biapy import run_biapy
from .calculate_CTF import calculate_CTF
from .config import FABGalConfig

#### Log options ####

import logging

# Set up logging 
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    datefmt="%H:%M:%S",
)

#  Silence Bioio warnings 
logging.getLogger("bioio").setLevel(logging.ERROR)

#  Get module logger 
logger = logging.getLogger(__name__)

#####################

def run_fabgal(cfg: FABGalConfig):
    """
    Main function to run Fab-Gal software

    Args:
        cfg (FABGalConfig): FABGal dataclass variable with all configuration options for running FAB-Gal.

    """
    ####### Load images #######
    
    myfiles = load_input(cfg.input_folder)

    ####### Iterate over files #######

    # Create results directory
    results_dir = Path(cfg.out_path) / f"Results_{cfg.experiment_name}"
    results_dir.mkdir(exist_ok = True, parents= True)

    # Create B-Gal results list
    bres = {}

    # Start B-Gal quantification message
    logger.info("Starting B-Gal quantification...")

    # Start iteration
    filenum = len(myfiles)
    for i, inf in enumerate(myfiles, start=1):

        print(f"Doing file {i} of {filenum}", end="\r", flush=True)

        ####### Open image #######
        try:
            img = BioImage(inf)
        except UnsupportedFileFormatError as e:
            raise UnsupportedFileFormatError("\nERROR: Image is not in TIFF format, convert it to TIFF using appropiate software (e.g. ImageJ)") from e
        except Exception as e:
            raise Exception(f"\nERROR: Could not open {inf.name}: {e}") from e
              
        ####### Generate input for BiaPy if needed #######

        if cfg.nuclei_ch is not None and cfg.run_biapy:

            # Create BiaPy input directory or check if it exists
            biapy_input = Path("biapy_input")
            biapy_input.mkdir(exist_ok=True)

            # Define path for generating BiaPy input image
            out_path = biapy_input / inf.name
            
            # Generate BiaPy input file
            generate_biapy_input(img, cfg.nuclei_ch, cfg.apply_subtract_background, cfg.sbg_rad, out_path)

        ####### Quantify B-Gal signal #######

        bres[sub(r'(?i)\.tiff?$', '', inf.name)] = calculate_bgal(img,
                                                                  cfg.bgal_ch,
                                                                  cfg.bgal_th,
                                                                  pxarea = cfg.pixel_area)


    ####### Output B-Gal calculations to file #######
    BGalres = results_dir / f"{cfg.experiment_name}_Raw_BGal_results.tsv"

    with BGalres.open('w') as bgal_f:

        # Create file header
        bgal_f.write("File\tNpxPos\tNpxTot\tAreaPos\tAreaTot\tPxArea\tBgal_RawIntDen\n")

        # Write calculations to file       
        for k, (npx, npxtot, areapos, areatot, pxarea, RawIntDen) in bres.items():
            bgal_f.write(f"{k}\t{npx}\t{npxtot}\t{areapos}\t{areatot}\t{pxarea}\t{RawIntDen}\n")

    # End of B-Gal quantification message
    logger.info("Finished B-Gal quantification")

    ####### Run BiaPy nuclei count if specified #######

    if cfg.nuclei_ch is not None and cfg.run_biapy:
        run_biapy(cfg)

    ####### Calculate CTF #######

    if cfg.nuclei_ch is not None:
        calculate_CTF(cfg)
    
    ####### Save config file #######

    cfg_dict = asdict(cfg)
    config_path = results_dir / "FABGal_config.yaml"

    with config_path.open("w") as f:
        yaml.safe_dump(
            asdict(cfg),
            f,
            sort_keys=False,
        )

    ####### End message #######
    logger.info("Finished analysis!")





