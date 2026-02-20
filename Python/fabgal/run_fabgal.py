from pathlib import Path
import logging
import warnings
from bioio import BioImage
from bioio_base.exceptions import UnsupportedFileFormatError
from dataclasses import asdict
import yaml

from .run_Bgal import run_Bgal
from .run_biapy import run_biapy
from .calculate_CTF import calculate_CTF
from .config import FABgalConfig
from .helpers import load_input,generate_biapy_input

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
warnings.filterwarnings("ignore", module="ome_types")
warnings.filterwarnings("ignore", module="pydantic")
warnings.filterwarnings("ignore", module="tifffile")

#  Get module logger 
logger = logging.getLogger(__name__)

#####################

def run_fabgal(cfg: FABgalConfig):
    """
    Main function to run Fab-Gal software

    Args:
        cfg (FABgalConfig): FABgal dataclass variable with all configuration options for running FAB-gal.

    """
    ####### Check config #######

    if cfg.Biapy_run and cfg.nuclei_ch:
        raise Exception("ERROR: Entering BiaPy prior info in BiaPy_run, but nuclei_ch is None. Please check FABgal config.")
    
    if cfg.is_rerun:
        if cfg.Bgal_run is None:
            if cfg.Biapy_run is None:
               raise Exception("ERROR: Current run is a rerun, but the run from which to get B-gal and/or BiaPy info is not indicated. Please check FABgal config.") 
    
    if not cfg.is_rerun:
        if cfg.Bgal_run is not None:
            raise Exception("ERROR: Entering B-gal prior info in Bgal_run, but is_rerun is False. Please check FABgal config.")
        if cfg.Biapy_run is not None:
            raise Exception("ERROR: Entering BiaPy prior info in Bgal_run, but is_rerun is False. Please check FABgal config.")


    ####### Create results directory #######

    results_dir = Path(cfg.out_path) / f"Results_{cfg.experiment_name}" / cfg.run_name
    results_dir.mkdir(exist_ok = True, parents= True)

    ####### Load images #######
    
    myfiles = load_input(cfg.input_folder)

    ####### Iterate over files #######

    ####### Run B-gal quantification ######

    if cfg.is_rerun and cfg.Bgal_run is not None :

        # Skip B-gal quantification
        logger.info(f"Using B-gal quantification values from run: {cfg.Bgal_run}")

        ####### If the user wants to run BiaPy again, we have to generate BiaPy input #######

        if not cfg.Biapy_run:
            logger.info(f"Generating BiaPy input images")

            # Create BiaPy input directory
            biapy_input = Path("biapy_input")
            biapy_input.mkdir(exist_ok=True)

            # Load input files
            myfiles = load_input(cfg.input_folder)

            # Iterate over files
            for inf in myfiles:

                # Open image
                try:
                    img = BioImage(inf)
                except UnsupportedFileFormatError as e:
                    raise UnsupportedFileFormatError("\nERROR: Image is not in TIFF format, convert it to TIFF using appropriate software (e.g. ImageJ)") from e
                except Exception as e:
                    raise Exception(f"\nERROR: Could not open {inf.name}: {e}") from e

                # Define path for generating BiaPy input image
                biapy_out_path = biapy_input / inf.name
                
                # Generate BiaPy input file
                generate_biapy_input(img, cfg.nuclei_ch, cfg.apply_subtract_background, cfg.sbg_rad, biapy_out_path)
            
            logger.info(f"Done!")

        #######################################################################################

    else:
        # Start B-Gal quantification message
        logger.info("Starting B-Gal quantification...")

        run_Bgal(cfg)

        # End of B-Gal quantification message
        logger.info("Finished B-Gal quantification")


    ####### Run BiaPy nuclei count if specified #######

    if cfg.is_rerun and cfg.Biapy_run is not None:
        # Skip B-gal quantification
        logger.info(f"Using BiaPy nuclei quantification from run: {cfg.Biapy_run}")
    elif cfg.nuclei_ch is None:
        # Skip BiaPy nuclei quantification (no nuclei info)
        logger.info("No nuclei channel, skipping nuclei quantification")
    else:

        # Start message
        logger.info("Starting BiaPy nuclei quantification...")

        run_biapy(cfg)

        # Final message
        logger.info("Finished BiaPy quantification")
    
    ####### Calculate CTF #######

    if cfg.nuclei_ch is not None:

        # Start message
        logger.info("Starting data processing and CTF calculation...")

        calculate_CTF(cfg)

        # End message
        logger.info("Finished data processing")
    
    ####### Save config file #######

    cfg_dict = asdict(cfg)
    config_path = results_dir / "FABgal_config.yaml"

    with config_path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(
            asdict(cfg),
            f,
            sort_keys=False,
        )

    ####### End message #######
    logger.info("Finished analysis!")





