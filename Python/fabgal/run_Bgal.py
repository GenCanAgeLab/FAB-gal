from pathlib import Path
from re import sub
from bioio import BioImage
from bioio_base.exceptions import UnsupportedFileFormatError
from bioio.writers import OmeTiffWriter


from .helpers import calculate_bgal,generate_biapy_input,load_input
from .config import FABgalConfig


#### Log options ####

import logging

logger = logging.getLogger(__name__)


def run_Bgal(cfg: FABgalConfig):
    """
    Run B-gal quantification analysis with stated parameters, returning a results file for all images.

    Args:
        cfg (FABgalConfig): FABgal dataclass variable with all configuration options for running FAB-gal.
    """
    ####### Load images #######

    results_dir = Path(cfg.out_path) / f"Results_{cfg.experiment_name}" / cfg.run_name
    myfiles = load_input(cfg.input_folder)
    
    ####### Iterate over files performing quantification #######

    # Create B-Gal results list
    bres = {}

    ################# Start iteration #####################
    
    filenum = len(myfiles)
    for i, inf in enumerate(myfiles, start=1):

        print(f"Doing file {i} of {filenum}", end="\r", flush=True)

        ####### Open image #######
        try:
            img = BioImage(inf)
        except UnsupportedFileFormatError as e:
            raise UnsupportedFileFormatError("\nERROR: Image is not in TIFF format, convert it to TIFF using appropriate software (e.g. ImageJ)") from e
        except Exception as e:
            raise Exception(f"\nERROR: Could not open {inf.name}: {e}") from e
            
        ####### Generate input for BiaPy if needed #######
        if cfg.is_rerun and cfg.Biapy_run is not None :
            pass # Skip BiaPy input generation (prior BiaPy run)
        elif cfg.nuclei_ch is None:
            pass # Skip BiaPy input generation (no nuclei ch info)
        else:
            # Create BiaPy input directory or check if it exists
            biapy_input = Path("biapy_input")
            biapy_input.mkdir(exist_ok=True)

            # Define path for generating BiaPy input image
            out_path = biapy_input / inf.name
            
            # Generate BiaPy input file
            generate_biapy_input(img, cfg.nuclei_ch, cfg.apply_subtract_background, cfg.sbg_rad, out_path)

        ####### Generate OME_metadata from image to load #######

        # Get image parameters to generate metadata
        shapes = [img.data.shape]
        dtypes = [img.data.dtype]
        dim_order = img.dims.order
        channel_names = img.channel_names if hasattr(img, 'channel_names') else None

        # Build_OME
        ome_metadata_variable = OmeTiffWriter.build_ome(
            shapes,
            dtypes,
            channel_names=[channel_names] if channel_names else None, # Needs to be a list of lists if provided
            image_name=["img"],     # List of names
            dimension_order=[dim_order] # List of strings
        )

        pxunit = ome_metadata_variable.images[0].pixels.physical_size_x_unit.value

        ####### Quantify B-Gal signal #######

        bres[sub(r'(?i)\.tiff?$', '', inf.name)] = calculate_bgal(img,
                                                                cfg.bgal_ch,
                                                                cfg.bgal_th,
                                                                pxarea = cfg.pixel_area,
                                                                pxunit = pxunit)
    
    ################# End iteration #####################


    ####### Output B-Gal calculations to file #######
    Bgalres = results_dir / "Raw_Bgal_results.tsv"

    with Bgalres.open('w', encoding="utf-8") as bgal_f:

        # Create file header
        bgal_f.write("File\tNpxPos\tNpxTot\tAreaPos\tAreaTot\tPxArea\tPxAreaUnits\tBgal_RawIntDen\tMean_Intens\n")

        # Write calculations to file       
        for k, (npx, npxtot, areapos, areatot, pxarea, pxareaunit, RawIntDen,MeanIntens) in bres.items():
            bgal_f.write(f"{k}\t{npx}\t{npxtot}\t{areapos}\t{areatot}\t{pxarea}\t{pxareaunit}\t{RawIntDen}\t{MeanIntens}\n")

