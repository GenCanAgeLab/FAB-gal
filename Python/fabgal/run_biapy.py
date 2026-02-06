from contextlib import redirect_stdout, redirect_stderr
from pathlib import Path
import logging
from biapy import BiaPy
import pandas as pd
from skimage import io
from skimage.color import label2rgb
import numpy as np
import shutil
from .config import FABGalConfig

#### Log options ####

import logging

logger = logging.getLogger(__name__)


def run_biapy(cfg: FABGalConfig):
    """
    Run BiaPy analysis to count nuclei with stated parameters, processes output, and generates final .

    Args:
        cfg (FABGalConfig): FABGal dataclass variable with all configuration options for running FAB-gal.
    """
    # Results dir for BiaPy
    results_dir = Path(cfg.out_path) / f"Results_{cfg.experiment_name}"
    biapy_result_dir = results_dir / f"biapy_output_run{cfg.run_id}"
    biapy_log = results_dir / f"biapy_{cfg.experiment_name}_run{cfg.run_id}.log"

    # If there is a previous execution of BiaPy with the same name, delete previous subtracted images
    sub_images_out = results_dir / f"subtracted_images_run{cfg.run_id}"
    if sub_images_out.exists():
        shutil.rmtree(sub_images_out)

    ########### Create and run the BiaPy job ###########

    # Start message
    logger.info("Starting BiaPy nuclei quantification...")

    # Catch stdout and stderr from BiaPy and redirect it
    with biapy_log.open("w") as f, redirect_stdout(f), redirect_stderr(f):

        biapy = BiaPy(cfg.config_file,
                    result_dir = biapy_result_dir,
                    name = cfg.experiment_name,
                    run_id = cfg.run_id,
                    gpu = cfg.gpu)
    
    
        biapy.run_job()
    
    # Final message
    logger.info("Finished BiaPy quantification")


    ########### Process BiaPy output ###########
    
    ## Define BiaPy output images folder
    biapyout = biapy_result_dir / cfg.experiment_name / "results" / f"{cfg.experiment_name}_{cfg.run_id}" / "per_image_instances"

    ## Define Nuclei stats output file
    biapystatsfile = results_dir / f"{cfg.experiment_name}_BiaPy_results.tsv"

    #### Process nuclei stats ####

    # List of csv with filename column
    biapy_csv = list(biapyout.glob("*full_stats.csv"))

    dfs = []
    for inf in biapy_csv:
        df = pd.read_csv(inf)
        df["File"] = inf.stem.replace("_full_stats","")
        dfs.append(df)

    # Concatenate dfs
    result = pd.concat(dfs, ignore_index=True)
    result = result.drop(columns=["conditions"])
    result.to_csv(str(biapystatsfile), sep="\t")

    # Delete original csv
    for inf in biapy_csv:
        inf.unlink()

    ########### Arrange intermediate files and images ###########
    
    #### Process output mask images ####

    if cfg.keep_images:
        for inf in biapyout.glob("*.tif"):
            img = io.imread(inf)
            img = label2rgb(img)*255
            img = img.astype(np.uint8)
            io.imsave(biapy_result_dir / f"{inf.stem}_mask.png",img, check_contrast = False)
            inf.unlink()
    
    #### Save subtracted images from BiaPy input (if generated) ####

    if cfg.apply_subtract_background:
        Path("biapy_input").rename(sub_images_out)

    ## If BiaPy ended correctly and sustract_background was not active, we delete biapy_input folder
    if not cfg.apply_subtract_background:
        shutil.rmtree("biapy_input")
        
    #### Move BiaPy config yml to results ####

    biapy_config = Path(biapy_result_dir / cfg.experiment_name / "config_files")

    for file in biapy_config.iterdir():
        file.rename(results_dir / "BiaPy_config.yaml")

    ### Delete extra files not used in the analysis ###

    if cfg.delete_intermediate_files:
        shutil.rmtree(biapy_result_dir / cfg.experiment_name)