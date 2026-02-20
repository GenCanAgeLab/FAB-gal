from contextlib import redirect_stdout, redirect_stderr
from pathlib import Path
import logging
from biapy import BiaPy
import pandas as pd
from skimage import io
from skimage.color import label2rgb
import numpy as np
import shutil
from .config import FABgalConfig

#### Log options ####

import logging

logger = logging.getLogger(__name__)


def run_biapy(cfg: FABgalConfig):
    """
    Run BiaPy analysis to count nuclei with stated parameters, processes output, and generates summary count info for all samples.

    Args:
        cfg (FABgalConfig): FABgal dataclass variable with all configuration options for running FAB-gal.
    """
    # Results dir for BiaPy
    results_dir = Path(cfg.out_path) / f"Results_{cfg.experiment_name}" / cfg.run_name
    
    biapy_log = results_dir / f"biapy_runInfo.log"

    ########### Create and run the BiaPy job ###########

    # Catch stdout and stderr from BiaPy and redirect it
    with biapy_log.open("w", encoding="utf-8") as f, redirect_stdout(f), redirect_stderr(f):

        biapy = BiaPy(cfg.config_file,
                    result_dir = results_dir,
                    name = "1",
                    run_id = 1,
                    gpu = cfg.gpu)
    
    
        biapy.run_job()
    
    ########### Process BiaPy output ###########
    
    ## Create BiaPy final output folder

    Path(results_dir / "BiaPy_output").mkdir(exist_ok=True)

    ## Define BiaPy output images folder
    biapyout = results_dir / "1" / "results" / "1_1" / "per_image_instances"

    ## Define Nuclei stats output file
    biapystatsfile = results_dir / "BiaPy_results.tsv"

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
    result.to_csv(str(biapystatsfile), sep="\t", encoding="utf-8")

    # Delete original csv
    for inf in biapy_csv:
        inf.unlink()

    ########### Arrange intermediate files and images ###########
    
    #### Process output mask images ####

    if cfg.keep_masks:
        for inf in biapyout.glob("*.tif"):
            img = io.imread(inf)
            img = label2rgb(img)*255
            img = img.astype(np.uint8)
            io.imsave(results_dir / "BiaPy_output" / f"{inf.stem}_mask.png",img, check_contrast = False)
            inf.unlink()
    
    #### Save subtracted images from BiaPy input (if generated) ####

    sub_images_out = results_dir / f"subtracted_images"

    if sub_images_out.exists():
        shutil.rmtree(sub_images_out)

    if cfg.apply_subtract_background:
        Path("biapy_input").replace(sub_images_out)

    ## If BiaPy ended correctly and sustract_background was not active, we delete biapy_input folder
    if not cfg.apply_subtract_background:
        shutil.rmtree("biapy_input")
        
    #### Move BiaPy config yml to results ####

    biapy_config = Path(results_dir / "1" / "config_files")

    for file in biapy_config.iterdir():
        file.replace(results_dir / "BiaPy_output" / "BiaPy_config.yaml")

    ### Delete extra files not used in the analysis ###

    shutil.rmtree(results_dir / "1")
