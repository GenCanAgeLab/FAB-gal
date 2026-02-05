from pathlib import Path
import pandas as pd
import logging
from .config import FABGalConfig



#### Log options ####

import logging

logger = logging.getLogger(__name__)



def calculate_CTF(cfg: FABGalConfig):
    """
    Process B-Gal calculations and BiaPy nuclei count to get the corrected total fluorescence (CTF) from each individual image. If an `image to individual`file is entered, also returns CTF calculations per individual, after adding the measurements of all image replicates from each individual.

    Args:
        cfg (FABGalConfig): FABGal dataclass variable with all configuration options for running FAB-Gal.
    """

    # Start message
    logger.info("Starting data processing and CTF calculation...")

    # Load BGal and nuclei stats
    results_dir = Path(cfg.out_path) / f"Results_{cfg.experiment_name}"
    biapystatsfile = results_dir / f"{cfg.experiment_name}_BiaPy_results.tsv"

    if not biapystatsfile.exists():
        raise FileNotFoundError("\nERROR: BiaPy output file cannot be found. Please check that BiaPy has been executed.")

    nucleidf = pd.read_table(biapystatsfile)
    bgaldf = pd.read_table(results_dir / f"{cfg.experiment_name}_Raw_BGal_results.tsv")

    # Filter nuclei below area threshold
    nucleidf_pxarea = pd.merge(nucleidf,bgaldf[['File','PxArea']],how='inner',on='File')
    nucleidf_pxarea['nucl_thr_pixel'] = cfg.nuclei_thr / nucleidf_pxarea.PxArea

    nucleidf_filt = nucleidf_pxarea[nucleidf_pxarea.area > nucleidf_pxarea.nucl_thr_pixel]

    # Count nuclei per image file 
    nucleitot = nucleidf_filt['File'].value_counts()

    # B-Gal background intensity
    computeCTF = True
    if cfg.backgr_val is not None:
        BGal_backgr = cfg.backgr_val
    elif cfg.backgr_img is not None:
        BGal_backgr = bgaldf.at[cfg.backgr_img, 'BGal_RawIntDen'] / bgaldf.at[cfg.backgr_img, 'NpxTot']
    else:
        computeCTF = False

    # Merge nuclei and B-Gal data
    resdf = pd.merge(bgaldf,nucleitot,how='outer',on='File')
    resdf = resdf.rename(columns = {'count':'NumNucl'})

    # Calculate CTF on individual images (only if background intensity is present)
    if computeCTF:
        CTFimg = resdf
        CTFimg['bgMF'] = BGal_backgr
        CTFimg['CTFnucl'] = (CTFimg.Bgal_RawIntDen - CTFimg.NpxPos * CTFimg.bgMF) / CTFimg.NumNucl
        CTFimg['CTFpix'] = (CTFimg.Bgal_RawIntDen - CTFimg.NpxPos * CTFimg.bgMF) / CTFimg.NpxTot
        CTFimg['CTFarea'] = (CTFimg.Bgal_RawIntDen - CTFimg.NpxPos * CTFimg.bgMF) / CTFimg.AreaTot

        CTFimg.to_csv(results_dir / f"{cfg.experiment_name}_results_perimage.tsv", sep="\t")

        # Load individual info (if present)
        if cfg.img_to_ind is not None:
            img_ind_df = pd.read_table(cfg.img_to_ind)
            img_ind_df["File"] = img_ind_df["File"].str.replace(r"(?i)\.tiff?$","",regex=True)
            resdf = pd.merge(resdf,img_ind_df,how='outer',on='File')

            CTFind = resdf.groupby(['Individual']).sum(numeric_only = True)
     
            CTFind['bgMF'] = BGal_backgr
            CTFind['CTFnucl'] = (CTFind.Bgal_RawIntDen - CTFind.NpxPos * CTFind.bgMF) / CTFind.NumNucl
            CTFind['CTFpix'] = (CTFind.Bgal_RawIntDen - CTFind.NpxPos * CTFind.bgMF) / CTFind.NpxTot
            CTFind['CTFarea'] = (CTFind.Bgal_RawIntDen - CTFind.NpxPos * CTFind.bgMF) / CTFind.AreaTot

            CTFind.to_csv(results_dir / f"{cfg.experiment_name}_results_perindividual.tsv")
    
    # End message
    logger.info("Finished data processing")