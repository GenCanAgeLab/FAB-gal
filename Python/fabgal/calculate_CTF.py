from pathlib import Path
import pandas as pd
from re import sub
import logging
from .config import FABgalConfig
from .helpers import choose_threshold



#### Log options ####

import logging

logger = logging.getLogger(__name__)



def calculate_CTF(cfg: FABgalConfig):
    """
    Process B-Gal calculations and BiaPy nuclei count to get the corrected total fluorescence (CTF) from each individual image. If an `image to individual`file is entered, also returns CTF calculations per individual, after adding the measurements of all image replicates from each individual.

    Args:
        cfg (FABgalConfig): FABgal dataclass variable with all configuration options for running FAB-gal.
    """

    ########## Load Bgal and nuclei stats ##########

    ### Define result dir ###

    results_dir = Path(cfg.out_path) / f"Results_{cfg.experiment_name}" / cfg.run_name

    ### Define B-gal and BiaPy results file ###

    if cfg.is_rerun:
        # B-gal results (prior or newly generated)
        if cfg.Bgal_run is not None:
            bgalquantfile = Path(cfg.out_path) / f"Results_{cfg.experiment_name}" / cfg.Bgal_run / "Raw_Bgal_results.tsv"
        else:
            bgalquantfile = results_dir / "Raw_Bgal_results.tsv"
        
        # BiaPy results (prior or newly generated)
        if cfg.Biapy_run is not None:
            biapystatsfile = Path(cfg.out_path) / f"Results_{cfg.experiment_name}" / cfg.Biapy_run / "BiaPy_results.tsv"
        else:
            biapystatsfile = results_dir / "BiaPy_results.tsv"
    else:
        # Newly generated results for both if is_rerun is false
        bgalquantfile = results_dir / "Raw_Bgal_results.tsv"
        biapystatsfile = results_dir / "BiaPy_results.tsv"

    ### Check if B-gal and BiaPy results file exist ###

    if not bgalquantfile.exists():
        if cfg.is_rerun and cfg.Bgal_run is not None:
            raise FileNotFoundError(f"\nERROR: Run {cfg.Bgal_run} cannot be found or B-gal output file is not located in it. Please review run {cfg.Bgal_run} results and check FABgal config.")
        else:
            raise FileNotFoundError("\nERROR: Cannot find newly generated B-gal output file. Please check that B-gal quantification ended correctly.")

    if not biapystatsfile.exists():
        if cfg.is_rerun and cfg.Biapy_run is not None:
            raise FileNotFoundError(f"\nERROR: Run {cfg.Biapy_run} cannot be found or BiaPy output file is not located in it. Please review run {cfg.Biapy_run} results and check FABgal config.")
        else:
            raise FileNotFoundError("\nERROR: Cannot find newly generated BiaPy output file. Please check that BiaPy ended correctly.")

    ### Load B-gal and BiaPy results ###

    nucleidf = pd.read_table(biapystatsfile)
    bgaldf = pd.read_table(bgalquantfile)

    # Filter nuclei below area threshold
    nucleidf_pxarea = pd.merge(nucleidf,bgaldf[['File','PxArea']],how='inner',on='File')

    ## If cfg.nuclei_thr is None, then choose interactively the nuclei thr
    if cfg.nuclei_thr is None:
        cfg.nuclei_thr = choose_threshold(nucleidf_pxarea)
        nucleidf_pxarea['nucl_thr_pixel'] = cfg.nuclei_thr / nucleidf_pxarea.PxArea
        nucleidf_filt = nucleidf_pxarea[nucleidf_pxarea.area > nucleidf_pxarea.nucl_thr_pixel]
    else:
        nucleidf_pxarea['nucl_thr_pixel'] = cfg.nuclei_thr / nucleidf_pxarea.PxArea
        nucleidf_filt = nucleidf_pxarea[nucleidf_pxarea.area > nucleidf_pxarea.nucl_thr_pixel]

    # Count nuclei per image file 
    nucleitot = nucleidf_filt['File'].value_counts()

    # B-Gal background intensity
    computeCTF = True
    if cfg.backgr_val is not None:
        Bgal_backgr = cfg.backgr_val
    elif cfg.backgr_img is not None:
        backgr_img = sub(r'(?i)\.tiff?$', '', cfg.backgr_img)
        Bgal_backgr = bgaldf.loc[bgaldf['File'] == backgr_img, 'Mean_Intens'].values[0]       
    else:
        computeCTF = False
        print("No B-gal background information supplied. Will not compute CTF.")

    # Merge nuclei and B-Gal data
    resdf = pd.merge(bgaldf,nucleitot,how='outer',on='File')
    resdf = resdf.rename(columns = {'count':'NumNucl'})

    # Calculate CTF on individual images (only if background intensity is present)
    if computeCTF:
        CTFimg = resdf
        CTFimg['bgMF'] = Bgal_backgr
        CTFimg['CTFnucl'] = (CTFimg.Bgal_RawIntDen - CTFimg.NpxPos * CTFimg.bgMF) / CTFimg.NumNucl
        CTFimg['CTFpix'] = (CTFimg.Bgal_RawIntDen - CTFimg.NpxPos * CTFimg.bgMF) / CTFimg.NpxTot
        CTFimg['CTFarea'] = (CTFimg.Bgal_RawIntDen - CTFimg.NpxPos * CTFimg.bgMF) / CTFimg.AreaTot

        CTFimg.to_csv(results_dir / f"CTF_perimage.tsv", sep="\t")

        # Load individual info (if present)
        if cfg.img_to_ind is not None:
            img_ind_df = pd.read_table(cfg.img_to_ind)
            img_ind_df["File"] = img_ind_df["File"].str.replace(r"(?i)\.tiff?$","",regex=True)
            resdf = pd.merge(resdf,img_ind_df,how='outer',on='File')

            CTFind = resdf.groupby(['Individual']).agg(
                NpxPos = ("NpxPos","sum"),
                NpxTot = ("NpxTot","sum"),
                AreaPos = ("AreaPos","sum"),
                AreaTot = ("AreaTot","sum"),
                PxArea = ("PxArea","first"), # All images from the same ind. must have the same PixArea
                Bgal_RawIntDen = ("Bgal_RawIntDen","sum"),
                Mean_Intens = ("Mean_Intens","mean"),
                NumNucl = ("NumNucl","sum"),
                NumImages=("File", "count")
            )
     
            CTFind['bgMF'] = Bgal_backgr
            CTFind['CTFnucl'] = (CTFind.Bgal_RawIntDen - CTFind.NpxPos * CTFind.bgMF) / CTFind.NumNucl
            CTFind['CTFpix'] = (CTFind.Bgal_RawIntDen - CTFind.NpxPos * CTFind.bgMF) / CTFind.NpxTot
            CTFind['CTFarea'] = (CTFind.Bgal_RawIntDen - CTFind.NpxPos * CTFind.bgMF) / CTFind.AreaTot

            CTFind.to_csv(results_dir / "CTF_perindividual.tsv", sep="\t")