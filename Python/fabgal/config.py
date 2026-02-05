from dataclasses import dataclass

###### Create pipeline config class ######

@dataclass
class FABGalConfig:

    # Run parameters
    delete_intermediate_files: bool
    apply_subtract_background: bool
    run_biapy: bool

    #File parameters
    input_folder: str
    experiment_name: str
    img_to_ind: str
    out_path: str

    # Image parameters
    nuclei_ch: int | None
    bgal_ch: int 
    pixel_area: float | None

    # β-Gal background parameters
    backgr_val: int | None 
    backgr_img: str | None

    # β-Gal threshold for positive signal
    bgal_th: int

    # Substract background parameters (for nuclei channel processing)
    sbg_rad: int

    # Biapy parameters
    config_file: str
    run_id: int      
    gpu: str
    keep_images: bool
    nuclei_thr: float