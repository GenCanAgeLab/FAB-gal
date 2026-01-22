# FAB-Gal: Fluorescence Analysis of B-Galactosidase activity

🔥 **NEWS** 🔥: There are still not news

FAB-Gal is a method that joins the simplicity, fast execution, low cost and versatility of the classical SA-B-Gal assay with the quantification and sensitivity of fluorescence assays. It exploits the far-red fluorescence of X-Gal product ([Levitsky et al., 2013](https://www.nature.com/articles/srep02937)) coupled with a count of cell nuclei powered by [BiaPy](https://github.com/BiaPyX/BiaPy) to quantify in a semiautomatic and unbiased manner SA-B-Gal activity. This github repository is dedicated to the image analysis part od the method. To read the full method go to **url_paper**

⚠️ **GENERAL NOTES** ⚠️:\
 In order to obtain comparable results input images of FAB-Gal must be taken with the same objective and using the same laser intensity and time exposure for the far-red channel\
FAB-Gal only accepts **TIFF** images as input, if your images have a different format please convert them to TIFF using softwares such as ImageJ.

## FAB-Gal desktop application for Windows, MacOS

We have developed a desktop application that enables a user-friendly usage of our analysis method designed for a reduce number of images. If you want to use our methos in a large set of images (more than 50) you may consider to use the Jupyter notebook that we have designed fo this purpose (see below).

Windows app

MacOs app

Linux app

### 🧩 How does it work ?

Explanaition of the method Jandro

### 📝 User guide

#### Start and file selection

First of all download the FAB-Gal desktop application (see above).

Once you have opened the app first of all select the directory were you store the images that you want to analyze, you can filter the images you want by file name.

""Pictures of the process""

#### Background correction

Before you start to analyze your experimental image, measure background intensity. You can easily do it using the app, you only have to load an image with no cells or debries and set the B-Gal threshold to 0. Thus, the app will prompt the mean background fluorescence at the bottom of the image, copy the value and paste it on its textbox.

""Pictures of the process""

#### Adjusting thresholds

After setting the mean fluorescence background you should load your negative control image, that is an image of cells with nucelar staining but no stained with X-Gal. You will use this image to adjust the nuceli and B-Gal thresholds, using the sliders located at the top. 

The nuclei threshold is used to compute the number of nuclei and is preadjust using an algorithm so it should be almost fine. However, you can test different nuclei thresholds and modify the blur until the mask adjusts well to the nuceli of your image. Moreover, after nuclei segmentention you can filer by particle size to eliminate debris bigger or smaller than the nuclei.

⚠️ **NOTE** ⚠️: a subtract background transformation is done by default for enhancing nuceli segmentation. You should adjust the radius of the nuclei to math the one of your cell

The B-Gal threshold is used to correct cell autofluorescence, so after loading your negative control image you should set it to a value that makes the mean signal close to 0

⚠️ **NOTE** ⚠️: is not desirable to achieve 0 in your negative control because usually this implies setting a threshold that is to restrictive and will make you lose information. It is advisable to use a positive control to determine a reasonable threshold and avoid this problem

#### Measuring and saving results

Once you have set all the parameters, you can simply measure the images you want by clicking on "Measure" or measure all images (or the ones you filtered) by clicking on "Run all images". Once the images have been measured you can export your results to a csv file.

⚠️ **NOTE** ⚠️: in case you are measuring SA-B-gal of tissues the procedure is the same but without the nuclei channel, because similar cellularity is assumed when comparing the same tissue zone.

#### Interpreting the results

The results file contains the folowing columns:

**File**: name of the image file\
**Nuclei**: number of nuclei of the image\
**NpxPos**: number of positive pixels for the B-Gal threshold\
**RawIntDen**: total intensity of the positive pixels for the B-Gal threshold\
**Total Px**: total pixels of image\
**Area**: total area of the image\
**Area units**: units of the "Area" column (cm, µmm, nm, ...)\
**bgMF**: mean fluorescence of the background\
**CTF/nuclei**: corrected total fluorescence per nuceli (for culture cells), this unit is normalized by number of cells and corrected for autofluorescence and background so it can be directly used. Even if the app calculates it for every image we recommend to compute it at your biological replicate level (e.g. for each well), that is mergign the data of tecnichal replicates (this mitigates the effect of outlier images with a low number of cells )\
**CTF/area(px)**: corrected total fluorescence per pixel (for tissues sections), this unit is normalized by number of cells and corrected for autofluorescence and background so it can be directly used\
**CTF/area(unit)**: corrected total fluorescence per physical unit (for tissues sections), this unit is normalized by number of cells and corrected for autofluorescence and background so it can be directly used\

⚠️ **NOTE** ⚠️: the columns with physical size information (Area, Area units, CTF/area(unit)) are optional, to obtain them you have to click on **"button_name"**, thus the app will try to extract them form the image meatadata if it fails you can always introduce them manually. Otherwise, these columns will appear as NA.

## 🧑‍💻 FAB-Gal Jupyter notebook for high-throughput applications

For those users willing to apply FAB-Gal to a large number of images or for users that have computanional skills we provide a Jupyter notebook for a faster and more customable execution of FAB-Gal.

### 🧩 How does it work?

FAB-Gal pipeline reads the TIFF images and calculates the raw integrated density, the positive and total pixels and the total area of the far-red channel after applying a threshold defined by the user using a negative control (unstained cells). Then, if working with culture cells, it applies a subtract background transformation to the nuclei channel (set by default but optional) and uses it as input for BiaPy, a deep learning algorithm that segmentates nuclei. Finally, these data are used to calculate CTF per nuclei (culture cells) or CTF per area (tissue sections) of each image. However, as mentioned before when working with technical replicates (e.g. several phtos of the same well) we recommend to compute CTF per nuclei at the biological replicate level (e.g each well) to mitigate the effect of possible outlier values proceeding from images with low numbers of cells.

$$
\begin{aligned}
\text{CTF per nuclei}
&=
\frac{\sum \text{RawIntDen} - \sum \text{Positive pixels} \times \text{Mean background}}
{\sum \text{nuclei}}
\end{aligned}
$$

**RawIntDen** = Raw integrated density, the sum of the intensity values of the positive pixels  
**Mean background** = Mean intensity of an area with no cells or debris

$$
\begin{aligned}
\text{CTF per area}
&=
\frac{\sum \text{RawIntDen} - \sum \text{Positive pixels} \times \text{Mean background}}
{\sum \text{Area}}
\end{aligned}
$$

**Area** = the total area of the image in pixels or in physical size (e.g. µm)

⚠️ **NOTE** ⚠️: if you have a GPU we kindly recommend using it beacuse it accelerates nuclei segmentation

### 📝 User guide 

#### Installation

Roiz complete this please

## Authors

| Name | Role | Affiliations |
|------|------|--------------|
| Antonio G. Tartiere | Creator, Implementation | <ul><li>Departamento de Bioquímica y Biología Molecular, Instituto Universitario de Oncología (IUOPA), Universidad de Oviedo, Oviedo, Spain</li><li>Instituto de Investigación Sanitaria del Principado de Asturias (ISPA), Oviedo, Spain</li></ul> |
| David Roiz-Valle | Implementation | <ul><li>Departamento de Bioquímica y Biología Molecular, Instituto Universitario de Oncología (IUOPA), Universidad de Oviedo, Oviedo, Spain</li><li>Instituto de Investigación Sanitaria del Principado de Asturias (ISPA), Oviedo, Spain</li></ul> |
| José M. P. Freije | Supervision | <ul><li>Departamento de Bioquímica y Biología Molecular, Instituto Universitario de Oncología (IUOPA), Universidad de Oviedo, Oviedo, Spain</li><li>Instituto de Investigación Sanitaria del Principado de Asturias (ISPA), Oviedo, Spain</li><li>Centro de Investigación Biomédica en Red de Cáncer (CIBERONC), Madrid, Spain</li></ul> |
| Alejandro P. Ugalde | Creator, Implementation, Software design | <ul><li>Departamento de Bioquímica y Biología Molecular, Instituto Universitario de Oncología (IUOPA), Universidad de Oviedo, Oviedo, Spain</li><li>Instituto de Investigación Sanitaria del Principado de Asturias (ISPA), Oviedo, Spain</li></ul> |

## Citation

Please note that FAB-Gal is based on a publication. If you use it successfully for your research please be so kind to cite our work:




