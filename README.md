# FAB-Gal: Fluorescence Analysis of B-Galactosidase activity

🔥 **NEWS** 🔥: There are still not news

FAB-Gal is a method that joins the simplicity, fast execution, low cost and versatility of the classical SA-B-Gal assay with the quantification and sensitivity of fluorescence assays. It exploits the far-red fluorescence of X-Gal product (Levitsky et al., 2013) coupled with a count of cell nuclei powered by BiaPy (https://github.com/BiaPyX/BiaPy) to quantify in a semiautomatic and unbiased manner SA-B-Gal activity. This github repository is dedicated to the image analysis part od the method. To read the full method go to **url_paper**

## FAB-Gal desktop application for Windows, MacOS

We have developed a desktop application that enables a user-friendly usage of our analysis method designed for a reduce number of images. If you want to use our methos in a large set of images (more than 50) you may consider to use the Jupyter notebook that we have designed fo this purpose (see below).

Windows app

MacOs app

Linux app

### How does it work ?

Explanaition of the method

### User guide

First of all download the FAB-Gal desktop application (see above).

Once you have opened the app first of all select the directory were you store the images that you want to analyze, you can filter the images you want by file name.

""Pictures of the process""

⚠️ **NOTE** ⚠️: the app only accepts **TIFF** imageges, if your images have a different format please convert them to TIFF using softwares such as ImageJ

Before you start to analyze your experimental image, measure background intensity. You can easily do it using the app, you only have to load an image with no cells or debries and set the B-Gal threshold to 0. Thus, the app will prompt the mean background fluorescence at the bottom of the image, copy the value and paste it on its textbox.

""Pictures of the process""

After setting the mean fluorescence background you should load your negative control image, that is an image of cells with nucelar staining but no stained with X-Gal. You will use this image to adjust the nuceli and B-gal thresholds, using the sliders located at the top. 

The nuclei threshold is used to compute the number of nuclei and is preadjust using an algorithm so it should be almost fine. However, you can test different nuclei thresholds and modify the blur until the mask adjusts well to the nuceli of your image.

⚠️ **NOTE** ⚠️: a subtract background transformation is done by default for enhancing nuceli segmentation. You should adjust the radius of the nuclei to math the one of your cell

## 🧑‍💻 FAB-Gal Jupyter notebook for high-throughput applications 




