// ================= USER PARAMETERS =================

// Input / Output
inputDir = getDirectory("Choose a Directory");
inputDir = replace(inputDir, '\\','/');
files = getFileList(inputDir);
outputDir = getDirectory("Choose a Directory");
outputDir = replace(outputDir, '\\','/');
roiDir    = outputDir + "ROIs/";
resDir    = outputDir + "Results/";

// Channels (1-based indexing)
cellChannel   = 3;   // channel used to create cell mask
nucleiChannel = 1;   // nuclei channel
measureChannels = newArray(2,3,4); // channels to measure

// Thresholds
cellThreshMin   = 60;
cellThreshMax   = 255;

nucleiThreshMin = 15;
nucleiThreshMax = 255;

// ================================================

// Create output folders
File.makeDirectory(outputDir);
File.makeDirectory(roiDir);
File.makeDirectory(resDir);

// Clear results
run("Clear Results");
//NumberofRows = Table.size("Summary");
//Table.deleteRows(0, NumberofRows-1, “Summary”);
run("Set Measurements...",
        "area mean integrated area_fraction display redirect=None decimal=3");

for (i=0; i<files.length; i++){
//for (i=0; i<4; i++){
	if (endsWith(files[i], '.czi')) {
		//Import the image
		run("Bio-Formats", "open=[" + inputDir + files[i] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack");
		rename(files[i]);
		run("Split Channels");
		// Cell mask
		selectWindow("C" + cellChannel + "-" + files[i]);
		run("Duplicate...", "title=cell_channel");
		setThreshold(cellThreshMin, cellThreshMax);  
		setOption("BlackBackground", false);
		run("Convert to Mask"); 
		run("Fill Holes"); 
		run("Close-"); 
		rename("cell_mask");
		// Nuclei mask
		selectWindow("C" + nucleiChannel + "-" + files[i]);
		run("Duplicate...", "title=cell_channel");
		run("Subtract Background...", "rolling=50");
		setThreshold(nucleiThreshMin, nucleiThreshMax);  
		setOption("BlackBackground", false);
		run("Convert to Mask"); 
		run("Fill Holes"); 
		run("Watershed");
		run("Analyze Particles...", "size=0-80 add");
		roiManager("Fill");
		roiManager("Reset");
		rename("nuclei_mask");
		// Subtract nuclei mask to cell mask
		imageCalculator("Subtract create", "cell_mask", "nuclei_mask");
    	rename("final_mask");
    	// Create and store ROI
    	setThreshold(122, 255);  
		setOption("BlackBackground", false);
		run("Convert to Mask");
		// Alternatively you can run make binary
		run("Analyze Particles...", "size=500-Infinity show=Nothing add");
		roiManager("Save", roiDir + files[i] + "_ROIs.zip");
		// Measure channels
		for (c = 0; c < measureChannels.length; c++) {
		selectWindow("C" + measureChannels[c] + "-" + files[i]);
        roiManager("Measure");
    	}
    	run("Close All");
    	roiManager("Reset")
	}
}

// Save results 
saveAs("Results", resDir +  "Intensity_measurements.csv");

print("DONE");
