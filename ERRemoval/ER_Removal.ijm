/*
 * This script will take a directory of tiffs where each has three channels, multi z, multi t
 * Channel 1 is ER, Channel 2 is PM (Channel 3 is DNA and not used).
 * It will measure the cell area in each section over time and the segmented ER within that area
 * this is as a single line in Results which needs further calculation to get the volumes at each t
 */

macro "ER Volume Measurements"	{
	dir = getDirectory("Choose a Directory ");
	list = getFileList(dir);
	setBatchMode(true);
	
	for (i = 0; i < list.length; i ++)	{
		inputPath = dir+list[i];
		if (endsWith(inputPath, ".tif") || endsWith(inputPath, ".tiff")) {
			if (endsWith(inputPath, "mask.tif")) {
				continue;
			}
			open(inputPath);
			segmentAndMeasure();
		}
	}
	setBatchMode(false);
}

function segmentAndMeasure()	{
	// 2D images are 15.3930 pixels per micron
	// 0.06496 um per pixel which is 1 px = 0.00422 um^2
	// each section is 2 um. So voxel is 0.00844 um^3
	
	inputId = getImageID(); // get active image 
	inputTitle = getTitle();
	dir = getInfo("image.directory");
	fName = getInfo("image.filename");
	fPath = dir + File.separator + fName;
	resultBaseName = File.getNameWithoutExtension(fPath);
	fPath = dir + File.separator + resultBaseName; // now set this as image path without extension (for appending)
	
	run("Duplicate...", "title=tempImg duplicate channels=1-2");
	tempTitle = "tempImg";
	// close original image
	selectImage(inputId);
	close();
	
	selectWindow(tempTitle);
	Stack.getDimensions(width, height, channels, slices, frames); 
	setBatchMode("hide");
	// do CLAHE on all frames all slices all channels
	for (i=1; i<frames+1; i++){ 
	        //Stack.setFrame(i); 
	        for (j=1; j<slices+1; j++) { 
	                //Stack.setSlice(j);
	                for (k=1; k<channels+1; k++) {
	                	Stack.setPosition(k,j,i);
	                	run("Enhance Local Contrast (CLAHE)", "blocksize=127 histogram=256 maximum=2 mask=*None*"); 
	                }
	        } 
	}
	selectWindow(tempTitle);
	run("Duplicate...", "title=tempCh1 duplicate channels=1");
	selectWindow(tempTitle);
	run("Duplicate...", "title=tempCh2 duplicate channels=2");
	
	// select stargazin channel, find mid point
	selectWindow("tempCh2");
	Stack.setFrame(floor(frames / 2));
	Stack.setSlice(floor(slices / 2));
	setAutoThreshold("Shanbhag dark no-reset");
	run("Convert to Mask", "method=Shanbhag background=Dark");
	// get rid of stray pixels
	run("Open", "stack");
	
	// select stargazin channel, find mid point
	selectWindow("tempCh1");
	Stack.setFrame(floor(frames / 2));
	Stack.setSlice(floor(slices / 2));
	setAutoThreshold("Shanbhag dark no-reset");
	run("Convert to Mask", "method=Shanbhag background=Dark");
	
	run("Set Measurements...", "area mean standard min integrated stack redirect=tempCh1 decimal=3");
	run("Clear Results");
	for (i=1; i<frames+1; i++){
		selectWindow("tempCh1");
		Stack.setFrame(i);
		selectWindow("tempCh2");
		Stack.setFrame(i);
		for (j=1; j<slices+1; j++) {
			selectWindow("tempCh1");
			Stack.setSlice(j);
			selectWindow("tempCh2");
			Stack.setSlice(j);
			run("Create Selection");
			if (getValue("selection.size") > 0) {
				run("Convex Hull");
				run("Measure");
			}
			run("Select None");
		}
	}
	
	// reset set measurements
	run("Set Measurements...", "area mean standard min integrated stack redirect=None decimal=3");
	run("Merge Channels...", "c1=tempCh2 c2=tempCh1 create");
	save(fPath + "_mask.tif");
	saveAs("results", fPath + "_results.txt");
	run("Close All");
	setBatchMode("show");
}