/*
 * Goal is to segment nuclei (in channel 1) in 3D and then quantify the fluorescence
 * of all channels in the 3D objects that were segmented.
 * 
 * Input: folder of multi-channel z-stacks
 * Process:	segment channel 1 and find 3D objects
 * 			then quantify the fluorescence in each and save csv of results for each channel
 * 			Background is automatically found and measured. 3D object map is also saved for verification.
 * Requires: installation of 3D Image Suite and ImageScience
 * Credit: Stephen Royle with initial help from Laura Cooper.
 */

// first macro definition
macro "Nuclei Segmentation 3D" {
	if (nImages > 0) exit ("Please close any images and try again.");
	inDir = getDirectory("Source Directory ");
	outDir = getDirectory("Destination Directory ");

	setBatchMode(true);
	processFolder(inDir, outDir);
	setBatchMode(false);
}

// builds file list
function processFolder(inputStr, outputStr) {
	list = getFileList(inputStr);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(endsWith(list[i], ".tif")) {
			processFile(inputStr, list[i], outputStr);
		}
	}
}

function processFile(inputStr, file, outputStr)	{
	fPath = inputStr + File.separator + file;
	open(fPath);

	title = getTitle();
	fname = File.nameWithoutExtension;
	selectWindow(title);
	dir = outputStr;
	if(checkImageIsOK(title) == -1) {
		run("Close All");
		return -1;
	}
	saveBackgroundValues(title, fname, dir);
	run("Duplicate...", "title=C1 duplicate channels=1");
	run("Median...", "radius=2 stack");
	setAutoThreshold("Li dark no-reset");
	setOption("BlackBackground", false);
	run("Convert to Mask", "method=Li background=Dark calculate");
	run("Fill Holes", "stack");
	saveAs("Tiff", dir + "mask_" + fname + ".tif");
	rename("C1");
	//Analyse Particles 
	run("3D Objects Counter", "threshold=128 slice=7 min.=500 max.=2240000 objects summary");
	//Start ROI manager
	run("3D Manager");
	Ext.Manager3D_SelectAll();
	Ext.Manager3D_Delete();
	Ext.Manager3D_AddImage();
	Ext.Manager3D_SelectAll();
	
	//Measure regions
	selectWindow(title);
	run("Split Channels");
	
	for (i = 1; i <= 4; i++) {
		chan="C"+i+"-";
		selectWindow(chan+title);
		Ext.Manager3D_Quantif();
		Ext.Manager3D_SaveResult("Q",dir + "quant_" + chan + fname + ".csv");
		Ext.Manager3D_CloseResult("Q");
	}
	selectWindow("Objects map of C1");
	saveAs("Tiff", dir + "obj_" + fname + ".tif");
	run("Close All");
}


function saveBackgroundValues(title, fname, dir) {
	selectWindow(title);
	Stack.getDimensions(width, height, channels, slices, frames);
	run("Duplicate...", "title=bg duplicate");
	currentWin = "bg";
	sumVar = 0;
	if (channels > 1) sumVar += 1;
	if (slices > 1) sumVar += 1;
	if (channels > 1) sumVar += 1;
	for (i = 0; i < sumVar - 1; i++) {
		run("Z Project...", "projection=[Sum Slices]");
		currentWin = getTitle();
	}
	selectWindow(currentWin);
	// find bg ROI
	inx = 10;
	iny = 10;
	boxw = 20;
	boxh = 20;
	for (i = 0; i < 9; i++) {
		ltx = getLTcoord("x",i,inx,boxw,width);
		lty = getLTcoord("y",i,inx,boxh,height);
		makeRectangle(ltx, lty, boxw, boxh);
		getStatistics(area, mean, min, max, std, histogram);
		if (i == 0) minMean = mean;
		minMean = minOf(minMean,mean);
		if (minMean == mean) best = i;
	}
	run("Clear Results");
	selectWindow("bg");
	run("Z Project...", "projection=[Average Intensity]");
	ltx = getLTcoord("x",best,inx,boxw,width);
	lty = getLTcoord("y",best,inx,boxh,height);
	for (i = 0; i < channels; i++) {
		Stack.setChannel(i + 1);
		makeRectangle(ltx, lty, boxw, boxh);
		run("Measure");
	}
	saveAs("Results", dir + "bg_" + fname + ".csv");
	// close background windows and leave original open
	selectWindow(title);
	close("\\Others");
}

function getLTcoord(XorY,ii, inset, box, img) {
	coord = 0;
	if(XorY == "x")	{
		if (ii%3 == 0) coord = inset;
		if (ii%3 == 1) coord = floor((img / 2) - (box / 2));
		if (ii%3 == 2) coord = img - inset - box;
	} else {
		if (floor(ii/3) == 0) coord = inset;
		if (floor(ii/3) == 1) coord = floor((img / 2) - (box / 2));
		if (floor(ii/3) == 2) coord = (img - inset) - box;
	}
	return coord;
}

function checkImageIsOK(winTitle) {
	okVar = 0;
	selectWindow(winTitle);
	Stack.getDimensions(width, height, channels, slices, frames);
	if(channels != 4 || slices == 1) okVar = -1;
	
	return okVar;
}