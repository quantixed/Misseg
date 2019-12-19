/*
 * Macro to process a folder (with subfolders) of D3D.dv files
 * The process step tests all automated thesholding procedures
 * for the files and saves a montage for each. Channels 4 and 3
 * are tested.
 */

#@ File (label = "Input directory:", style = "directory") input
#@ File (label = "Output directory:", style = "directory") output

// housekeeping
run("Clear Results");

setBatchMode(true);
inputStr = input;
processFolder(inputStr);
setBatchMode(false);

function processFolder(inputStr) {
	list = getFileList(inputStr);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(inputStr + File.separator + list[i]))
			processFolder(inputStr + File.separator + list[i]);
		if(endsWith(list[i], "D3D.dv")) {
			processFile(inputStr, list[i]);
		}
	}
}

function processFile(inputStr, file)	{
	//open(input + File.separator + file);
	fPath = inputStr + file;
	s = "open=[" + fPath + "] autoscale color_mode=Grayscale rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT";
	run("Bio-Formats Importer", s);
	imgName = getTitle();
	getDimensions(width, height, channels, slices, frames);
	if(slices < 5) {
		run("Close");
		return 0;
	}
	baseName = replace(imgName, ".dv", "");
	midStack = floor(slices / 2);
	selectWindow(imgName);
	run("Duplicate...", "duplicate channels=4 slices=" + midStack);
	rename("kts");
	run("Auto Threshold", "method=[Try all] white");
	selectWindow("Montage");
	fPath = output + File.separator + baseName + "_ktSeg.tif";
	saveAs("Tiff", fPath);
	close();
	selectWindow("kts");
	close();
	selectWindow(imgName);
	run("Duplicate...", "duplicate channels=3 slices=" + midStack);
	rename("knstrn");
	run("Auto Threshold", "method=[Try all] white");
	selectWindow("Montage");
	fPath = output + File.separator + baseName + "_knstrnSeg.tif";
	saveAs("Tiff", fPath);
	close();
	selectWindow("knstrn");
	close();
	selectWindow(imgName);
	close();
}
