/*
 * Macro to process a folder (with subfolders) of D3D.dv files
 * 
 * This will save out an RGB stack (24-bit) of each file in the destination directory
 * Assumes 4-channel dv stack to start with
 * 
 */

#@ File (label = "Input directory:", style = "directory") input
#@ File (label = "Output directory:", style = "directory") output

if (nImages > 0) exit ("Please close any images and try again.");
// get choices for RGB export
fourChannel = newArray("Ch1","Ch2","Ch3","Ch4","*None*");
Dialog.create("dv2 to TIFF conversion");
Dialog.addMessage("Multi-channel RGB export:");
Dialog.addChoice("Red", fourChannel, "*None*");
Dialog.addChoice("Green", fourChannel, "*None*");
Dialog.addChoice("Blue", fourChannel, "*None*");
Dialog.show();
rCh = Dialog.getChoice();
gCh = Dialog.getChoice();
bCh = Dialog.getChoice();

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
			processFile(inputStr, list[i], rCh, gCh, bCh);
		}
	}
}

function processFile(inputStr, file, rCh, gCh, bCh)	{
	//open(input + File.separator + file);
	fPath = inputStr + file;
	s = "open=[" + fPath + "] autoscale color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT";
	run("Bio-Formats Importer", s);
	id0 = getImageID();
	imgName = getTitle();
	getDimensions(width, height, channels, slices, frames);
	if(slices < 5) {
		run("Close");
		return 0;
	}
	baseName = replace(imgName, ".dv", "");
	// split channels
	run("Split Channels");
	// predicted names for the split channels
	c1 = "[C1-"+imgName+"]";
	c2 = "[C2-"+imgName+"]";
	c3 = "[C3-"+imgName+"]";
	c4 = "[C4-"+imgName+"]";
	// convert the RGB choice into the correct named windows - with [] formatting
	mergeString = "red="+rCh+" green="+gCh+" blue="+bCh+" create ignore";
	mergeString = replace(mergeString,"Ch1",c1);
	mergeString = replace(mergeString,"Ch2",c2);
	mergeString = replace(mergeString,"Ch3",c3);
	mergeString = replace(mergeString,"Ch4",c4);
	run("Merge Channels...", mergeString);
	id1 = getImageID();
	selectImage(id1);
	run("RGB Color", "slices");
	fPath = output + File.separator + baseName + ".tif";
	saveAs("Tiff", fPath);
	close();
	// close all open windows
	imgArray = newArray(nImages);
	for (i=0; i<nImages; i++)	{
		selectImage(i+1);
		close();
	}
}