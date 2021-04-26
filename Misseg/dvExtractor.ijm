/*
 * Deconvolved dv2 images as stacks, multiple channels.
 * We want to extract to TIFF for analysis in IMOD and Igor
 */

macro "Save dv As Single Channel Tiff" {
	if (nImages > 0) exit ("Please close any images and try again.");
	// make dialog to get choices
	Dialog.create("dv2 to TIFF conversion");
	Dialog.addMessage("Single channel export:");
	Dialog.addCheckbox("Channel 1", false);
	Dialog.addCheckbox("Channel 2", false);
	Dialog.addCheckbox("Channel 3", false);
	Dialog.addCheckbox("Channel 4", false);
	Dialog.addCheckbox("All channels as one hyperstack", false);
	Dialog.show();
	doCh1 = Dialog.getCheckbox();
	doCh2 = Dialog.getCheckbox();
	doCh3 = Dialog.getCheckbox();
	doCh4 = Dialog.getCheckbox();
	hStack = Dialog.getCheckbox();
	// choices taken, now find locations
	dir1 = getDirectory("Choose Source Directory ");
	dir2 = getDirectory("Choose Destination Directory ");
	list = getFileList(dir1);

	dvnum = 0;
	// How many d3d.dv files do we have? Directory could contain other directories.
	for (i=0; i<list.length; i++) {
		if (indexOf(toLowerCase(list[i]), "d3d.dv")>0) {
			dvnum=dvnum+1;
		}
	}
	tifflist = newArray(dvnum);
	dvlist = newArray(dvnum);
	j = 0;
	for (i=0; i<list.length; i++) {
		if (indexOf(toLowerCase(list[i]), "d3d.dv")>0) {
			dvlist[j]=list[i];
			tifflist[j]=replace(list[i],".dv",".tif");
			j=j+1;
		}
	}
	// set up folders in destination
	if(doCh1 == true) File.makeDirectory(dir2+"ch1"+File.separator);
	if(doCh2 == true) File.makeDirectory(dir2+"ch2"+File.separator);
	if(doCh3 == true) File.makeDirectory(dir2+"ch3"+File.separator);
	if(doCh4 == true) File.makeDirectory(dir2+"ch4"+File.separator);
	if(hStack == true) File.makeDirectory(dir2+"all"+File.separator);
	setBatchMode(true);
	for (i=0; i<dvlist.length; i++) {
		showProgress(i+1, list.length);
		s = "open=["+dir1+dvlist[i]+"] autoscale color_mode=Grayscale rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT";
		run("Bio-Formats Importer", s);
		id0 = getImageID();
		if(doCh1 == true)	{
			specificDir = dir2+"ch1"+File.separator;
			run("Duplicate...", "duplicate channels=1");
			id1 = getImageID();
			selectImage(id1);
			saveAs("tiff", specificDir+tifflist[i]);
			close();
			selectImage(id0);
		}
		if(doCh2 == true)	{
			specificDir = dir2+"ch2"+File.separator;
			run("Duplicate...", "duplicate channels=2");
			id1 = getImageID();
			selectImage(id1);
			saveAs("tiff", specificDir+tifflist[i]);
			close();
			selectImage(id0);
		}
		if(doCh3 == true)	{
			specificDir = dir2+"ch3"+File.separator;
			run("Duplicate...", "duplicate channels=3");
			id1 = getImageID();
			selectImage(id1);
			saveAs("tiff", specificDir+tifflist[i]);
			close();
			selectImage(id0);
		}
		if(doCh4 == true)	{
			specificDir = dir2+"ch4"+File.separator;
			run("Duplicate...", "duplicate channels=4");
			id1 = getImageID();
			selectImage(id1);
			saveAs("tiff", specificDir+tifflist[i]);
			close();
			selectImage(id0);
		}
		if(hStack == true)	{
			specificDir = dir2+"all"+File.separator;
			selectImage(id0);
			saveAs("tiff", specificDir+tifflist[i]);
		}
		selectImage(id0);
		close();
	}
	setBatchMode(false);
}

macro "Save dv As MultiChannel Tiff" {
	if (nImages > 0) exit ("Please close any images and try again.");
	// get choices for RGBW export
	fourChannel = newArray("Ch1","Ch2","Ch3","Ch4","*None*");
	Dialog.create("dv2 to TIFF conversion");
	Dialog.addMessage("Multi-channel RGB export:");
	Dialog.addChoice("Red", fourChannel, "*None*");
	Dialog.addChoice("Green", fourChannel, "*None*");
	Dialog.addChoice("Blue", fourChannel, "*None*");
	Dialog.addChoice("White", fourChannel, "*None*");
	Dialog.show();
	rCh = Dialog.getChoice();
	gCh = Dialog.getChoice();
	bCh = Dialog.getChoice();
	wCh = Dialog.getChoice();
	// sanity check user choice
	if(rCh == gCh && gCh == bCh && bCh == wCh) exit ("Incorrect choice.");
	// choices taken, now find locations
	dir1 = getDirectory("Choose Source Directory ");
	dir2 = getDirectory("Choose Destination Directory ");
	list = getFileList(dir1);

	dvnum = 0;
	// How many d3d.dv files do we have? Directory could contain other directories.
	for (i=0; i<list.length; i++) {
		if (indexOf(toLowerCase(list[i]), "d3d.dv")>0) {
			dvnum=dvnum+1;
		}
	}
	tifflist = newArray(dvnum);
	dvlist = newArray(dvnum);
	j = 0;
	for (i=0; i<list.length; i++) {
		if (indexOf(toLowerCase(list[i]), "d3d.dv")>0) {
			dvlist[j]=list[i];
			tifflist[j]=replace(list[i],".dv",".tif");
			j=j+1;
		}
	}
	// set up folders in destination
	setBatchMode(true);
	for (i=0; i<dvlist.length; i++) {
		showProgress(i+1, list.length);
		s = "open=["+dir1+dvlist[i]+"] autoscale color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT";
		run("Bio-Formats Importer", s);
		id0 = getImageID();
		// split channels
		run("Split Channels");
		// predicted names for the split channels
		c1 = "[C1-"+dvlist[i]+"]";
		c2 = "[C2-"+dvlist[i]+"]";
		c3 = "[C3-"+dvlist[i]+"]";
		c4 = "[C4-"+dvlist[i]+"]";
		// convert the RGBW choice into the correct named windows - with [] formatting
		mergeString = "red="+rCh+" green="+gCh+" blue="+bCh+" gray="+wCh+" keep ignore";
		mergeString = replace(mergeString,"Ch1",c1);
		mergeString = replace(mergeString,"Ch2",c2);
		mergeString = replace(mergeString,"Ch3",c3);
		mergeString = replace(mergeString,"Ch4",c4);
		run("Merge Channels...", mergeString);
		id1 = getImageID();
		selectImage(id1);
		saveAs("tiff", dir2+tifflist[i]);
		close();
		// close all open windows
		imgArray = newArray(nImages);
		for (i=0; i<nImages; i++)	{
			selectImage(i+1);
			close();
		}
	}
	setBatchMode(false);
}

macro "Save A Segmented Channel From a dv" {
	if (nImages > 0) exit ("Please close any images and try again.");
	// make dialog to get choices
	Dialog.create("dv2 to TIFF conversion");
	Dialog.addMessage("Single channel export:");
	Dialog.addCheckbox("Channel 1", false);
	Dialog.addCheckbox("Channel 2", false);
	Dialog.addCheckbox("Channel 3", false);
	Dialog.addCheckbox("Channel 4", false);
	Dialog.show();
	doCh1 = Dialog.getCheckbox();
	doCh2 = Dialog.getCheckbox();
	doCh3 = Dialog.getCheckbox();
	doCh4 = Dialog.getCheckbox();
	// choices taken, now find locations
	dir1 = getDirectory("Choose Source Directory ");
	dir2 = getDirectory("Choose Destination Directory ");
	list = getFileList(dir1);

	dvnum = 0;
	// How many d3d.dv files do we have? Directory could contain other directories.
	for (i=0; i<list.length; i++) {
		if (indexOf(toLowerCase(list[i]), "d3d.dv")>0) {
			dvnum=dvnum+1;
		}
	}
	tifflist = newArray(dvnum);
	dvlist = newArray(dvnum);
	j = 0;
	for (i=0; i<list.length; i++) {
		if (indexOf(toLowerCase(list[i]), "d3d.dv")>0) {
			dvlist[j]=list[i];
			tifflist[j]=replace(list[i],".dv",".tif");
			j=j+1;
		}
	}
	// set up folders in destination
	if(doCh1 == true) File.makeDirectory(dir2+"ch1"+File.separator);
	if(doCh2 == true) File.makeDirectory(dir2+"ch2"+File.separator);
	if(doCh3 == true) File.makeDirectory(dir2+"ch3"+File.separator);
	if(doCh4 == true) File.makeDirectory(dir2+"ch4"+File.separator);
	setBatchMode(true);
	for (i=0; i<dvlist.length; i++) {
		showProgress(i+1, list.length);
		s = "open=["+dir1+dvlist[i]+"] autoscale color_mode=Grayscale rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT";
		run("Bio-Formats Importer", s);
		id0 = getImageID();
		if(doCh1 == true)	{
			specificDir = dir2+"ch1"+File.separator;
			run("Duplicate...", "duplicate channels=1");
			id1 = getImageID();
			selectImage(id1);
			doTheSegmentation();
			saveAs("tiff", specificDir+tifflist[i]);
			close();
			selectImage(id0);
		}
		if(doCh2 == true)	{
			specificDir = dir2+"ch2"+File.separator;
			run("Duplicate...", "duplicate channels=2");
			id1 = getImageID();
			selectImage(id1);
			doTheSegmentation();
			saveAs("tiff", specificDir+tifflist[i]);
			close();
			selectImage(id0);
		}
		if(doCh3 == true)	{
			specificDir = dir2+"ch3"+File.separator;
			run("Duplicate...", "duplicate channels=3");
			id1 = getImageID();
			selectImage(id1);
			doTheSegmentation();
			saveAs("tiff", specificDir+tifflist[i]);
			close();
			selectImage(id0);
		}
		if(doCh4 == true)	{
			specificDir = dir2+"ch4"+File.separator;
			run("Duplicate...", "duplicate channels=4");
			id1 = getImageID();
			selectImage(id1);
			doTheSegmentation();
			saveAs("tiff", specificDir+tifflist[i]);
			close();
			selectImage(id0);
		}
		selectImage(id0);
		close();
	}
	setBatchMode(false);
}

function doTheSegmentation() {
	// previous method
//	run("Gaussian Blur...", "sigma=2 stack");
//	setAutoThreshold("Otsu dark");
//	run("Convert to Mask", "method=Otsu background=Default calculate");
	// new method
	run("Gamma...", "value=3 stack");
//	run("Gaussian Blur...", "sigma=2 stack");
	run("Convert to Mask", "method=MaxEntropy background=Dark calculate");
}

macro "Get Measurements From All dv Files" {
	if (nImages > 0) exit ("Please close any images and try again.");
	// choices taken, now find locations
	dir1 = getDirectory("Choose Source Directory ");
	list = getFileList(dir1);

	dvnum = 0;
	// How many d3d.dv files do we have? Directory could contain other directories.
	for (i=0; i<list.length; i++) {
		if (indexOf(toLowerCase(list[i]), "d3d.dv")>0) {
			dvnum=dvnum+1;
		}
	}
	tifflist = newArray(dvnum);
	dvlist = newArray(dvnum);
	j = 0;
	for (i=0; i<list.length; i++) {
		if (indexOf(toLowerCase(list[i]), "d3d.dv")>0) {
			dvlist[j]=list[i];
			j=j+1;
		}
	}
	setBatchMode(true);
	for (i=0; i<dvlist.length; i++) {
		showProgress(i+1, list.length);
		s = "open=["+dir1+dvlist[i]+"] autoscale color_mode=Grayscale rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT";
		run("Bio-Formats Importer", s);
		id0 = getImageID();
		makeLine(0, 1, 0, 0);
		run("Measure");
		close();
	}
	setBatchMode(false);
}
