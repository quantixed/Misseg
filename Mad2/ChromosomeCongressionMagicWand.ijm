/*
 * This is a script to process outputs from our older method for segmenting misaligned chromosomes
 * and metaphase plate. A cropped version has been segmented and then ROISETs created using the magic wand(!)
 * for each chromosome and plate. Previously we would then make a csv of each roiset and process in igor.
 * Now we want to do something similar to the ilastik workflow and read from the rois to measure fluorescence
 * as well as output the coords.
 */

macro "Process chromosomes" {
	setOption("ExpandableArrays", true);
	print("\\Clear");
	dir1 = getDirectory("Original Image Directory ");
	dir2 = getDirectory("Chromosome Outputs Directory ");
	list1 = getFileList(dir1);

	// list of (processed) tiffs one folder deep
	tifPathList = processFolder(dir2);
	tifNameList = newArray(tifPathList.length);
	
	for (i = 0; i < tifPathList.length; i++) {
		tifNameList[i] = File.getNameWithoutExtension(tifPathList[i]);
	}
	
	setBatchMode(true);
	
	// Identify .tif files
	origTifList = newArray;
	
	for (i = 0; i < list1.length; i ++) {		
 		if (endsWith(toLowerCase(list1[i]), ".tif")) {
 			origTifList[origTifList.length] = list1[i];
 		}
	}

	list1 = origTifList;
	lookup = newArray(list1.length);
	match = newArray(list1.length);
	
	for (i = 0; i < list1.length; i ++) {		
 		tifname = list1[i];
		rootstr = File.getNameWithoutExtension(dir1 + tifname) + "_tto"; // ensure we find the segmented version
 		match[i] = 0;
 		// find the segmented version
 		for (j = 0; j < tifNameList.length; j ++) {
 			if (tifNameList[j] == rootstr) {
 				lookup[i] = j;
 				// rather than break we will check for multiple (or no) matches
 				match[i] += 1;
 			}
 		}
 	}

	// now send original and matching processed tif to function
	// roisets can be found using the directory of processed tif from path
	// report if there were multiple (or no) matches
 	for (i = 0; i < list1.length; i ++) {		
		if(match[i] == 0) {
			print("no match for " + list1[i]);
		} else if (match[i] > 1) {
			print("multiple matches for " + list1[i]);
			processOutput(dir1 + list1[i],tifPathList[lookup[i]]);
		} else {
			print("match for " + list1[i]);
			processOutput(dir1 + list1[i],tifPathList[lookup[i]]);
		}
 	}
 	
 	setBatchMode(false);
}

function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	tiflist = newArray;
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i])) {
			sublist = getFileList(input + File.separator + list[i]);
			for (j = 0; j < sublist.length; j ++) {
 				if (endsWith(sublist[j], ".tif")) {
 					tiflist[tiflist.length] = input + list[i] + sublist[j];
 				}
 			}
		}
	}
	
	return tiflist;
}

function processOutput(originalPath,processedPath) {
	open(processedPath);
	rename("processed");
	run("Select None");

	open(originalPath);
	rename("original");
	run("Select None");
	
	selectWindow("processed");
	dir = getInfo("image.directory"); // has file.separator
	list = findROISets(dir);
	fName = getInfo("image.filename");
	fPath = dir + fName;
	resultBaseName = File.getNameWithoutExtension(fPath);
	fPath = dir + resultBaseName; // now set this as image path without extension (for appending)
	
	selectWindow("processed");
	run("Set Measurements...", "area mean standard min integrated stack redirect=None decimal=3");
	
	for(i = 0; i < list.length; i ++) {
		// error check
		if(list[i] == 0) {
			print("Error");
			Array.print(list);
			break;
		}
		retrieveData(fPath,list[i],i);
	}
	run("Close All");
	
	return 0;
}

function retrieveData(fPath,roiPath, ii) {
	selectWindow("processed");
	roiManager("reset");
	run("Clear Results");
	roiManager("open", roiPath);
	sliceList = newArray;
	frameList = newArray;
	
	if (File.exists(fPath + ii + "_coords.csv")) {
		File.delete(fPath + ii + "_coords.csv");
	}
	f = File.open(fPath + ii + "_coords.csv");
	print(f, "roi,obj,xcoord,ycoord,zpos,tpos");

	nROI = roiManager("count");
	
	for(i = 0; i < nROI; i ++) {
		roiManager("select", i);
		run("Measure");
		Roi.getCoordinates(xpoints, ypoints);
		sliceList[i] = getResult('Slice', i);
		frameList[i] = getResult('Frame', i);
		
		for (j = 0; j < xpoints.length; j++) {
			print(f, i + "," + ii + "," + xpoints[j] + "," + ypoints[j] + "," + sliceList[i] + "," + frameList[i]);
		}
	}
	
	File.close(f);
	
	selectWindow("original");
	Stack.getDimensions(width, height, channels, slices, frames);
	nROI = roiManager("count");
	run("Clear Results");
	if (nROI == 0) {
		run("Close All");
		return -1;
	}
	
	f = File.open(fPath + ii + "_orig.csv");
	print(f, "roi,obj,area,mean,max,ch,slice,frame");
	
	for(i = 0; i < nROI; i ++) {
		roiManager("select", i);
		
		for(j = 0; j < channels; j ++) {
			Stack.setPosition(j + 1, sliceList[i], frameList[i]);
			print(f, i + "," + ii + "," + getValue("Area") + "," + getValue("Mean") + "," + getValue("Max") + "," + (j + 1) + "," + getValue("Slice") + "," + getValue("Frame"));
		}
	}
	File.close(f);
	run("Clear Results");
	
	f = File.open(fPath + ii + "_enlg.csv");
	print(f, "roi,obj,area,mean,max,ch,slice,frame");
	
	for(i = 0; i < nROI; i ++) {
		roiManager("select", i);
		run("Enlarge...", "enlarge=3 pixel");
		
		for(j = 0; j < channels; j ++) {
			Stack.setPosition(j + 1, sliceList[i], frameList[i]);
			print(f, i + "," + ii + "," + getValue("Area") + "," + getValue("Mean") + "," + getValue("Max") + "," + (j + 1) + "," + getValue("Slice") + "," + getValue("Frame"));
		}
	}
	File.close(f);
	run("Clear Results");

	return 0;
}

function findROISets(dir) {
	list = getFileList(dir);
	roilist = newArray;
	for (i = 0; i < list.length; i++) {
		name = list[i];
		if (endsWith(name, ".zip")) {
			roilist[roilist.length] = dir + name;
		}
	}
	// we now have a list of roisets only
	arrangedlist = newArray(roilist.length);
	chrCount = 0;
	for (i = 0; i < roilist.length; i++) {
		name = roilist[i];
		if (indexOf(toLowerCase(name), "chr") > 0) {
			chrCount += 1;
			arrangedlist[chrCount] = roilist[i];
		} else {
			arrangedlist[0] = roilist[i];
		}
	}
	return arrangedlist;
}
