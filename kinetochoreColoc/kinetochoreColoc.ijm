/*
 * Macro to process a folder (with subfolders) of D3D.dv files
 * The process step finds kinetochores (3D blobs) in cells if they
 * are close to the chromosomes. A second step finds 3D blobs of another
 * channel. XYZ positions of both blob sets are saved in two files together
 * with the log file to show how the blobs were detected.
 */

#@ File (label = "Input directory:", style = "directory") input
#@ File (label = "Output directory:", style = "directory") output

// housekeeping
setTheOC();
run("Clear Results");

setBatchMode(true);
inputStr = input;
processFolder(inputStr);
setBatchMode(false);


function setTheOC()	{
	string3doc = "volume surface nb_of_obj._voxels nb_of_surf._voxels"
	string3doc += " integrated_density mean_gray_value std_dev_gray_value";
	string3doc += " median_gray_value minimum_gray_value maximum_gray_value";
	string3doc += " centroid mean_distance_to_surface std_dev_distance_to_surface";
	string3doc += " median_distance_to_surface centre_of_mass bounding_box";
	string3doc += " close_original_images_while_processing_(saves_memory)";
	string3doc += " dots_size=5 font_size=10 show_numbers white_numbers";
	string3doc += " store_results_within_a_table_named_after_the_image_(macro_friendly)";
	run("3D OC Options", string3doc + " redirect_to=none");
}

function processFolder(inputStr) {
	list = getFileList(inputStr);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(inputStr + File.separator + list[i]))
			processFolder(inputStr + File.separator + list[i]);
		if(endsWith(list[i], "D3D.dv")) {
			processFile3DOC(inputStr, list[i]);
		}
	}
}

function processFile3DOC(inputStr, file)	{
	//open(input + File.separator + file);
	fPath = inputStr + file;
	s = "open=[" + fPath + "] autoscale color_mode=Grayscale rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT";
	run("Bio-Formats Importer", s);
	imgName = getTitle();
	getDimensions(width, height, channels, slices, frames);
	if(slices > 5) {
		// make a mask of the chromosomes
		run("Duplicate...", "duplicate channels=1");
		rename("dapi");
		selectWindow("dapi");
		setAutoThreshold("Otsu dark");
		run("Convert to Mask", "method=Otsu background=Dark");
		run("Dilate", "stack");
		run("Dilate", "stack");
		run("Divide...", "value=255 stack");
		// now extract the kinetochores and find a level for 3DOC
		selectWindow(imgName);
		run("Duplicate...", "duplicate channels=4");
		rename("kts");
			midStack = floor(slices / 2);
			setSlice(midStack);
			setAutoThreshold("MaxEntropy dark");
			getThreshold(lower,upper);
		// create a stack where the kts within the chromosome mask are shown
		imageCalculator("Multiply create stack", "kts","dapi");
		selectWindow("Result of kts");
		rename("ktsMod");
		close("dapi");
		close("kts");
		selectWindow("ktsMod");
		stringOC = "threshold=" + lower + " slice=" + midStack;
		stringOC += " min.=100 max.=2000 exclude_objects_on_edges statistics summary";
		run("3D Objects Counter", stringOC);
		close("ktsMod");
		// copy knstrn channel
		run("Duplicate...", "duplicate channels=3");
		rename("knstrn");
			setSlice(midStack);
			setAutoThreshold("MaxEntropy dark");
			getThreshold(lower,upper);
		selectWindow("knstrn");
		stringOC = "threshold=" + lower + " slice=" + midStack;
		stringOC += " min.=100 max.=2000 exclude_objects_on_edges statistics summary";
		run("3D Objects Counter", stringOC);
		close("knstrn");
		close(imgName);
		// save the results
		baseName = replace(imgName, ".dv", "");
		selectWindow("Statistics for ktsMod");
		fPath = output + File.separator + baseName + "_ktsMod.csv";
		saveAs("Results", fPath);
		run("Close");
		selectWindow("Statistics for knstrn");
		fPath = output + File.separator + baseName + "_knstrn.csv";
		saveAs("Results", fPath);
		run("Close");
		selectWindow("Log");
		fPath = output + File.separator + baseName + "_log.txt";
		saveAs("Text", fPath);
		run("Close");
	} else {
		run("Close");
	}	
}
