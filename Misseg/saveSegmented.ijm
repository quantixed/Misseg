/*
 * This macro is for making segmented images
 * It was designed so that a user would segment each frame of a stack
 * saving the ROIs to the manager. This macro makes a copy of the stack
 * with all values as 255 and then fills all ROIs with 0
 */

Stack.getDimensions(width, height, channels, slices, frames);
win = getTitle();
dir = getDirectory("image");
close();
outputDir = dir + "output" + File.separator;
if (!File.exists(outputDir))
	File.makeDirectory(outputDir);
roiSet = replace(win, ".tif", ".zip");
roiManager("Save", outputDir + roiSet);
newName = replace(win, ".tif", "_a.tif");
//newImage(newName, "8-bit grayscale-mode", width, height, channels, slices, frames);
newImage(newName, "8-bit white", width, height, slices);
selectWindow(newName);
setColor(0);
roiManager("Fill"); // fill 0 in a 255 stack
// save with original name in output directory
saveAs("TIFF", outputDir + win);
close();
roiManager("reset");