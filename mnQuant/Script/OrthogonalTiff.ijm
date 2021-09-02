macro "Make Ortho Composite"	{
	if (nImages > 1) exit ("Use a single image or single montage");
	if (nImages == 0)	exit("No image open");
	
	imageID = getImageID();
	title = getTitle();
	dir = getDirectory("image");
	getDimensions(width, height, channels, slices, frames);

	xx = 184;
	yy = 42;
	zz = 10;
	grout = 4;

	Dialog.create("Specify Orthogonal view");
	Dialog.addMessage("Specify XYZ cuts");
	Dialog.addNumber("X (px)", xx);
	Dialog.addNumber("Y (px)", yy);
	Dialog.addNumber("Z (px)", zz);
	Dialog.addNumber("Grout size in frame (px)", grout);
	Dialog.show();
	xx = Dialog.getNumber();
	yy = Dialog.getNumber();
	zz = Dialog.getNumber();
	grout = Dialog.getNumber();
	
	selectImage(imageID);
	run("Orthogonal Views");
	Stack.setOrthoViews(xx, yy, zz);
	Stack.getOrthoViewsIDs(XY, YZ, XZ);
	// need to find width (and height) of YZ (XZ)
	selectImage(YZ);
	getDimensions(YZwidth, YZheight, YZchannels, YZslices, YZframes);
	//print(XY,YZ,XZ);
	
	newImage("ortho", "16-bit composite-mode", width + grout + YZwidth, height + grout + YZwidth, channels, 1, 1);
	if (bitDepth() == 8) setColor(255);
	if (bitDepth() == 24) setColor(255,255,255);
	if (bitDepth() == 16) setColor(65535);
	run("Select All");
	
	for (i=0; i<channels; i++)	{
		Stack.setChannel(i + 1);
		run("Fill", "slice");
	}
	
	for (i=0; i<channels; i++)	{
		selectImage(XY);
		Stack.setChannel(i + 1);
		getMinAndMax(min, max);
		run("Copy");
		selectWindow("ortho");
		Stack.setChannel(i + 1);
		makeRectangle(0, 0, width, height);
		run("Paste");
	
		selectImage(YZ);
		run("Copy");
		selectWindow("ortho");
		Stack.setChannel(i + 1);
		makeRectangle(width + grout, 0, YZwidth, height);
		run("Paste");
	
		selectImage(XZ);
		run("Copy");
		selectWindow("ortho");
		Stack.setChannel(i + 1);
		makeRectangle(0, height + grout, width, YZwidth);
		run("Paste");
		setMinAndMax(min, max);
	}
	selectImage(imageID);
	Stack.stopOrthoViews;
	close();

	// now burn on the lines hard coded as 10 x 2
	selectWindow("ortho");
	run("Select None");
	depth = floor(YZwidth * (zz / slices));
	for (i=0; i<channels; i++)	{
		Stack.setChannel(i + 1);
		// add ticks, left, right, top, bottom
		// XY
		fillRect(0, yy - 1, 10, 2);
		fillRect(width - 10, yy - 1, 10, 2);
		fillRect(xx - 1, 0, 2, 10);
		fillRect(xx - 1, height - 10, 2, 10);
		// YZ
		fillRect(width + grout, yy - 1, 10, 2);
		fillRect(width + grout + YZwidth - 10, yy - 1, 10, 2);
		fillRect(width + grout + depth, 0, 2, 10);
		fillRect(width + grout + depth, height - 10, 2, 10);
		// XZ
		fillRect(0, height + grout + depth, 10, 2);
		fillRect(width - 10, height + grout + depth, 10, 2);
		fillRect(xx - 1, height + grout, 2, 10);
		fillRect(xx - 1, height + grout + YZwidth - 10, 2, 10);
		
	}
	if (bitDepth() == 8) setColor(0);
	if (bitDepth() == 24) setColor(0,0,0);
	if (bitDepth() == 16) setColor(0);
	for (i=0; i<channels; i++)	{
		Stack.setChannel(i + 1);
		fillRect(width, height + grout + depth, grout + depth + 2, 2);
		fillRect(width + grout + depth, height, 2, grout + depth + 2);
	}
	selectWindow("ortho");
	run("Select None");
	save(dir + title + "_ortho.tif");
	close();
}