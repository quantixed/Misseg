#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// CSVs from ROI Monager
// Naming of files is important

////////////////////////////////////////////////////////////////////////
// Menu items
////////////////////////////////////////////////////////////////////////
Menu "Macros"
	"Congression Analysis...",  ROIAnalysis()
	"Start Over", CleanSlate()
	End
End

////////////////////////////////////////////////////////////////////////
// Master functions and wrappers
////////////////////////////////////////////////////////////////////////
Function ROIAnalysis()
	PreLoader()
End

// we can stop beween these two if User input needed.
//Function TheLoader()
//	LoadCSVs()
//	ProcessAllModels()
//	PlotSamples(100)
//	FormatPlots()
//	CollectAllMeasurements()
//	MakeTheLayouts("cond",5,3, alphaSort = 1, saveIt = 0)
//	MakeTheLayouts("p",5,3, rev = 1, saveIt = 0)
End

////////////////////////////////////////////////////////////////////////
// Main functions
////////////////////////////////////////////////////////////////////////
Function PreLoader()
	NewPath/O/Q/M="Please find chromosome ROI folder" expDiskFolder1
	if (V_flag != 0)
		DoAlert 0, "Disk folder error"
		return -1
	endif
	String FileList1 = IndexedFile(expDiskFolder1,-1,".csv")
	if(ItemsInList(FileList1) == 0)
		DoAlert 0, "No point (csv) files found"
		return -1
	endif
	Wave/T fileName1Wave = ListToTextWave(FileList1,";")
	MoveWave fileName1Wave, root:fileName1Wave
	
	NewPath/O/Q/M="Please find plate ROI folder" expDiskFolder2
	if (V_flag != 0)
		DoAlert 0, "Disk folder error"
		return -1
	endif
	String FileList2 = IndexedFile(expDiskFolder2,-1,".csv")
	if(ItemsInList(FileList2) == 0)
		DoAlert 0, "No point (csv) files found"
		return -1
	endif
	Wave/T fileName2Wave = ListToTextWave(FileList2,";")
	MoveWave fileName2Wave, root:fileName2Wave
	
	LoadCSVs(fileName1Wave, "expDiskFolder1", "chr", 9, 2, 3, 0.0649645943)
	LoadCSVs(fileName2Wave, "expDiskFolder2", "plt", 9, 2, 3, 0.0649645943)
	DoCalcs()
End

Function LoadCSVs(fWaveToLoad, pathNameString, folderString, zMax, zStep, tStep, pxSize)
	Wave/T/Z fWaveToLoad
	String pathNameString,folderString
	Variable zMax, zStep, tStep, pxSize
	
	Variable nFiles = numpnts(fWaveToLoad)
	String thisFile
	
	Variable i
	
	// set up folder
	NewDataFolder/O/S $("root:" + folderString)
	String wList = "coords0;coords1;zW;tW;"
 	
	for (i = 0; i < nFiles; i += 1)
		ThisFile = fWaveToLoad[i]
		LoadWave/A=coords/J/O/Q/K=1/L={0,1,0,0,0}/P=$pathNameString ThisFile
		Wave xW = $("coords0")
		Make/O/N=(DimSize(xW,0))/D zW, tW // make another 2 columns
		zW = mod(str2num(ThisFile[0,3]) - 1, zMax) * zStep
		tW = floor((str2num(ThisFile[0,3]) - 1) / zMax) * tStep
		Concatenate/O/NP=1/KILL wList, $("xy_" + num2str(i))
		Wave matA = $("xy_" + num2str(i))
		matA[][0,1] *= pxSize
	endfor
	
	Concatenate/O/NP=0 WaveList("xy*",";",""), $("root:"+folderString+"_mat")
	SetDataFolder root:
End

Function DoCalcs()
	WAVE/Z chr_mat
	WAVE/Z plt_mat
	
	Variable tStep = 30
	Variable nFrames = 41
	
	Make/O/N=(nFrames)/D c1c2distance, c1edgeDistance, timeW = p * 3
	
	Variable i
	
	for(i = 0; i < nFrames; i += 1)
		Wave tempC = XYZatT(chr_Mat, timeW[i])
		Wave tempP = XYZatT(plt_Mat, timeW[i])
		if(sum(tempC) == 0 && sum(tempP) == 0)
			c1c2distance[i] = NaN
			c1edgeDistance[i] = NaN
		elseif(sum(tempC) == 0)
			c1c2distance[i] = 0
			c1edgeDistance[i] = 0
		else
			c1c2distance[i] = Centroid2Centroid(tempC,tempP)
			c1edgeDistance[i] = Centroid2Edge(tempC,tempP)
		endif
		KillWaves/Z tempC, tempP
	endfor
	
	KillWindow/Z p_dist
	Display/N=p_dist
	AppendToGraph/W=p_dist c1c2distance vs timeW
	AppendToGraph/W=p_dist c1edgedistance vs timeW
	Label/W=p_dist left "Distance (μm)"
	Label/W=p_dist bottom "Time (min)"
	ModifyGraph/W=p_dist lstyle(c1edgeDistance)=3
	Legend/W=p_dist/C/N=text0/J/F=0/X=0.00/Y=0.00 "\\s(c1c2distance) Centroid-Centroid\r\\s(c1edgeDistance) Centroid-Edge"
End

STATIC Function/WAVE XYZatT(matA,tVar)
	Wave matA
	Variable tVar
	
	MatrixOP/O index = col(matA,3)
	index[] = (index[p] == tVar) ? p : NaN 
	Wavetransform zapNaNs index
	if(numpnts(index) == 0)
		Make/O/D/N=(1,3)/FREE resultMat = 0 // no data
	else
		Make/O/D/N=(numpnts(index),3)/FREE resultMat = matA[index[p]][q]
	endif
	
	return resultMat
End

STATIC Function Centroid2Centroid(matA,matB)
	Wave matA, matB
	
	MatrixOP/O/FREE cA = averagecols(matA)
	MatrixOP/O/FREE cB = averagecols(matB)
	
	MatrixOP/O/FREE result = sqrt(sumrows((cA - cB) * (cA - cB)))
	
	return result[0]
End

STATIC Function Centroid2Edge(matA,matB)
	Wave matA, matB
	
	Variable nRows = DimSize(matB,0)
	MatrixOP/O/FREE centroid1 = averagecols(matA)
	MatrixOP/O/FREE cA = rowRepeat(centroid1,nRows)
	
	MatrixOP/O/FREE result = sqrt(sumrows((cA - matB) * (cA - matB)))
	
	return WaveMin(result)
End



//Function MakeObjectContourWaves(id,pxSize)
//	Variable id // this is the row number from fileNameFWave - will be used as a UID
//	Variable pxSize // nm per pixel
//	
//	Concatenate/O/KILL wavelist("wave*",";",""), matA
//	WaveStats/Q/RMD=[][0] matA
//	// Scale the coordinates to real values
//	matA[][2,4] *= pxSize
//	Variable nObjects = V_max + 1 // objects in IMOD is 1-based
//	Variable nContours, contourVar
//	String wName
//	
//	Variable i,j
//	
//	for (i = 0; i < nObjects; i += 1)
//		MatrixOP/O filtObj = col(matA,0)
//		filtObj[] = (filtObj[p] == i) ? matA[p][1] : NaN
//		WaveTransform zapnans filtObj
//		FindDuplicates/RN=uniqueContours filtObj
//		nContours = numpnts(uniqueContours)
//		// zero-indexed list of contours in this object
//		for (j = 0; j < nContours; j += 1)
//			contourVar = uniqueContours[j]
//			// find the rows that correspond to each contour
//			Duplicate/O/FREE matA,matB
//			matB[][] = (matB[p][0] == i && matB[p][1] == contourVar) ? matB[p][q] : NaN
//			MatrixOp/O/FREE xW = col(matB,2)
//			MatrixOp/O/FREE yW = col(matB,3)
//			// no need to take the z column here
//			WaveTransform zapnans xW
//			WaveTransform zapnans yW
//			// Now make ObjectContour waves
//			wName = "vs_" + num2str(id) + "_" + num2str(i) + "_" + num2str(contourVar)
//			Concatenate/O/NP=1 {xW,yW}, $wName
//			// close contour if it's not already closed 
//			Wave cW = $wName // contour wave
//			if(cW[DimSize(cW,0)-1][0] != cW[0][0] || cW[DimSize(cW,0)-1][1] != cW[0][1])
//				InsertPoints DimSize(cW,0), 1, cW
//				cW[DimSize(cW,0)-1][] = cW[0][q]
//			endif
//			// delete if it's just 0-3 pixels?
//			if(DimSize(cW,0) < 4)
//				KillWaves cW
//			endif
//		endfor
//	endfor
//	KillWaves/Z filtObj,UniqueContours,MatA
//End
//
//// this function goes into each datafolder and run some code on the contours in there
//Function ProcessAllModels()
//	SetDataFolder root:data:	// relies on earlier load
//	DFREF dfr = GetDataFolderDFR()
//	String folderName
//	Variable numDataFolders = CountObjectsDFR(dfr, 4)
//	SetUpWindows(numDataFolders)
//	
//	Variable i
//		
//	for(i = 0; i < numDataFolders; i += 1)
//		folderName = GetIndexedObjNameDFR(dfr, 4, i)
//		SetDataFolder ":" + folderName + ":"
//		FindVesicleCentresAndPlotOut()
//		SetDataFolder root:data:
//	endfor
//	SetDataFolder root:
//End
//
//STATIC Function SetUpWindows(nWindows)
//	Variable nWindows
//	
//	String plotName
//	Variable i
//	
//	for(i = 0; i < nWindows; i += 1)
//		plotName = "cond_" + num2str(i) + "_allVsPlot"
//		KillWindow/Z $plotName
//		Display/N=$plotName/HIDE=1
//		plotName = "cond_" + num2str(i) + "_allRotVsPlot"
//		KillWindow/Z $plotName
//		Display/N=$plotName/HIDE=1
//	endfor
//End
//
//Function FindVesicleCentresAndPlotOut()
//	// finding the centre of each vesicle and placing coords into a wave called VsWave
//	String wList = WaveList("Vs_*",";","")
//	// in case the rotation has already been done
//	wList = RemoveFromList(WaveList("Vs_*_r",";",""),wList)
//	Variable nWaves = ItemsInList(wList)
//	Variable nCol = 3 // x y and z
//	Make/O/N=(nWaves,nCol) Img_VsCentre
//	Make/O/N=(nWaves) Img_VsMinAxis,Img_VsMajAxis,Img_VsPerimeter,Img_VsArea,Img_VsMean2R
//	String currentDF = GetDataFolder(0)
//	String wName, plotName
//	
//	Variable i,j
//	
//	for (i = 0; i < nWaves; i += 1)
//		wName = StringFromList(i, wList)
//		Wave w0 = $wName
//		nCol = dimsize(w0,1)
//		for (j = 0; j < nCol; j += 1)
//			WaveStats/Q/M=1/RMD=[][j] w0
//			Img_VsCentre[i][j] = V_avg
//		endfor
//		plotName = currentDF + "_allVsPlot"
//		AppendToGraph/W=$plotName w0[][1] vs w0[][0]
//		// offset to origin
//		ModifyGraph/W=$plotName offset($wName)={-Img_VsCentre[i][0],-Img_VsCentre[i][1]}
//		// find eigenvectors and rotate vesicle coords (also offset them)
//		Wave w1 = FindEV(w0)
//		plotName = currentDF + "_allRotVsPlot"
//		AppendToGraph/W=$plotName w1[][1] vs w1[][0]
//		Img_VsMinAxis[i] = VesicleAxisLength(w1,1)
//		Img_VsMajAxis[i] = VesicleAxisLength(w1,0)
//		Img_VsPerimeter[i] = FindLengthOfXYCoords(w1)
//		MatrixOp/O/FREE w1c0 = col(w1,0)
//		MatrixOp/O/FREE w1c1 = col(w1,1)
//		Img_VsArea[i] = PolygonArea(w1c0,w1c1)
//		// use polar coords to get mean radius, convert to diam
//		Make/O/N=(dimsize(w1,0)-1)/FREE rW
//		rW[] = sqrt(w1[p][0]^2 + w1[p][1]^2)
//		Img_VsMean2R[i] = 2 * mean(rW)
//	endfor
//	if(numpnts(Img_VsArea) > 1)
//		MatrixOp/O/NTHR=0 Img_VsAspectRatio = Img_VsMinAxis / Img_VsMajAxis
//		MatrixOp/O/NTHR=0 Img_VsCircularity = (4 * pi * Img_VsArea) / (Img_VsPerimeter * Img_VsPerimeter)
//	elseif(numpnts(Img_VsArea) == 1)
//		MatrixOp/O/NTHR=0 Img_VsAspectRatio = Img_VsMinAxis / Img_VsMajAxis
//		MatrixOp/O/FREE/NTHR=0 tempMat = Img_VsPerimeter * Img_VsPerimeter
//		MatrixOp/O/NTHR=0 Img_VsCircularity = 4 * pi * tempMat
//	else
//		Print "No vesicles in", currentDF
//	endif
//End
//
/////	@param	m1	2D wave of xy coords
//Function/WAVE FindEV(m1)
//	Wave m1
//	MatrixOp/O/FREE xCoord = col(m1,0)
//	MatrixOp/O/FREE yCoord = col(m1,1)
//	
//	// translate to origin
//	Variable offX = mean(xCoord)
//	Variable offY = mean(yCoord)
//	xCoord[] -= offX
//	yCoord[] -= offY
//	// do PCA. Rotated points are in M_R
//	PCA/ALL/SEVC/SRMT/SCMT xCoord,yCoord
//	WAVE M_R
//	String mName = NameOfWave(m1) + "_r"
//	Duplicate/O M_R, $mName
//	Wave m2 = $mName
//	Return m2
//End
//
/////	@param	m1	2D wave of xy coords
/////	@param	colNo	column number to use for search
//STATIC Function VesicleAxisLength(m1,colNo)
//	Wave m1
//	Variable colNo
//	MatrixOp/O/FREE m1c0 = col(m1,0)
//	MatrixOp/O/FREE m1c1 = col(m1,1)
//	Variable V_Value,len
//	if(colNo == 1)
//		FindLevel/Q/EDGE=1/P m1c0, 0
//		len = abs(m1c1(V_LevelX))
//		FindLevel/Q/EDGE=2/P m1c0, 0
//		len += abs(m1c1(V_LevelX))
//	else
//		FindLevel/Q/EDGE=1/P m1c1, 0
//		len = abs(m1c0(V_LevelX))
//		FindLevel/Q/EDGE=2/P m1c1, 0
//		len += abs(m1c0(V_LevelX))
//	endif
//	return len
//End
//
/////	@param	m1	2D wave of xy coords
//STATIC Function FindLengthOfXYCoords(m1)
//	Wave m1
//	// make new 2D wave of xy coords
//	Duplicate/O/FREE m1,tempDist
//	// offset to zero
//	tempDist[][0] -= m1[0][0]
//	tempDist[][1] -= m1[0][1]
//	// Differentiate, backward difference
//	Differentiate/METH=2 tempDist
//	// find norm, cumlative distance
//	MatrixOp/O/FREE/NTHR=0 tempNorm = sqrt(sumRows(tempDist * tempDist))
//	tempNorm[0] = 0 // first point is garbage
//	// return the sum of distances
//	return sum(tempNorm)
//End
//
//Function PlotSamples(rVs)
//	Variable rVs // requested sample size to plot out
//	
//	SetDataFolder root:data:
//	DFREF dfr = GetDataFolderDFR()
//	Variable numDataFolders = CountObjectsDFR(dfr, 4)
//	String wList, wName, plotName
//	Variable nVs = rVs
//	Variable nWaves
//	
//	Variable i,j
//	
//	// check if we need to alter nVs
//	for(i = 0; i < numDataFolders; i += 1)
//		SetDataFolder $("root:data:cond_" + num2str(i) + ":")
//		wList = WaveList("vs_*_r",";","")
//		nWaves = ItemsInList(wList)
//		nVs = min(nVs,nWaves)
//	endfor
//	
//	if (nVs < rVs)
//		Print "Condition(s) had less than", rVs, "contours. Plotting a sample of", nVs, "contours."	
//	endif
//	
//	for(i = 0; i < numDataFolders; i += 1)
//		SetDataFolder $("root:data:cond_" + num2str(i) + ":")
//		plotName = "cond_" + num2str(i) + "_allRotVsSamplePlot"
//		KillWindow/Z plotName
//		Display/N=$plotName/HIDE=1
//		// get list of all rotated contours (they are precentred so easier to plot)
//		wList = WaveList("vs_*_r",";","")
//		nWaves = ItemsInList(wList)
//		Make/O/N=(nWaves)/FREE rowW = p
//		StatsSample/N=(nVs) rowW
//		WAVE/Z W_Sampled
//		
//		for(j = 0; j < nVs; j += 1)
//			wName = StringFromList(W_Sampled[j],wList)
//			AppendToGraph/W=$plotName $wName[][1] vs $wName[][0]
//		endfor
//		KillWaves W_Sampled
//	endfor
//	SetDataFolder root:
//End
//
//Function FormatPlots()
//	String allGraphs = WinList("cond*",";","WIN:1")
//	Variable nWindows = ItemsInList(allGraphs)
//	String PlotName, boxString
//	Wave/Z/T labelWave = root:labelWave
//	Wave/Z colorWave = root:colorWave
//	Variable condRow
//	
//	Variable i
//	
//	for(i = 0; i < nWindows; i += 1)
//		plotName = StringFromList(i,allGraphs)
//		// 10% opacity
//		ModifyGraph/W=$plotName rgb=(0,0,0,6554)
//		SetAxis/W=$plotName left -300,300
//		SetAxis/W=$plotName bottom -300,300
//		ModifyGraph/W=$plotName width={Aspect,1}
//		ModifyGraph/W=$plotName mirror=1
//		condRow = str2num(plotName[5]) // a hack but will work
//		// add text colour for condition
//		boxString = "\\K(" + num2str(colorWave[condRow][0]) + "," + num2str(colorWave[condRow][1]) + "," + num2str(colorWave[condRow][2])
//		boxString += ")" + labelWave[condRow]
//		TextBox/W=$plotName/C/N=text0/F=0/A=RB/X=0.00/Y=0.00 boxString
//	endfor
//End
//
//
//STATIC Function CollectAllMeasurements()
//	SetDataFolder root:data:	// relies on earlier load
//	DFREF dfr = GetDataFolderDFR()
//	Variable numDataFolders = CountObjectsDFR(dfr, 4)
//	String wList = ""
//	
//	Variable i,j
//	// assemble a string of semi-colon separated targets in the data folder
//	// we'll do this by manually assembling the string to ensure the cond_* folders are in numerical order
//	for(i = 0; i < numDataFolders; i += 1)
//		wList += "root:data:cond_" + num2str(i) + ":thisWave;"
//	endfor
//	
//	SetDataFolder root:
//	// we need to concatenate these waves into root in a way that lets us build violin or box plots
//	String targetWaveList = "Img_VsArea;Img_VsAspectRatio;Img_VsCircularity;Img_VsMajAxis;Img_VsMinAxis;Img_VsMean2R;Img_VsPerimeter;"
//	Variable nTargets = ItemsInList(targetWaveList)
//	String fullName, targetName, tList
//	// we need to find the condition will the most contours
//	tList = ReplaceString("thisWave",wList,"Img_VsArea")
//	Variable mostContours = 0
//	
//	for(i = 0; i < numDataFolders; i +=1)
//		fullName = StringFromList(i, tList)
//		Wave testW = $fullName
//		mostContours = max(mostContours,numpnts(testW))
//	endfor
//	
//	String matName
//	Variable nRows
//	
//	for(i = 0; i < nTargets; i += 1)
//		// we pick a target
//		targetName = StringFromList(i,targetWaveList)
//		// assemble a list of those targets from each condition
//		tList = ReplaceString("thisWave",wList,targetName)
//		
//		for(j = 0; j < numDataFolders; j +=1)
//			// data folder path to the target for each condition
//			fullName = StringFromList(j, tList)
//			Wave targetWave = $fullName
//			nRows = numpnts(targetWave)
//			matName = "vb_" + num2str(j) + "_" + targetName
//			Make/O/N=(mostContours,numDataFolders) $matName = NaN
//			Wave matA = $matName
//			matA[0,nRows - 1][j] = targetWave[p]
//		endfor
//	endfor
//
//	String plotName
//	
//	for(i = 0; i < nTargets; i += 1)
//		plotName =  "p_" + StringFromList(i,targetWaveList)
//		KillWindow/Z $plotName
//		Display/N=$plotName
//		for(j = 0; j < numDataFolders; j += 1)
//			matName = "vb_" + num2str(j) + "_" + StringFromList(i,targetWaveList)
//			Wave matA = $matName
//			BuildBoxOrViolinPlot(matA,plotName,j)
//		endfor
//		SetAxis/A/N=1/E=1/W=$plotName left
//		ModifyGraph/W=$plotName toMode=-1
//		ModifyGraph/W=$plotName margin(left)=40
//	endfor
//	
//	// Label y-axes
//	Label/W=p_Img_VsArea left "Area (nm\S2\M)"
//	Label/W=p_Img_VsminAxis left "Minor axis (nm)"
//	Label/W=p_Img_VsmajAxis left "Major axis (nm)"
//	Label/W=p_Img_VsAspectRatio left "Aspect Ratio"
//	Label/W=p_Img_VsCircularity left "Circularity"
//	Label/W=p_Img_VsMean2R left "Diameter (nm)"
//	Label/W=p_Img_VsPerimeter left "Perimeter (nm)"
//End
//
//// This function will make a "multicolumn" boxplot or violinplot (Igor >8 only) 
/////	@param	matA	matrix of points to be appended
/////	@param	plotName	string to tell igor which graph window to work on
/////	@param	ii	variable to indicate which condition (for coloring)
//STATIC Function BuildBoxOrViolinPlot(matA,plotName,ii)
//	WAVE matA
//	String plotName
//	Variable ii
//	
//	String wName = NameOfWave(matA)
//	Wave/T/Z labelWave = root:labelWave
//	Wave/Z colorWave = root:colorWave
//	//  This works because all matrices passed to this function have the same dimensions
//	Variable nCells = DimSize(matA,0)
//	if(nCells < 100)
//		AppendBoxPlot/W=$plotName matA vs labelWave
//		ModifyBoxPlot/W=$plotName trace=$wName,markers={19,-1,19},markerSizes={2,2,2}
//		ModifyBoxPlot/W=$plotName trace=$wName,whiskerMethod=4
//	else
//		AppendViolinPlot/W=$plotName matA vs labelWave
//		ModifyViolinPlot/W=$plotName trace=$wName,ShowMean,MeanMarker=19,CloseOutline
//		ModifyViolinPlot/W=$plotName trace=$wName,DataMarker=19
//	endif
//	Variable alphaLevel = DecideOpacity(nCells)
//	ModifyGraph/W=$plotName rgb($wName)=(colorWave[ii][0],colorWave[ii][1],colorWave[ii][2],alphaLevel)
//End
//
////////////////////////////////////////////////////////////////////////
// Utility functions
////////////////////////////////////////////////////////////////////////
Function CleanSlate()
	String fullList = WinList("*", ";","WIN:71")
	Variable allItems = ItemsInList(fullList)
	String name
	Variable i
 
	for(i = 0; i < allItems; i += 1)
		name = StringFromList(i, fullList)
		KillWindow/Z $name		
	endfor
	
	KillDataFolder/Z root:data:
		
	// Kill waves in root
	KillWaves/A/Z
	// Look for data folders and kill them
	DFREF dfr = GetDataFolderDFR()
	allItems = CountObjectsDFR(dfr, 4)
	for(i = 0; i < allItems; i += 1)
		name = GetIndexedObjNameDFR(dfr, 4, i)
		KillDataFolder $name		
	endfor
End

STATIC Function KillTheseWaves(wList)
	String wList
	Variable nWaves = ItemsInList(wList)
	String wName
	
	Variable i
	
	for(i = 0; i < nWaves; i += 1)
		wName = StringFromList(i, wList)
		Wave w0 = $wName
		KillWaves/Z w0
	endfor
End
// Colours are taken from Paul Tol SRON stylesheet
// Colours updated. Brighter palette for up to 6 colours, then palette of 12 for > 6
// Define colours
StrConstant SRON_1 = "0x4477aa;"
StrConstant SRON_2 = "0x4477aa;0xee6677;"
StrConstant SRON_3 = "0x4477aa;0xccbb44;0xee6677;"
StrConstant SRON_4 = "0x4477aa;0x228833;0xccbb44;0xee6677;"
StrConstant SRON_5 = "0x4477aa;0x66ccee;0x228833;0xccbb44;0xee6677;"
StrConstant SRON_6 = "0x4477aa;0x66ccee;0x228833;0xccbb44;0xee6677;0xaa3377;"
StrConstant SRON_7 = "0x332288;0x88ccee;0x44aa99;0x117733;0xddcc77;0xcc6677;0xaa4499;"
StrConstant SRON_8 = "0x332288;0x88ccee;0x44aa99;0x117733;0x999933;0xddcc77;0xcc6677;0xaa4499;"
StrConstant SRON_9 = "0x332288;0x88ccee;0x44aa99;0x117733;0x999933;0xddcc77;0xcc6677;0x882255;0xaa4499;"
StrConstant SRON_10 = "0x332288;0x88ccee;0x44aa99;0x117733;0x999933;0xddcc77;0x661100;0xcc6677;0x882255;0xaa4499;"
StrConstant SRON_11 = "0x332288;0x6699cc;0x88ccee;0x44aa99;0x117733;0x999933;0xddcc77;0x661100;0xcc6677;0x882255;0xaa4499;"
StrConstant SRON_12 = "0x332288;0x6699cc;0x88ccee;0x44aa99;0x117733;0x999933;0xddcc77;0x661100;0xcc6677;0xaa4466;0x882255;0xaa4499;"

/// @param hex		variable in hexadecimal
Function hexcolor_red(hex)
	Variable hex
	return byte_value(hex, 2) * 2^8
End

/// @param hex		variable in hexadecimal
Function hexcolor_green(hex)
	Variable hex
	return byte_value(hex, 1) * 2^8
End

/// @param hex		variable in hexadecimal
Function hexcolor_blue(hex)
	Variable hex
	return byte_value(hex, 0) * 2^8
End

/// @param data	variable in hexadecimal
/// @param byte	variable to determine R, G or B value
STATIC Function byte_value(data, byte)
	Variable data
	Variable byte
	return (data & (0xFF * (2^(8*byte)))) / (2^(8*byte))
End

/// @param	nRow	variable for number of conditions
/// @param	wName	string to name the resulting colorwave
/// @param	[alpha]	optional variable for alpha column, n.b. 16-bit integer
Function MakeColorWave(nRow, wName, [alpha])
	Variable nRow
	String wName
	Variable alpha
	
	// Pick colours from SRON palettes
	String pal
	if(nRow == 1)
		pal = SRON_1
	elseif(nRow == 2)
		pal = SRON_2
	elseif(nRow == 3)
		pal = SRON_3
	elseif(nRow == 4)
		pal = SRON_4
	elseif(nRow == 5)
		pal = SRON_5
	elseif(nRow == 6)
		pal = SRON_6
	elseif(nRow == 7)
		pal = SRON_7
	elseif(nRow == 8)
		pal = SRON_8
	elseif(nRow == 9)
		pal = SRON_9
	elseif(nRow == 10)
		pal = SRON_10
	elseif(nRow == 11)
		pal = SRON_11
	else
		pal = SRON_12
	endif
	
	Variable color
	String colorWaveFullName = "root:" + wName
	if(ParamisDefault(alpha) == 1)
		Make/O/N=(nRow,3) $colorWaveFullName
		WAVE w = $colorWaveFullName
	else
		Make/O/N=(nRow,4) $colorWaveFullName
		WAVE w = $colorWaveFullName
		w[][3] = alpha
	endif
	
	Variable i
	
	for(i = 0; i < nRow; i += 1)
		// specify colours
		color = str2num(StringFromList(mod(i, 12),pal))
		w[i][0] = hexcolor_red(color)
		w[i][1] = hexcolor_green(color)
		w[i][2] = hexcolor_blue(color)
	endfor
End

STATIC Function MakeTheLayouts(prefix,nRow,nCol,[iter, filtVar, rev, alphaSort, saveIt, orient])
	String prefix
	Variable nRow, nCol
	Variable iter	// this is if we are doing multiple iterations of the same layout
	Variable filtVar // this is the object we want to filter for
	Variable rev // optional - reverse plot order
	Variable alphaSort // optional - do alphanumeric sort
	Variable saveIt
	Variable orient //optional 1 = landscape, 0 or default is portrait
	if(ParamIsDefault(filtVar) == 0)
		String filtStr = prefix + "_*_" + num2str(filtVar) + "_*"	// this is if we want to filter for this string from the prefix
	endif
	
	String layoutName = "all"+prefix+"Layout"
	DoWindow/K $layoutName
	NewLayout/N=$layoutName
	String allList = WinList(prefix+"*",";","WIN:1") // edited this line from previous version
	String modList = allList
	Variable nWindows = ItemsInList(allList)
	String plotName
	
	Variable i
	
	if(ParamIsDefault(filtVar) == 0)
		modList = "" // reinitialise
		for(i = 0; i < nWindows; i += 1)
			plotName = StringFromList(i,allList)
			if(stringmatch(plotName,filtStr) == 1)
				modList += plotName + ";"
			endif
		endfor
	endif
	
	if(ParamIsDefault(alphaSort) == 0)
		if(alphaSort == 1)
			modList = SortList(modList)
		endif
	endif
	
	nWindows = ItemsInList(modList)
	Variable PlotsPerPage = nRow * nCol
	String exString = "Tile/A=(" + num2str(ceil(PlotsPerPage/nCol)) + ","+num2str(nCol)+")"
	
	Variable pgNum=1
	
	for(i = 0; i < nWindows; i += 1)
		if(ParamIsDefault(rev) == 0)
			if(rev == 1)
				plotName = StringFromList(nWindows - 1 - i,modList)
			else
				plotName = StringFromList(i,modList)
			endif
		else
			plotName = StringFromList(i,modList)
		endif
		AppendLayoutObject/W=$layoutName/PAGE=(pgnum) graph $plotName
		if(mod((i + 1),PlotsPerPage) == 0 || i == (nWindows -1)) // if page is full or it's the last plot
			if(ParamIsDefault(orient) == 0)
				if(orient == 1)
					LayoutPageAction size(-1)=(842,595), margins(-1)=(18, 18, 18, 18)
				endif
			else
				// default is for portrait
				LayoutPageAction/W=$layoutName size(-1)=(595, 842), margins(-1)=(18, 18, 18, 18)
			endif
			ModifyLayout/W=$layoutName units=0
			ModifyLayout/W=$layoutName frame=0,trans=1
			Execute /Q exString
			if (i != nWindows -1)
				LayoutPageAction/W=$layoutName appendpage
				pgNum += 1
				LayoutPageAction/W=$layoutName page=(pgNum)
			endif
		endif
	endfor
	
	String fileName
	// if anthing is passed here we save an iteration, otherwise usual name
	if(!ParamIsDefault(iter))
		fileName = layoutName + num2str(iter) + ".pdf"
	else
		fileName = layoutName + ".pdf"
	endif
	// if anthing is passed here we save the filtered version
	if(ParamIsDefault(filtVar) == 0)
		fileName = ReplaceString(".pdf",fileName, "_" + num2str(filtVar) + ".pdf")
	endif
	if(ParamIsDefault(saveIt) == 0)
		if(saveIt == 1)
			SavePICT/O/WIN=$layoutName/PGR=(1,-1)/E=-2/W=(0,0,0,0) as fileName
		endif
	else
		// default is to save
		SavePICT/O/WIN=$layoutName/PGR=(1,-1)/E=-2/W=(0,0,0,0) as fileName
	endif
End

STATIC Function TidyAndSave(prefix)
	String prefix
	String layoutName = "all"+prefix+"Layout"
	// go to first page
	LayoutPageAction/W=$layoutName page=(1)
	// build the key
	WAVE/Z/T labelWave = root:labelWave
	WAVE/Z colorWave = root:colorWave
	Variable cond = numpnts(labelWave)
	String boxString = ""
	
	Variable i
	
	for(i = 0; i < cond; i += 1)
		// add text colour for condition
		boxString += "\\K(" + num2str(colorWave[i][0]) + "," + num2str(colorWave[i][1]) + "," + num2str(colorWave[i][2])
		boxString += ")" + labelWave[i]
		if (i < cond - 1)
			boxString += "\r"
		endif
	endfor
	TextBox/W=$layoutName/C/N=text0/F=0/A=RB/X=5.00/Y=5.00 boxString
	String fileName = layoutName + ".pdf"
	SavePICT/O/WIN=$layoutName/PGR=(1,-1)/E=-2/W=(0,0,0,0) as fileName
End

STATIC Function DecideOpacity(nTrace)
	Variable nTrace
	Variable alpha
	if(nTrace < 10)
		alpha = 1
	elseif(nTrace < 50)
		alpha = 0.5
	elseif(nTrace < 100)
		alpha = 0.3
	else
		alpha = 0.2
	endif
	alpha = round(65535 * alpha)
	return alpha
End

// for axis scaling
///	@param	value	this is the input value that requires rounding up
///	@param	roundto	round to the nearest...
STATIC Function RoundFunction(value,roundTo)
	Variable value, roundTo
	
	value /= roundTo
	Variable newVal = ceil(value)
	newVal *= roundTo
	return newVal
End

//////////////////////////////////////////////////////////////////////////
//// Panel functions
//////////////////////////////////////////////////////////////////////////
//STATIC Function Unblinder2User(condWave)
//	WAVE/T condWave
//	DoWindow/F SelectPanel
//	if (V_Flag != 0)
//		return 0
//	endif
//	Make/O/N=(numpnts(condWave)) root:testSel
//	NewPanel/N=SelectPanel/K=1/W=(350,125,650,325) as "Put conditions in order"
//	ListBox list0,pos={1,2},size={298,157},proc=ListBoxProc_DragNDropLB
//	ListBox list0,listWave=root:condWave,selWave=root:testSel,mode= 1,selRow= 1
//	Button DoIt,pos={100,177},size={100,20},proc=ButtonProc,title="Do It"
//End
// 
//Function ListBoxProc_DragNDropLB(lba) : ListBoxControl
//	STRUCT WMListboxAction &lba
// 
//	Variable row = lba.row
//	Variable col = lba.col
//	WAVE/T/Z listWave = lba.listWave
//	WAVE/Z selWave = lba.selWave
// 
//	switch( lba.eventCode )
//		case -1: // control being killed
//			break
//		case 1: // mouse down
//			Variable/G V_MouseDownRow = row
//			break
//		case 2: // mouse up
//			//quit if V_MouseDownRow isn't a string or numeric variable
//			//otherwise create a pointer to it.
//			if( exists( "V_MouseDownRow" ) == 2 )
//				NVAR V_MouseDownRow
//			else
//				break
//			endif
//			if(row != V_MouseDownRow)						// dragged?
//				String item = listWave[V_MouseDownRow]
//				DeletePoints V_MouseDownRow, 1, listWave	// do swap
//				InsertPoints row, 1, listWave
//				listWave[row] = item
//			endif
//			KillVariables V_MouseDownRow	// cleanup variable
//			break
//		case 3: // double click
//			break
//		case 4: // cell selection
//		case 5: // cell selection plus shift key
//			break
//		case 6: // begin edit
//			break
//		case 7: // finish edit
//			break
//	endswitch
// 	//print "code 1: ", lba.eventcode
//	return 1
//End
//
//// Define button
//Function ButtonProc(ctrlName) : ButtonControl
//	String ctrlName
//		strswitch(ctrlName) 
//			case "DoIt" :
//				// run the rest of the program
//				KillWindow/Z SelectPanel
//				OptionalAliases(0)
//				//break
//		endswitch
//End
//
//// This is the Alias check to make a label wave
//Function OptionalAliases(loadVar)
//	Variable loadVar // 0 is not loaded (need to load), 1 is already loaded
//	Make/O/N=(1) root:loadChoice = loadVar
//	Wave/T/Z condWave = root:condWave
//	Variable cond = numpnts(condWave)
//	MakeColorWave(cond, "colorWave")
//	Wave/Z colorWave = root:colorWave
//	// duplicate the condWave to make a tentative labelWave
//	Duplicate/O condWave, root:labelWave
//	Wave/T/Z labelWave = root:labelWave
//	
//	DoWindow/K AliasCheck
//	NewPanel/N=AliasCheck/K=1/W=(40, 40, 460, 150+30*cond)
//	// labelling of columns
//	DrawText/W=AliasCheck 10,30,"Name"
//	DrawText/W=AliasCheck 160,30,"Alias (a nice name for the plot labels)"
//	DrawText/W=AliasCheck 10,100+30*cond,"Liposome Size Analysis"
//	// do it button
//	Button DoIt,pos={280,70+30*cond},size={100,20},proc=DoItButtonProc,title="Do It"
//	// insert rows
//	String buttonName1a,buttonName1b,buttonName2a,buttonName2b,boxName0,boxName1,boxName2
//	Variable i
//	
//	for(i = 0; i < cond; i += 1)
//		boxName0 = "box0_" + num2str(i)
//		boxName1 = "box1_" + num2str(i)
//		// row label
//		DrawText/W=AliasCheck 10,68+i*30,num2str(i+1)
//		// condition label
//		SetVariable $boxName0,pos={30,53+i*30},size={100,14},value= condWave[i], title=" "
//		// file or dir box
//		SetVariable $boxName1,pos={160,53+i*30},size={220,14},value= labelWave[i], title=" "
//		SetDrawEnv fillfgc=(colorWave[i][0],colorWave[i][1],colorWave[i][2])
//		DrawOval/W=AliasCheck 130,50+i*30,148,68+i*30
//	endfor
//End
//
//// define buttons
//Function DoItButtonProc(ctrlName) : ButtonControl
//	String ctrlName
// 	
// 	WAVE/Z loadChoice = root:loadChoice
// 	WAVE/T/Z labelWave = root:labelWave
//	Variable okvar = 0
//	
//	strswitch(ctrlName)	
//		case "DoIt" :
//			// check CondWave
//			okvar = WaveChecker(labelWave)
//			if (okvar == -1)
//				Print "Error: Not all conditions have a name."
//				break
//			endif
//			okvar = NameChecker(labelWave)
//			if (okvar == -1)
//				Print "Error: Two conditions have the same name."
//				break
//			else
//				KillWindow/Z AliasCheck
//				// now execute main program
//				KillWaves/Z testSel
//				if(loadChoice[0] == 0)
//					TheLoader()
//				elseif(loadChoice[0] == 1)
//					TheProcessor()
//				endif
//				KillWaves/Z loadChoice
//				return 0
//			endif
//	endswitch	
//End
//
//STATIC function WaveChecker(TextWaveToCheck)
//	Wave/T TextWaveToCheck
//	Variable nRows = numpnts(TextWaveToCheck)
//	Variable len
//	
//	Variable i
//	
//	for(i = 0; i < nRows; i += 1)
//		len = strlen(TextWaveToCheck[i])
//		if(len == 0)
//			return -1
//		elseif(numtype(len) == 2)
//			return -1
//		endif
//	endfor
//	return 1
//End
//
//STATIC function NameChecker(TextWaveToCheck)
//	Wave/T TextWaveToCheck
//	Variable nRows = numpnts(TextWaveToCheck)
//	Variable len
//	
//	Variable i,j
//	
//	for(i = 0; i < nRows; i += 1)
//		for(j = 0; j < nRows; j += 1)
//			if(j > i)
//				if(cmpstr(TextWaveToCheck[i], TextWaveToCheck[j], 0) == 0)
//					return -1
//				endif
//			endif
//		endfor
//	endfor
//	return 1
//End
//
//
//STATIC function CheckThatAllDataExists(condW, exptW, fileList)
//	Wave/T condW, exptW
//	String fileList
//	
//	Variable nCond = numpnts(condW)
//	Variable nExpt = numpnts(exptW)
//	Variable nFiles = ItemsInList(fileList)
//	String stringToMatch
//	Variable matches = 0
//	
//	Variable i, j, k
//	
//	for(i = 0; i < nCond; i += 1)
//		for(j = 0; j < nExpt; j += 1)
//			stringToMatch = condW[i] + "_" + exptW[j] + "_*"
//			for(k = 0; k < nFiles; k += 1)
//				if(StringMatch(StringFromList(k,fileList), stringToMatch) == 1)
//					matches += 1
//					break
//				endif
//			endfor
//		endfor
//	endfor
//	
//	if(matches == nCond * nExpt)
//		return 0
//	else
//		return -1
//	endif
//End
//
//STATIC Function DoStatsAndLabel(m0,plotName)
//	Wave m0
//	String plotName
//	
//	String wName = NameOfWave(m0)
//	Variable groups = DimSize(m0,1)
//	Variable reps = DimSize(m0,0)
//	if(reps < 3)
//		Print "Less than three repeats, so no stats added to superplot"
//		return -1
//	endif
//	String pStr, boxName, lookup
//	Variable pVal, i
//	if(groups == 2)
//		Make/O/N=(reps)/FREE w0,w1
//		w0[] = m0[p][0]
//		w1[] = m0[p][1]
//		KillWaves/Z m0
//		StatsTTest/Q w0,w1
//		Wave/Z W_StatsTTest
//		pVal = W_StatsTTest[%P]
//		pStr = FormatPValue(pVal)
//		TextBox/C/N=text0/W=$plotName/F=0/A=MT/X=0.00/Y=0.00 "p = " + pStr
//	elseif(groups > 2)
//		SplitWave m0
//		StatsDunnettTest/Q/WSTR=S_WaveNames
//		WAVE/Z M_DunnettTestResults
//		for(i = 1; i < groups; i += 1)
//			boxName = "text" + num2str(i)
//			lookup = "0_vs_" + num2str(i)
//			pStr = FormatPValue(M_DunnettTestResults[%$(lookup)][%P])
//			TextBox/C/N=$boxName/W=$plotName/F=0/A=MT/X=(((i - (groups/2 - 0.5))/(groups / 2))/2 * 100)/Y=0.00 pStr
//		endfor
//		KillTheseWaves(S_WaveNames)
//	else
//		return -1
//	endif
//end
//
//STATIC Function/S FormatPValue(pValVar)
//	Variable pValVar
//	
//	String pVal = ""
//	String preStr,postStr
//	
//	if(pValVar > 0.05)
//		sprintf pVal, "%*.*f", 2,2, pValVar
//	else
//		sprintf pVal, "%*.*e", 1,1, pValVar
//	endif
//	if(pValVar > 0.99)
//		// replace any p ~= 1 with p > 0.99
//		pVal = "> 0.99"
//	elseif(pValVar == 0)
//		// replace any p = 0 with p < 1e-24
//		pVal = "< 1e-24"
//	endif
//	if(StringMatch(pVal,"*e*") == 1)
//		preStr = pVal[0,2]
//		if(StringMatch(pVal[5],"0") == 1)
//			postStr = pVal[6]
//		else
//			postStr = pVal[5,6]
//		endif
//		pVal = preStr + " x 10\S\u2212" + postStr
//	endif
//	return pVal
//End