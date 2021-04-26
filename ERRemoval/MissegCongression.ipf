#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// CSVs from ROI Monager
// Naming of files is important

////////////////////////////////////////////////////////////////////////
// Menu items
////////////////////////////////////////////////////////////////////////
Menu "Macros"
	"Congression Analysis...", ROIAnalysis()
	"Save PXP", SaveIt()
	"Start Over", CleanSlate()
	End
End

////////////////////////////////////////////////////////////////////////
// Master functions and wrappers
////////////////////////////////////////////////////////////////////////
Function ROIAnalysis()
	PreLoader()
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

// after code is run, we can save to the directory containing chromosome/plate directories
// use name of that directory for pxp
// check first that the code has been run
Function SaveIt()
	WAVE/Z c1c2distance
	if(!WaveExists(c1c2distance))
		DoAlert 0, "No data to save. Run main program and try again"
		return 0
	endif
	PathInfo expDiskFolder1
	String parentDiskFolderPath = ParseFilePath(1,S_path,":",1,0)
	NewPath/O/Q parentDiskFolder, parentDiskFolderPath
	String pxpName = ParseFilePath(0,parentDiskFolderPath,":",1,0) + ".pxp"
	SaveExperiment/P=parentDiskFolder as pxpName
End

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
