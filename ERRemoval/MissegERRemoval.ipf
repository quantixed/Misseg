#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma DefaultTab={3,20,4}		// Set default tab width in Igor Pro 9 and later
#include <Waves Average>

// text files are outputs from ER_Removal.ijm
// this procedure will (currently) read in the files and do the calculations
// needs further work to take all the information out

////////////////////////////////////////////////////////////////////////
// Menu items
////////////////////////////////////////////////////////////////////////
Menu "Misseg"
	Submenu "ER Removal"
		"Analyse Removal...", /Q, ERAnalysis(0)
		"Analyse Removal Multilevel...", /Q, ERAnalysis(1)
	End
End

////////////////////////////////////////////////////////////////////////
// Master functions and wrappers
////////////////////////////////////////////////////////////////////////
Function ERAnalysis(opt)
	Variable opt
	CleanSlate()
	ERPreLoader(opt)
	GraphAll()
	MakeTheLayouts("p_", 6, 2, rev = 1, saveIt = 0)
End

////////////////////////////////////////////////////////////////////////
// Main functions
////////////////////////////////////////////////////////////////////////
Function ERPreLoader(opt)
	Variable opt
	
	NewPath/O/Q/M="Please find folder with measurement files" expDiskFolder1
	if (V_flag != 0)
		DoAlert 0, "Disk folder error"
		return -1
	endif
	String FileList1
	
	if (opt == 0)
		// this option assumes all files are in the top folder (expDiskFolder1)
		FileList1 = IndexedFile(expDiskFolder1,-1,".txt")
		if(ItemsInList(FileList1) == 0)
			DoAlert 0, "No txt files found"
			return -1
		endif
		Wave/T fileName1Wave = ListToTextWave(FileList1,";")
		MoveWave fileName1Wave, root:fileName1Wave
		LoadMeasurements(fileName1Wave, "expDiskFolder1", "data", 9, 2, 3, 0.0649645943)
	elseif (opt == 1)
		// this option assumes all files are in subfolders in (expDiskFolder1)
		// it assumes that there are one or more condition folders each containing exp folders
		// files are in experiment folders
		PathInfo /S expDiskFolder1
		String topDiskFolderName = S_path
		
		String dirList = IndexedDir(expDiskFolder1,-1,0)
		Variable nDirs = ItemsInList(dirList)
		
		Variable dirLoop, subDirLoop, fileLoop, nSubDirs, nFiles
		String condDiskFolderName, condDataFolderName
		String expDiskFolderName
		String subDirList, tmplist
		
		// need to make root:data: first
		NewDataFolder/O root:data
		
		for (DirLoop = 0; DirLoop < nDirs; DirLoop += 1)
			fileList1 = ""
			condDiskFolderName = StringFromList(dirLoop, dirList) // something like "wt"
			// now we need to get the experiment folders within the condition data folder
			condDiskFolderName = topDiskFolderName + condDiskFolderName + ":"
			NewPath/O/Q condDiskFolder, condDiskFolderName
			subDirList = IndexedDir(condDiskFolder,-1,0)
			nSubDirs = ItemsInList(subDirList)
			// find txt files in each subdir
			for (subDirLoop = 0; subDirLoop < nSubDirs; subDirLoop += 1)
				expDiskFolderName = StringFromList(subDirLoop, subDirList)
				expDiskFolderName = condDiskFolderName + expDiskFolderName + ":"
				NewPath/O/Q expDiskFolder, expDiskFolderName
				tmpList = IndexedFile(expDiskFolder,-1,".txt")
				tmpList = ListMatch(tmpList,"*_results.txt")
				nFiles = ItemsInList(tmpList)
				for (fileLoop = 0; fileLoop < nFiles; fileLoop += 1)
					fileList1 += expDiskFolderName + StringFromList(fileLoop, tmpList) + ";"
				endfor
			endfor
			Wave/T fileName1Wave = ListToTextWave(FileList1,";")
			condDataFolderName = "data:" + StringFromList(dirLoop, dirList)
			LoadMeasurements(fileName1Wave, "expDiskFolder1", condDataFolderName, 9, 2, 3, 0.0649645943)
			MoveWave fileName1Wave, $("root:fileNameWave_" + StringFromList(dirLoop, dirList))
		endfor
	endif
End

Function LoadMeasurements(fWaveToLoad, pathNameString, folderString, zMax, zStep, tStep, pxSize)
	Wave/T/Z fWaveToLoad
	String pathNameString,folderString
	Variable zMax, zStep, tStep, pxSize
	
	Variable vxSize = pxSize^2 * zStep
	
	Variable nFiles = numpnts(fWaveToLoad)
	String thisFile
	
	Variable i,j,k
	
	// set up folder
	NewDataFolder/O/S $("root:" + folderString)
	String wList = "XWave;AreaW;MeanW;StdDev;MinW;MaxW;IntDen;RawIntDen;Slice;Frame;"
	String hList = "row;theArea;theMean;StdDev;Min;Max;IntDen;RawIntDen;Slice;Frame;"
 	Variable nCol = ItemsInList(wList)
 	Variable maxT, maxSlice, matRow
 	
	for (i = 0; i < nFiles; i += 1)
		ThisFile = fWaveToLoad[i]
		LoadWave/A/J/D/O/Q/K=1/L={0,1,0,0,0}/P=$pathNameString/W ThisFile
		Concatenate/O/NP=1/KILL wList, $("mat_" + num2str(i))
		Wave matA = $("mat_" + num2str(i))
		// apply headers as dimlabels
		for (j = 0; j < nCol; j += 1)
			SetDimLabel 1,j,$(StringFromList(j,hList)),matA	
		endfor
		// now find the maximum frame number
		MatrixOp/O/FREE theTime = col(matA,9)
		maxT = WaveMax(theTime)
		MatrixOp/O/FREE theSlice = col(matA,8)
		maxSlice = WaveMax(theSlice)
		matRow = DimSize(matA,0)
		// make some waves to hold volume data for this cell
		Make/O/N=(maxT)/D $("cellVol_" + num2str(i)), $("ERVol_" + num2str(i))
		Wave cellW = $("cellVol_" + num2str(i))
		Wave erW = $("ERVol_" + num2str(i))

		// find biggest area for each slice at each timepoint
		Make/O/N=(matRow)/FREE biggest=0
		for (j = 0; j < maxT; j += 1)
			for (k = 0; k < maxSlice; k += 1)
				Make/O/N=(matRow)/FREE index
				MatrixOp/O/FREE thisArea = col(matA,1)
				index[] = (theTime[p] == j + 1 && theSlice[p] == k + 1) ? p : NaN
				thisArea[] = (theTime[p] == j + 1 && theSlice[p] == k + 1) ? thisArea[p] : NaN
				Wavetransform zapNaNs index
				Wavetransform zapNaNs thisArea
				if(numpnts(index) > 0)
					Sort/R thisArea, thisArea,index
					biggest[index[0]] = 1
				endif
			endfor
		endfor
		// we now have a wave called biggest where biggest for each Z/T combo is 1, all else is 1
		for (j = 0; j < maxT; j += 1)
			Make/O/N=(matRow)/FREE index
			index[] = (theTime[p] == j + 1 && biggest == 1) ? p : NaN
			Wavetransform zapNaNs index
			if(numpnts(index) == 0)
				cellW[j] = 0
				erW[j] = 0
			else
				// RawIntDen is the sum of pixels * 255 of ER in each slice at T
				Make/O/D/N=(numpnts(index))/FREE RawIntDen = matA[index[p]][7]
				erW[j] = (sum(RawIntDen) / 255) * vxSize
				Make/O/D/N=(numpnts(index))/FREE theArea = matA[index[p]][1]
				cellW[j] = (sum(theArea)) * zStep
			endif
		endfor
		SetScale/P x 0,tStep,"", erW
		SetScale/P x 0,tStep,"", cellW
		MatrixOp/O $("volRatio_" + num2str(i)) = erW / cellW
		SetScale/P x 0,tStep,"", $("volRatio_" + num2str(i))
	endfor
	SetDataFolder root:
End

Function GraphAll()
	// this is a wrapper
	GraphThem("p_cell_","cellVol_*","Cell Volume")
	GraphThem("p_er_","ERVol_*","ER Volume")
	GraphThem("p_ratio_","volRatio_*","ER Ratio")
	NormaliseThem()
	GraphThem("p_ncell_","cellVol_*_n","Cell Volume (normalized)")
	GraphThem("p_ner_","ERVol_*_n","ER Volume (normalized)")
	GraphThem("p_nratio_","volRatio_*_n","ER Ratio (normalized)")
End

Function GraphThem(plotName,wavesearch,yLab)
	String plotName, wavesearch, yLab
	SetDataFolder root:data
	
	DFREF dfr = GetDataFolderDFR()
	Variable allItems = CountObjectsDFR(dfr, 4)
	Wave/Z colorWave = root:colorWave
	if(!WaveExists(colorWave))
		MakeColorWave(allItems, "colorWave", alpha = round(65535 / 2))
		Wave/Z colorWave = root:colorWave
	endif
	
	String dfName, wList, thisPlotName, wName
	Variable nWaves
	Variable i,j
	
	for(i = 0; i < allItems; i += 1)
		dfName = GetIndexedObjNameDFR(dfr, 4, i)
		SetDataFolder $("root:data:" + dfName)
		thisPlotName = plotName + dfName
		wList = WaveList(wavesearch,";","")
		nWaves = ItemsInList(wList)
		Display/N=$thisPlotName/HIDE=1
		for(j = 0; j < nWaves; j += 1)
			wName = StringFromList(j,wList)
			AppendToGraph/W=$thisPlotName $wName
		endfor
		SetAxis/A/N=1/E=1/W=$thisPlotName left
		Label/W=$thisPlotName left yLab
		Label/W=$thisPlotName bottom "Time (min)"
		ModifyGraph/W=$thisPlotName rgb=(colorWave[i][0],colorWave[i][1],colorWave[i][2],colorWave[i][3])
		ModifyGraph/W=$thisPlotName margin(left)=54
		AddOrUpdateAverage(thisPlotName)
	endfor
	SetDataFolder root:
End

Function AddOrUpdateAverage(plotName)
	String plotName
	
	String tList = TraceNameList(plotName,";",1)
	Variable nTraces = ItemsInList(tList)
	String tName
	String wList = ""
	Variable AverageWavePresent = 0
	
	Variable i
	
	for(i = 0; i < nTraces; i += 1)
		tName = StringFromList(i, tList)
		if(stringmatch(tName,"W_Ave*") == 0)
			Wave w = TraceNameToWaveRef(plotName,tName)
			wList += GetWavesDataFolder(w,2) + ";"
		else
			AverageWavePresent = 1
		endif
	endfor
	
	SetDataFolder root:
	String avName = ReplaceString("p_",plotName,"W_Ave_")
	String errName = ReplaceString("Ave", avName, "Err")
	fWaveAverage(wList, "", 3, 1, avName, errName)
	if(AverageWavePresent == 0)
		AppendToGraph/W=$plotName $avName
		ErrorBars/W=$plotName $avName SHADE= {0,0,(0,0,0,0),(0,0,0,0)},wave=($ErrName,$ErrName)
		ModifyGraph/W=$plotName lsize($avName)=2,rgb($avName)=(0,0,0)
	endif
End

Function NormaliseThem()
	SetDataFolder root:data
	
	DFREF dfr = GetDataFolderDFR()
	Variable allItems = CountObjectsDFR(dfr, 4)
	
	String dfName, wList, wName, newWName
	Variable nWaves, normVal
	Variable i,j
	
	for(i = 0; i < allItems; i += 1)
		dfName = GetIndexedObjNameDFR(dfr, 4, i)
		SetDataFolder $("root:data:" + dfName)
		wList = WaveList("*",";","")
		wList = RemoveFromList(WaveList("mat*",";",""),wList)
		nWaves = ItemsInList(wList)
		for(j = 0; j < nWaves; j += 1)
			wName = StringFromList(j,wList)
			Wave w0 = $wName
			newWName = wName + "_n"
			Duplicate/O w0 $newWName
			Wave w1 = $newWName
			normVal = (w0[0] + w0[1]) / 2 // use this rather than mean due to uncertainty of scaling
			w1[] = w0[p] / normVal
		endfor
	endfor
	SetDataFolder root:
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

// specific to this pxp
// no error checking
Function RemoveCaseFromPlots(cond,cell)
	String cond
	Variable cell
	
	String plotName
	
	plotName = "p_cell_" + cond
	RemoveFromGraph/W=$plotName $("cellVol_" + num2str(cell))
	AddOrUpdateAverage(plotName)
	plotName = "p_er_" + cond
	RemoveFromGraph/W=$plotName $("ERVol_" + num2str(cell))
	AddOrUpdateAverage(plotName)
	plotName = "p_ratio_" + cond
	RemoveFromGraph/W=$plotName $("volRatio_" + num2str(cell))
	AddOrUpdateAverage(plotName)

	plotName = "p_ncell_" + cond
	RemoveFromGraph/W=$("p_ncell_" + cond) $("cellVol_" + num2str(cell) + "_n")
	AddOrUpdateAverage(plotName)
	plotName = "p_ner_" + cond
	RemoveFromGraph/W=$("p_ner_" + cond) $("ERVol_" + num2str(cell) + "_n")
	AddOrUpdateAverage(plotName)
	plotName = "p_nratio_" + cond
	RemoveFromGraph/W=$("p_nratio_" + cond) $("volRatio_" + num2str(cell) + "_n")
	AddOrUpdateAverage(plotName)
End