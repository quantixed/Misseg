#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// Data one the onset of three phases NEB, metaphase, anaphase is recorded (frame number)
// in 3 columns per sheet (each sheet is a different condition)
// Data should start at A1

////////////////////////////////////////////////////////////////////////
// Menu items
////////////////////////////////////////////////////////////////////////
Menu "Macros"
	"Mitotic Progression...", /Q, MitoticProgression()
End


////////////////////////////////////////////////////////////////////////
// Master functions and wrappers
////////////////////////////////////////////////////////////////////////
Function MitoticProgression()
	CleanSlate()
	LoadDataFromExcel()
	CalculateTimings(3)
	MakeHistograms()
End

////////////////////////////////////////////////////////////////////////
// Main functions
////////////////////////////////////////////////////////////////////////
Function LoadDataFromExcel()
	// each experimental condition needs to be a separate sheet
	String sheet,wList,wName
	Variable i,j
	
	XLLoadWave/J=1
	if(V_Flag)
		DoAlert 0, "The user pressed Cancel"
	endif
	Variable nSheets = ItemsInList(S_value)
	NewPath/O/Q path1, S_path
	String colList = "NEB;Metaphase;Anaphase;"
	String nameList
	Make/O/T/N=(nSheets) condWave
	Variable nExp,nCell,Length
	String newName
	
	for(i = 0; i < nSheets; i += 1)
		sheet = StringFromList(i,S_Value)
		condWave[i] = sheet
		nameList = PrefixStringList(sheet,colList)
		XLLoadWave/S=sheet/D/W=1/NAME=nameList/O/K=0/P=path1 S_fileName
		// now check lengths. Any NaNs at end of wave will be truncated by XLLoadWave. Concatenate will fail
		wList = WaveList(sheet + "_*",";","")
		nExp = ItemsInList(wList)
		nCell = 0
		for(j = 0; j < nExp; j += 1)
			wName = StringFromList(j,wList)
			Wave w0 = $wName
			nCell = max(numpnts(w0),nCell) // find max length
		endfor
		// correct the lengths
		for(j = 0; j < nExp; j += 1)
			wName = StringFromList(j,wList)
			Wave w0 = $wName
			length = numpnts(w0)
			if (length < nCell)
				InsertPoints length, (nCell - length), w0 // add NaNs if necessary
			endif
		endfor	
		// make matrix
		Concatenate/O/NP=1/KILL nameList, $(sheet + "_frames")
	endfor
	// Print "***\r The sheet", sheet, "was loaded from", S_path,"\r  ***"
End

///	@param	prefix	string for inserting before each item in list
///	@param	strList	the string list to be modified.
Function/S PrefixStringList(prefix,strList)
	String prefix,strList
	String newString = ""
	
	Variable i
	
	for(i = 0; i < ItemsInList(strList); i += 1)
		newString += prefix + "_" + StringFromList(i,strList) + ";"
	endfor
	
	return newString
End

///	@param	rate	variable, minutes per frame
Function CalculateTimings(rate)
	Variable rate
	WAVE/Z/T condWave
	
	Variable nCond = numpnts(condWave)
	Variable nCell
	String wName, newName
	Variable i
	
	for(i = 0; i < nCond; i += 1)
		wName = condWave[i] + "_frames"
		Wave w0 = $wName
		nCell = DimSize(w0,0)
		// NEB-to-metaphase
		newName = ReplaceString("_frames",wName,"_NM")
		Make/O/N=(nCell) $newName
		Wave w1 = $newName
		w1[] = (w0[p][1] - w0[p][0]) * rate
		w1[] = (w1[p] < 0) ? NaN : w1[p] // correct any negative values
		// metaphase-to-anaphase
		newName = ReplaceString("_frames",wName,"_MA")
		Make/O/N=(nCell) $newName
		Wave w1 = $newName
		w1[] = (w0[p][2] - w0[p][1]) * rate
		w1[] = (w1[p] < 0) ? NaN : w1[p]
		// NEB-to-anaphase
		newName = ReplaceString("_frames",wName,"_NA")
		Make/O/N=(nCell) $newName
		Wave w1 = $newName
		w1[] = (w0[p][2] - w0[p][0]) * rate
		w1[] = (w1[p] < 0) ? NaN : w1[p]
		// Make 1D wave of NEB onset (used for plotting durations vs onset)
		newName = ReplaceString("_frames",wName,"_onset")
		Make/O/N=(nCell) $newName
		Wave w1 = $newName
		w1[] = w0[p][0] * rate
	endfor
End

Function MakeHistograms()
	WAVE/Z/T condWave
	
	Variable nCond = numpnts(condWave)
	MakeColorWave(nCond,"colorWave",alpha = 32639)
	WAVE/Z colorWave
	
	String plotList = "NM;MA;NA;"
	String labelList = "NEB-Meta;Meta-Ana;NEB-Ana;"
	String plotName, wName, histName
	Variable nPlots = ItemsInList(plotList)
	DoWindow/K summaryLayout
	NewLayout/N=summaryLayout
	LayoutPageAction/W=summaryLayout size(-1)=(595, 842), margins(-1)=(18, 18, 18, 18)
	
	Variable maxLength = 0
	String legendStr = "", labelStr = ""
	Variable normVar
	
	Variable i,j
	
	for(i = 0; i < nPlots; i += 1)
		plotname = "p_" + StringFromList(i,plotList)
		KillWindow/Z $plotName
		Display/N=$plotName
		maxLength = 0
		legendStr = ""
		
		for(j = 0; j < nCond; j += 1)
			wName = condWave[j] + "_" + StringFromList(i,plotList)
			Wave w0 = $wName
			maxLength = max(WaveMax(w0),maxLength)
		endfor
		
		for(j = 0; j < nCond; j += 1)
			wName = condWave[j] + "_" + StringFromList(i,plotList)
			Wave w0 = $wName
			histName = wName + "_hist"
			Make/O/N=(ceil((maxLength) / 3) + 2) $histName
//			Histogram/CUM/P/B={0,3,ceil((maxLength + 1) / 3)} w0,$histName
			Histogram/CUM/B={0,3,ceil((maxLength) / 3) + 2} w0,$histName
			Wave w1 = $histName
//			normVar = WaveMax(w1) // this normalises to max of cumHist
			normVar = numpnts(w0)
			w1 /= normVar
			AppendToGraph/W=$plotName $histName
			ModifyGraph/W=$plotName rgb($histName)=(colorWave[j][0],colorWave[j][1],colorWave[j][2],colorWave[j][3])
			if(j > 0)
				legendStr += "\r"
			endif
			legendStr += "\\s(" + histName + ") " + condWave[j]
			// add to labelStr
			labelStr += "\r" + NameOfWave(w0) + " progression: " + num2str(w1[numpnts(w1) - 1]) + " for " + num2str(normVar) + "cells."
			FindLevel/Q w1, (w1[numpnts(w1) - 1] * 0.5)
			labelStr += "T-half: " + num2str(V_LevelX)
		endfor
		FormatHisto(plotName, legendStr, StringFromList(i,labelList))
		AppendLayoutObject/W=summaryLayout graph $plotName
	endfor
	
	// make a plot of both
	plotName = "p_NMNA"
	KillWindow/Z $plotName
	Display/N=$plotName
	for(i = 0; i < nCond; i += 1)
		wName = condWave[i] + "_NM"
		histName = wName + "_hist"
		AppendToGraph/W=$plotName $histName
		ModifyGraph/W=$plotName rgb($histName)=(colorWave[i][0],colorWave[i][1],colorWave[i][2],colorWave[i][3])
		wName = condWave[i] + "_NA"
		histName = wName + "_hist"
		AppendToGraph/W=$plotName $histName
		ModifyGraph/W=$plotName rgb($histName)=(colorWave[i][0],colorWave[i][1],colorWave[i][2],colorWave[i][3])
	endfor
	// lengendStr should still be correct
	FormatHisto(plotName, legendStr, "")
	AppendLayoutObject/W=summaryLayout graph $plotName
	
	// now add duration vs onset
	for(i = 0; i < nPlots; i += 1)
		plotname = "q_" + StringFromList(i,plotList)
		KillWindow/Z $plotName
		Display/N=$plotName
		
		for(j = 0; j < nCond; j += 1)
			wName = condWave[j] + "_" + StringFromList(i,plotList)
			AppendToGraph/W=$plotName $wName vs $(condWave[j] + "_onset")
			ModifyGraph/W=$plotName rgb($wName)=(colorWave[j][0],colorWave[j][1],colorWave[j][2],colorWave[j][3])
		endfor
		ModifyGraph/W=$plotName mode=3,marker=19,mrkThick=0
		SetAxis/W=$plotName/A/N=1/E=1 left
		SetAxis/W=$plotName/A/N=1/E=1 bottom
		Label/W=$plotName left "Duration (min)"
		Label/W=$plotName bottom "Onset (min)"
		TextBox/W=$plotName/C/N=text0/F=0/S=3/A=RT/X=0.00/Y=0.00 StringFromList(i,labelList)
		AppendLayoutObject/W=summaryLayout graph $plotName
	endfor
	
	
	// format layout
	ModifyLayout/W=summaryLayout units=0
	ModifyLayout/W=summaryLayout frame=0,trans=1
	Execute /Q "Tile/A=(6,3)"
	
	TextBox/W=summaryLayout/C/N=text0/F=0/A=LB/X=0.00/Y=0.00 labelStr
End

STATIC Function FormatHisto(plotName, legendStr, labelStr)
	String plotName, legendStr, labelStr
	
	SetAxis/A/N=1/W=$plotName left
	SetAxis/A/N=1/W=$plotName bottom
	Label/W=$plotName bottom "Time (min)"
	Label/W=$plotName left "Cumulative frequency"
	TextBox/W=$plotName/C/N=text0/F=0/S=3/A=RB/X=0.00/Y=0.00 legendStr
	TextBox/W=$plotName/C/N=text1/F=0/S=3/A=LT/X=0.00/Y=0.00 labelStr
	ModifyGraph/W=$plotName lsize=2
End

////////////////////////////////////////////////////////////////////////
// Utility functions
////////////////////////////////////////////////////////////////////////
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

/// @param	cond	variable for number of conditions
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

STATIC Function CleanSlate()
	String fullList = WinList("*", ";","WIN:7")
	Variable allItems = ItemsInList(fullList)
	String name
	Variable i
 
	for(i = 0; i < allItems; i += 1)
		name = StringFromList(i, fullList)
		KillWindow/Z $name		
	endfor
	
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