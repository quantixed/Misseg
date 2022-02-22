#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3				// Use modern global access method and strict wave access
#pragma DefaultTab={3,20,4}		// Set default tab width in Igor Pro 9 and later
#include "PXPUtils"
#include <Waves Average>

// This version of the ipf works with ROIs generated from 
// a script called ChromosomeCongresionMagicWand.ijm
// Now, the outputs are nested in folders and there are 3 separate csvs (coords, orig, enlg)
// for each obj. 0 is the plate 1 and greater are the chromosomes.
// In a change to the previous version we need to deal with this change and remove the
// linkage sections of code (snakes etc.)

////////////////////////////////////////////////////////////////////////
// Menu items
////////////////////////////////////////////////////////////////////////

Menu "Macros"
	"Run Chromosome Congression", /Q, ChromConWrapper()
End

////////////////////////////////////////////////////////////////////////
// Master functions and wrappers
////////////////////////////////////////////////////////////////////////

Function ChromConWrapper()
	PXPUtils#CleanSlate()
	SetUpParamWave()
	PreLoader()
	PXPUtils#MakeTheLayouts("p",6,4, alphaSort = 1, saveIt = 0)
	FurtherAnalysisEngine()
End

Function PreLoader()	
	NewPath/O/Q/M="Please find top folder with output subfolders" expDiskFolder
	if (V_flag != 0)
		DoAlert 0, "Disk folder error"
		return -1
	endif
	
	PathInfo expDiskFolder
	String pathToEDF = S_path
	
	String dirList = IndexedDir(expDiskFolder,-1,0)
	Wave/T dirNameWave = ListToTextWave(dirList,";")
	MoveWave dirNameWave, root:dirNameWave
	Variable nDir = numpnts(dirNameWave)
	String fileList

	// set up data folder
	NewDataFolder/O/S root:data
	Variable i

	for (i = 0; i < nDir; i += 1)
		NewDataFolder/O/S $("root:data:cell" + num2str(i))
		NewPath/O/Q subDiskFolder, pathToEDF + dirNameWave[i]
		fileList = IndexedFile(subDiskFolder,-1,".csv")
		if(ItemsInList(fileList) == 0)
			continue
		endif
		Wave/T csvNameWave = ListToTextWave(FileList,";")
		MoveWave csvNameWave, fileNameWave
		LoadFiles(fileList)
		ChromConEngine()
	endfor
	
	SetDataFolder root:
	
	return 0
End

Function LoadFiles(fList)
	String fList
	
	fList = SortList(fList)
	String fName, labelStr, objStr, newFName
	
	String subList = ListMatch(fList,"*coords.csv")
	Variable nFiles = ItemsInList(subList)
	
	Variable i
	
	for(i = 0; i < nFiles; i += 1)
		fName = StringFromList(i, subList)
		objStr = fName[strsearch(fName,"_coords.csv",0) - 1]
		LoadWave/A/O/D/J/W/Q/P=subDiskFolder fName
		labelStr = "roi;obj;xcoord;ycoord;zpos;tpos;"
		Concatenate/O/NP=1/KILL labelStr, $("coords_" + objStr)
		Wave m0 = $("coords_" + objStr)
		LabelWaveDimensions(m0,1,labelStr)
		// convert data to real units and make centroids for plotting
		ProcessCoords(m0)
		// now load measurements - note that we now take maxW instead of rawintden
		labelStr = "roi;obj;areaW;meanW;maxW;ch;slice;frame;"
		newFName = ReplaceString("_coords.csv",fName,"_orig.csv")
		LoadWave/A/O/D/J/W/Q/P=subDiskFolder newFName
		Concatenate/O/KILL/NP=1 labelStr, $("measOrig_" + objStr)
		Wave m1 = $("measOrig_" + objStr)
		LabelWaveDimensions(m1,1,labelStr)
		
		newFName = ReplaceString("_coords.csv",fName,"_enlg.csv")
		LoadWave/A/O/D/J/W/Q/P=subDiskFolder newFName
		Concatenate/O/KILL/NP=1 labelStr, $("measEnlg_" + objStr)
		Wave m2 = $("measEnlg_" + objStr)
		LabelWaveDimensions(m2,1,labelStr)
	endfor
End


Function ChromConEngine()
	// at this point we have the coords (and the scaled coords) plus two sets of measurements
	// for all objects. Plate is 0 and chromosomes are 1 upwards.
	// Still working within a datafolder...
	CalculateChromosomePlateDistances()	
	// Plot out the distances
	MakeDistancesPlot()
	// retrieve and plot intensities
	ChromosomeFluorescenceIntensities()
End

Function FurtherAnalysisEngine()
	PlotVsDistance()
	AverageTheWaves("ch1")
	AverageTheWaves("ch2")
	MakeDistancePlotsFromLastChromosomeTraces("q_ch1_con_raw")
	MakeDistancePlotsFromLastChromosomeTraces("q_ch1_mis_raw")
	CompareDistanceMeans()
	SurvivalStyleDistancePlot()
	MakeAndPlotFits("q_ch1_con_raw")
	MakeAndPlotFits("q_ch1_mis_raw")
	MakeAndPlotFits("q_ch2_con_raw")
	MakeAndPlotFits("q_ch2_mis_raw")
	PXPUtils#MakeTheLayouts("q",6,4, alphaSort = 1, saveIt = 0)
End


////////////////////////////////////////////////////////////////////////
// Main functions
////////////////////////////////////////////////////////////////////////

STATIC Function SetUpParamWave()
	Make/O/N=(5) paramWave = {0.1819,0.1819,2,3,3}
	SetDimLabel 0, 0, xx, paramWave // microns
	SetDimLabel 0, 1, yy, paramWave // microns
	SetDimLabel 0, 2, zz, paramWave // microns
	SetDimLabel 0, 3, tt, paramWave // minutes
	SetDimLabel 0, 4, cc, paramWave // number of channels
End

Function ProcessCoords(matA)
	Wave matA // 2D matrix of data from csv load
	
	Wave/Z paramWave = root:paramWave
	
	String newName = ReplaceString("coords_",NameOfWave(matA),"coordsReal_")
	Duplicate/O/RMD=[][2,4] matA, $newName
	Wave w = $newName
	SetScale/P y 0,1,"", w
	w[][0,1] *= paramWave[%xx]
	w[][2] *= paramWave[%zz]
End

Function CalculateChromosomePlateDistances()
	Wave/Z paramWave = root:paramWave
	
	String wList = WaveList("coords_*",";","")
	wList = RemoveFromList("coords_0",wList)
	Variable nWaves = ItemsInList(wList)
	nWaves = ItemsInList(wList)
	Wave/Z coW0 = $(GetDataFolder(1) + "coords_0")
	Wave/Z coRW0 = $(GetDataFolder(1) + "coordsReal_0")
	if(!WaveExists(coW0))
		return -1
	endif
	
	// find t limits (in frames) for plate
	MatrixOp/O/FREE allT = col(coW0,5)
	FindDuplicates/RN=uT allT
	Sort uT, uT // just in case they are out of sequence
	Variable nT = numpnts(uT)
	// to store time point (in real units)
	Make/O/N=(nT) tReal_0
	tReal_0[] = (uT[p] - 1) * paramWave[%tt] // real t
	// we'll to store centroid for this object here (in real units)
	Make/O/N=(nT,3) $(GetDataFolder(1) + "centroid_0")
	Wave cW0 = $(GetDataFolder(1) + "centroid_0")
	String wName
	
	Variable i,j
	
	// find the centroid and time for plate
	for(i = 0; i < nT; i += 1)
		// extract the xyz real coords for this frame
		Duplicate/O/FREE coRW0, tempMat
		tempMat[][] = (coW0[p][%tPos] == uT[i]) ? coRW0[p][q] : NaN
		MatrixOp/O/FREE tempSMat = zapnans(tempMat)
		if(numpnts(tempSMat) > 0)
			// convert back to 2D wave
			Redimension/N=(numpnts(tempSMat)/DimSize(coRW0,1),DimSize(coRW0,1)) tempSMat
			// find centroid for object in this frame
			MatrixOp/O/FREE resultW = sumcols(tempSMat) / numrows(tempSMat)
			cW0[i][] = resultW[0][q]
		else
			cW0[i][] = NaN
		endif
	endfor
	
	for(i = 0; i < nWaves; i += 1)
		wName = StringFromList(i,wList)
		Wave coW1 = $wName
		Wave coRW1 = $(ReplaceString("coords_",wName,"coordsReal_"))
		// we'll to store centroid for this object here (in real units)
		Make/O/N=(nT,3) $(ReplaceString("coords_",wName,"centroid_"))
		Wave cW1 = $(ReplaceString("coords_",wName,"centroid_"))
	
		for(j = 0; j < nT; j += 1)
			// extract the xyz real coords for this frame
			Duplicate/O/FREE coRW1, tempMat
			tempMat[][] = (coW1[p][%tPos] == uT[j]) ? coRW1[p][q] : NaN
			MatrixOp/O/FREE tempSMat = zapnans(tempMat)
			if(numpnts(tempSMat) > 0)
				// convert back to 2D wave
				Redimension/N=(numpnts(tempSMat)/DimSize(coRW1,1),DimSize(coRW1,1)) tempSMat
				// find centroid for object in this frame
				MatrixOp/O/FREE resultW = sumcols(tempSMat) / numrows(tempSMat)
				cW1[j][] = resultW[0][q]
			else
				cW1[j][] = NaN
			endif
		endfor
		MatrixOp/O $(ReplaceString("coords_",wName,"distance_")) = sqrt(sumrows((cW1 - cW0) * (cW1 - cW0)))
		Wave dW = $(ReplaceString("coords_",wName,"distance_"))
		dW[] = (numtype(dW[p]) == 2) ? 0 : dW[p]
	endfor
End

Function MakeDistancesPlot()
	String wList = WaveList("distance_*",";","")
	Variable nWaves = ItemsInList(wList)
	WAVE/Z tReal_0
	if(!WaveExists(tReal_0))
		return -1
	endif
	
	String ii = ReplaceString("cell",GetDataFolder(0),"")
	String plotName = "p_" + ii + "_chrDist"
	KillWindow/Z $plotName
	Display/N=$plotName/HIDE=1
	
	Variable i
	
	for(i = 0; i < nWaves; i += 1)
		AppendToGraph/W=$plotName $(StringFromList(i,wList)) vs tReal_0
	endfor
	
	Label/W=$plotName left "Distance (µm)"
	Label/W=$plotName bottom "Time (min)"
//	SetAxis/W=$plotName/A/N=1/E=1 left
	SetAxis/W=$plotName left 0,10
	ModifyGraph/W=$plotName rgb=(0,0,0)
End


Function ChromosomeFluorescenceIntensities()
	Wave/Z paramWave = root:paramWave
	
	String wList = WaveList("coords_*",";","")
	Variable nWaves = ItemsInList(wList)
	WAVE/Z uT, tReal_0
	if(!WaveExists(tReal_0))
		return -1
	endif
	
	String ii = ReplaceString("cell",GetDataFolder(0),"")
	String plotName, wName, newName
	
	Variable i,j
	
	// graphs to display chromosome intensities
	for(i = 0; i < paramWave[%cc]; i += 1)
		plotName = "p_" + ii + "_chrCh" + num2str(i+1)
		KillWindow/Z $plotName
		Display/N=$plotName/HIDE=1
	endfor
	
	// get measurements and display
	for(i = 0; i < nWaves; i += 1)
		wName = StringFromList(i,wList)
		Wave measureBig = $(ReplaceString("coords_",wName,"measEnlg_"))
		// we need mean for ch2 and ch3 but max for ch1
		MatrixOp/O/FREE diffMeanW = col(measureBig,3)
		MatrixOp/O/FREE diffMaxW = col(measureBig,4)
		
		for(j = 0; j < paramWave[%cc]; j += 1)
			newName = (ReplaceString("coords_",wName,"diff_ch" + num2str(j + 1) + "_"))
			Make/O/N=(numpnts(uT)) $newName
			Wave w = $newName
			if(j + 1 == 1)
				w[] = GetMeasurement(diffMaxW,measureBig,j+1,uT[p])
			else
				w[] = GetMeasurement(diffMeanW,measureBig,j+1,uT[p])
			endif
			plotName = "p_" + ii + "_chrCh" + num2str(j+1)
			AppendToGraph/W=$plotName w vs tReal_0
			if(StringMatch(wName,"*_0"))
				ModifyGraph/W=$plotName lstyle($newName)=2
			endif
		endfor
	endfor
	
	Make/O/N=(3,3) colorWave = {{0,166,81},{237,28,36},{34,118,185}}
	colorWave *= 257 // convert to 16-bit
	MatrixTranspose colorWave
	
	for(i = 0; i < paramWave[%cc]; i += 1)
		plotName = "p_" + ii + "_chrCh" + num2str(i+1)
		Label/W=$plotName left "Intensity"
		Label/W=$plotName bottom "Time (min)"
		SetAxis/W=$plotName/A/N=1/E=3 left
		ModifyGraph/W=$plotName rgb=(colorWave[i][0],colorWave[i][1],colorWave[i][2])
	endfor
End

STATIC Function GetMeasurement(dW,m0,chn,tVar)
	Wave dW, m0
	Variable chn, tVar
	
	Duplicate/O/FREE dW, tempW
	tempW[] = (m0[p][%ch] == chn && m0[p][%frame] == tVar) ? dW[p] : NaN
	WaveTransform zapnans tempW
	if(numpnts(tempW) > 0)
		return mean(tempW)
	else
		return 0
	endif
End


STATIC Function LabelWaveDimensions(m0,dimVar,labelStr)
	Wave m0
	Variable dimVar
	String labelStr
	Variable nLabel = ItemsInList(labelStr)
	String str
	
	Variable i
	
	for(i = 0; i < nLabel; i += 1)
		str = StringFromList(i,labelStr)
		SetDimLabel dimVar,i,$str,m0
	endfor
End

STATIC Function/S GetDimLabelsFromWave(m0,dimVar)
	Wave m0
	Variable dimVar
	String labelStr = ""
	Variable nLabel = DimSize(m0,dimVar)
	
	Variable i
	
	for(i = 0; i < nLabel; i += 1)
		labelStr += GetDimLabel(m0,dimVar,i) + ";"
	endfor
	
	return labelStr
End

Function GenerateExamineWindow(dataset)
	Variable dataset
	Wave paramWave = root:paramWave
	
	SetDataFolder "root:data:cell" + num2str(dataset)
	
	KillWindow/Z examine
	Display/N=examine/W=(35,53,625,665)
	MakeDisplayWavesAndAdd()
//	SetAxis/W=examine bottom 0,200*paramWave[%xx]
//	SetAxis/W=examine left 200*paramWave[%xx],0
	SetAxis/W=examine bottom 0,200
	SetAxis/W=examine left 200,0
	ControlBar/W=examine 50
	ModifyGraph/W=examine width={Aspect,1}
	ModifyGraph/W=examine tick=3,mirror=1,noLabel=2,standoff=0
	ModifyGraph/W=examine margin=7
	SetWindow examine hook(myhook) = hook
	WAVE/Z coords_0
	WaveStats/Q/RMD=[][5] coords_0
	Variable maxT = V_Max
	Slider sld, limits={1,maxT - 1,1}, size={400,20}, vert=0, proc=ActionProcName
	
	UpdateFrame(0)
	
	return 0
End

Function updateFrame(liveSlice)
	Variable liveSlice
	DFREF dfr = $("root:exam:")
	Variable index = 0
	String tName
	
	do
		tName = GetIndexedObjNameDFR(dfr,1,index)
		if(strlen(tName) == 0)
			break
		endif
		if(stringmatch(tName,"live_*"))
			Wave lW = $("root:exam:" + tName)
			Wave w = $("root:exam:" + ReplaceString("live_",tName,""))
			lw[][] = (w[p][2] == liveSlice) ? w[p][q] : NaN
		endif
		index += 1
	while(1)
End

Function hook(s)
	STRUCT WMWinHookStruct &s

	Variable hookResult = 0
	switch(s.eventCode)
		// Could do something here with the displayed data - or delete this hook?
	endswitch
	
	return hookResult
End

Function ActionProcName(sa) : SliderControl
	STRUCT WMSliderAction &sa

	sa.blockReentry=1
	if(sa.eventcode == 9)
//		ModifyImage/W=examine Image plane=(sa.curval)
		updateFrame(sa.curval)
	endif
End

Function MakeDisplayWavesAndAdd()
	String wList = WaveList("coords_*",";","")
	Variable nWaves = ItemsInList(wList)
	PXPUtils#MakeColorWave(nWaves,"colorWave")
	WAVE/Z colorWave = root:colorWave
	NewDataFolder/O root:exam
	
	Variable i
	
	for(i = 0; i < nWaves; i += 1)
		Wave oW = $(StringFromList(i,wList))
		Duplicate/O/RMD=[][2,4] oW, $("root:exam:snk_" + num2str(i))
		Wave w = $("root:exam:snk_" + num2str(i))
		SetScale/P y 0,1,"", w
		w[][2] = oW[p][5]
		Duplicate/O w, $("root:exam:live_snk_" + num2str(i))
		Wave ls = $("root:exam:live_snk_" + num2str(i))
		AppendToGraph/W=examine ls[][1] vs ls [][0]
		if(i > 0)
			ModifyGraph/W=examine lstyle($("live_snk_" + num2str(i)))=1
		endif
		ModifyGraph/W=examine rgb($("live_snk_" + num2str(i)))=(colorWave[i][0],colorWave[i][1],colorWave[i][2])
	endfor
End

Function PlotVsDistance()
	SetDataFolder root:data:
	
	String windowList = "q_ch1VsDist;q_ch2VsDist;"
	MakeGraphWindows(windowList)
	
	String df1, wList, wName, newName
	Variable nCell, nWaves
	
	DFREF dfr0 = GetDataFolderDFR()
	nCell = CountObjectsDFR(dfr0, 4)
	
	Variable i,j
	
	for(i = 0; i < nCell; i += 1)
		df1 = GetIndexedObjNameDFR(dfr0, 4, i)
		SetDataFolder $("root:data:" + df1 + ":")
		wList = WaveList("distance_*",";","")
		nWaves = ItemsInList(wList)
		
		for(j = 0; j < nWaves; j += 1)
			wName = StringFromList(j,wList)
			// ch1
			Wave plateW = $("diff_ch1_0")
			Wave chrW = $(ReplaceString("distance_",wName,"diff_ch1_"))
			newName = ReplaceString("distance_",wName,"contrast_ch1_")
			MatrixOp/O $newName = chrW / plateW
			AppendToGraph/W=q_ch1VsDist $newName vs $wName
			
			// ch2
			Wave plateW = $("diff_ch2_0")
			Wave chrW = $(ReplaceString("distance_",wName,"diff_ch2_"))
			newName = ReplaceString("distance_",wName,"contrast_ch2_")
			MatrixOp/O $newName = chrW / plateW
			AppendToGraph/W=q_ch2VsDist $newName vs $wName
		endfor
	endfor
		
	for (String plotName : ListToTextWave(windowList,";"))
		ModifyGraph/W=$plotName mode=3,marker=19,msize=1,mrkThick=0
		Label/W=$plotName left "Chromosome / Plate (log\\B2\\M)"
		Label/W=$plotName bottom "Distance (µm)"
		SetAxis/A/N=1/E=2/W=$plotName left
		ModifyGraph/W=$plotName log(left)=2
	endfor
	// free wave
	Make/O/N=(3,3)/FREE colorWave = {{0,166,81},{237,28,36},{34,118,185}}
	colorWave *= 257 // convert to 16-bit
	MatrixTranspose colorWave
	ModifyGraph/W=q_ch1VsDist rgb=(colorWave[0][0],colorWave[0][1],colorWave[0][2],32768)
	ModifyGraph/W=q_ch2VsDist rgb=(colorWave[1][0],colorWave[1][1],colorWave[1][2],32768)
	
	SetDataFolder root:
End

Function AverageTheWaves(STRING chStr)
	SetDataFolder root:data:
	Wave paramWave = root:paramWave
	Variable tStep = paramWave[%tt]
		
	String df1, wList, wName, newName
	Variable nCell, nWaves
	
	String windowList = "q_X_con_raw;q_X_mis_raw;q_X_con_norm;q_X_mis_norm;"
	windowList = ReplaceString("_X_",windowList,"_"+chStr+"_")
	MakeGraphWindows(windowList)
	
	DFREF dfr0 = GetDataFolderDFR()
	nCell = CountObjectsDFR(dfr0, 4)
	Make/O/N=(nCell,2) $"root:congressionW"
	Wave congressionW = $"root:congressionW"
	
	Variable i,j
	
	// go into each folder and look for contrast_ch1 waves
	for(i = 0; i < nCell; i += 1)
		df1 = GetIndexedObjNameDFR(dfr0, 4, i)
		SetDataFolder $("root:data:" + df1 + ":")
		wList = WaveList("contrast_" + chStr + "_*",";","")
		nWaves = ItemsInList(wList)
		// 3 column classification wave. Col 0 = congressed?, col 1 = congression p from start, col 2 = congression p from anaphase
		Make/O/N=(nWaves,3)/FREE classW
		
		for(j = 0; j < nWaves; j += 1)
			wName = StringFromList(j,wList)
			newName = ReplaceString("contrast_",wName,"rev_")
			Duplicate/O $wName, $newName
			Wave w = $newName
			Reverse w
			if(w[0] == 0)
				// it congressed
				classW[j][0] = 1
				FindLevel/EDGE=1/Q w, 0
				classW[j][1] = V_levelX
				classW[j][2] = x2pnt(w,V_levelX) + 1
				DeletePoints 0, x2pnt(w,V_levelX), w
			else
				// it missegregated
				classW[j][] = 0
			endif
			// rescale from 0 and scale to (negative) time.
			SetScale/P x 0,-tStep,"", w
			// replace 0 "drop outs" with NaN
			w[1,*] = (w[p] == 0) ? NaN : w[p]
			// do the same for raw trace
			Wave contrastW = $wName
			contrastW[] = (contrastW[p] == 0) ? NaN : contrastW[p]
		endfor
		// we need to make a version of tReal_0 that is a countdown to anaphase, this will help make raw traces more comparable
		Duplicate/O $"tReal_0", $"aReal_0"
		Wave tR = $"tReal_0"
		Wave aR = $"aReal_0"
		aR[] = 0 - tR[numpnts(aR) - 1 - p]
		// only take the one that congressed last i.e. discard earliest congression point chromosomes
		if(nWaves == 1)
			if(classW[0][0] == 0)
				// misaligned
				AppendToGraph/W=$(StringFromList(1,windowList)) $wName vs $("aReal_0")
				AppendToGraph/W=$(StringFromList(3,windowList)) $newName
			else
				// congressed
				AppendToGraph/W=$(StringFromList(0,windowList)) $wName vs $("aReal_0")
				AppendToGraph/W=$(StringFromList(2,windowList)) $newName
			endif
			congressionW[i][] = classW[0][q * 2] // col 0 and 1 are 0 and 2
		elseif(nWaves > 1)
			WaveStats/M=1/RMD=[][1]/Q classW
			wName = StringFromList(V_minRowLoc,wList) // select correct wave
			newName = ReplaceString("contrast_",wName,"rev_")
			if(classW[V_minRowLoc][0] == 0)
				// misaligned
				AppendToGraph/W=$(StringFromList(1,windowList)) $wName vs $("aReal_0")
				AppendToGraph/W=$(StringFromList(3,windowList)) $newName
			else
				// congressed
				AppendToGraph/W=$(StringFromList(0,windowList)) $wName vs $("aReal_0")
				AppendToGraph/W=$(StringFromList(2,windowList)) $newName
			endif
			congressionW[i][] = classW[V_minRowLoc][q * 2]
		endif
	endfor
	
	SetDataFolder root:
	
	for (String plotName : ListToTextWave(windowList,";"))
		Label/W=$plotName left "Chromosome / Plate (log\\B2\\M)"
		Label/W=$plotName bottom "Time (min)"
		SetAxis/A/N=1/E=2/W=$plotName left
		ModifyGraph/W=$plotName log(left)=2
		if(cmpstr(chStr,"ch1") == 0)
			ModifyGraph/W=$plotName rgb=(0,166*257,81*257,32768)
		else
			ModifyGraph/W=$plotName rgb=(237*257,28*257,36*257,32768)
		endif
		AddMeanToGraph(plotName)
		if(StringMatch(plotName,"*_norm") == 1)
			SetAxis/W=$plotName bottom -60,0
		else
			SetAxis/W=$plotName bottom -120,0
		endif
	endfor
	
	plotName = "q_" + chStr + "_mean"
	KillWindow/Z $plotName
	Display/N=$plotName/HIDE=1
	AppendToGraph/W=$plotName $("ave_" + chStr + "_con_raw"),$("ave_" + chStr + "_mis_raw")
	ModifyGraph/W=$plotName lsize=2
	String tList = TraceNameList(plotName,";",1)
	if(cmpstr(chStr,"ch1") == 0)
		ModifyGraph/W=$plotName rgb($StringFromList(0,tList))=(0,166*257,81*257)
		ModifyGraph/W=$plotName rgb($StringFromList(1,tList))=(0,103*257,51*257)
	else
		ModifyGraph/W=$plotName rgb($StringFromList(0,tList))=(237*257,28*257,36*257)
		ModifyGraph/W=$plotName rgb($StringFromList(1,tList))=(152*257,0*257,8*257)
	endif
	ErrorBars $StringFromList(0,tList) SHADE= {0,0,(0,0,0,0),(0,0,0,0)},wave=($(ReplaceString("ave_",StringFromList(0,tList),"err_")),$(ReplaceString("ave_",StringFromList(0,tList),"err_")))
	ErrorBars $StringFromList(1,tList) SHADE= {0,0,(0,0,0,0),(0,0,0,0)},wave=($(ReplaceString("ave_",StringFromList(1,tList),"err_")),$(ReplaceString("ave_",StringFromList(1,tList),"err_")))
	ModifyGraph/W=$plotName log(left)=2
	SetAxis/W=$plotName bottom -70,0
	Label/W=$plotName left "Chromosome / Plate (log\\B2\\M)"
	Label/W=$plotName bottom "Time (min)"
End

// Igor9 style
Function MakeGraphWindows(String plotList)
   
	for (String plotName : ListToTextWave(plotList,";"))
		KillWindow/Z $plotName
		Display/N=$plotName/HIDE=1
	endfor
End

Function AddMeanToGraph(plotName)
	String plotName
	String xList = "", yList = ""
	Variable nTraces = ItemsInList(TraceNameList(plotName, ";", 1))
	
	Variable i
	
	for(i = 0; i < nTraces; i += 1)
		yList += GetWavesDataFolder(WaveRefIndexed(plotName, i, 1), 2) + ";"
		if(strlen(NameOfWave(WaveRefIndexed(plotName,i,2))) > 0)
			xList += GetWavesDataFolder(WaveRefIndexed(plotName, i, 2), 2) + ";"
		endif
	endfor
	String avName = ReplaceString("q_",plotName,"ave_")
	String errName = ReplaceString("q_",plotName,"err_")
	if(ItemsInList(yList) == ItemsInList(xList))
		fWaveAverage(yList, xList, 1, 1, avName, errName)
	else
		fWaveAverage(yList, "", 1, 1, avName, errName)
	endif
	AppendToGraph/W=$PlotName $avName		
	ErrorBars/W=$plotName $avName SHADE= {0,0,(0,0,0,0),(0,0,0,0)},wave=($ErrName,$ErrName)
	ModifyGraph/W=$plotName lsize($avName)=2,rgb($avName)=(0,0,0)
End

Function MakeDistancePlotsFromLastChromosomeTraces(plotName)
	String plotName
	SetDataFolder root:
	// get a list of traces from the original window
	String xList = "", yList = ""
	String tList = TraceNameList(plotName, ";", 1)
	Variable nTraces = ItemsInList(tList)
	
	Variable i
	
	for(i = 0; i < nTraces; i += 1)
		if(!stringmatch(StringFromList(i,tList),"ave_*"))
			yList += GetWavesDataFolder(WaveRefIndexed(plotName, i, 1), 2) + ";"
			if(strlen(NameOfWave(WaveRefIndexed(plotName,i,2))) > 0)
				xList += GetWavesDataFolder(WaveRefIndexed(plotName, i, 2), 2) + ";"
			endif
		endif
	endfor
	
	// these waves are all called contrast_ch1_* and aReal_0
	
	plotName = ReplaceString("ch1",plotName,"dist")
	KillWindow/Z $plotName
	Display/HIDE=1/N=$plotName
	
	nTraces = ItemsInList(yList)
	String wName,dName,newName
	
	for(i = 0; i < nTraces; i += 1)
		// these names are full path - we should be in root:
		wName = StringFromList(i,yList)
		// corresponding distance wave should be
		dName = ReplaceString("contrast_ch1",wName,"distance")
		// we will make
		newName = ReplaceString("ch1",wName,"dist")
		Duplicate/O $dName, $newName
		Wave w = $newName
		w[] = (w[p] == 0) ? NaN : w[p]
		AppendToGraph/W=$plotName w vs $(StringFromList(i,xList))
	endfor
	
	ModifyGraph/W=$plotName rgb=(0,0,0,32768)
	Label/W=$plotName left "Distance (µm)"
	Label/W=$plotName bottom "Time (min)"
	AddMeanToGraph(plotName)
End

Function CompareDistanceMeans()
	String plotName = "q_dist_mean"
	KillWindow/Z $plotName
	Display/N=$plotName/HIDE=1 $("ave_dist_con_raw"),$("ave_dist_mis_raw")
	ModifyGraph/W=$plotName lsize=2
	ModifyGraph/W=$plotName rgb($"ave_dist_con_raw")=(16448,16448,16448)
	ModifyGraph/W=$plotName rgb($"ave_dist_mis_raw")=(0,0,0)
	ErrorBars $"ave_dist_con_raw" SHADE= {0,0,(0,0,0,0),(0,0,0,0)},wave=($("err_dist_con_raw"),$("err_dist_con_raw"))
	ErrorBars $"ave_dist_mis_raw" SHADE= {0,0,(0,0,0,0),(0,0,0,0)},wave=($("err_dist_mis_raw"),$("err_dist_mis_raw"))
	SetAxis/W=$plotName/A/N=1/E=1 left
	SetAxis/W=$plotName bottom -70,0
	Label/W=$plotName left "Distance (µm)"
	Label/W=$plotName bottom "Time (min)"
End

Function SurvivalStyleDistancePlot()
	SetDataFolder root:
	String plotName = "q_dist_survival"
	KillWindow/Z $plotName
	Display/N=$plotName/HIDE=1
	WAVE/Z congressionW
	if(!WaveExists(congressionW))
		return -1
	endif
	WAVE/Z paramWave
	// two groups, congressed (con = 1), missegregated (mis = 0)
	Duplicate/O congressionW, congressionW_0,congressionW_1
	// con first
	congressionW_0[][] = (congressionW[p][0] == 0) ? congressionW[p][q] : NaN
	MatrixOp/O congressionW_0 = col(congressionW_0,1)
	WaveTransform zapnans congressionW_0
	// by definition cogressionW_0 will all = 0 and that is correct
	congressionW_1[][] = (congressionW[p][0] == 1) ? congressionW[p][q] : NaN
	MatrixOp/O congressionW_1 = col(congressionW_1,1)
	WaveTransform zapnans congressionW_1
	// data needs + 1 frame to be relative to anaphase, needs to be scaled to negative minutes
	congressionW_1[] = (congressionW_1[p] + 1) * -paramWave[%tt]
	Variable firstCon = WaveMin(congressionW_1) - 3
	Variable nBins = (0 - firstCon) / 3
	// histograms. 1 - cumulative, normed histo is a survival curve
	Make/N=(nBins)/O congressionW_0_Hist,congressionW_1_Hist
	Histogram/CUM/B={firstCon,3,nBins}/P congressionW_0,congressionW_0_Hist
	congressionW_0_hist[] = 1 - congressionW_0_hist[p]
	Histogram/CUM/B={firstCon,3,nBins}/P congressionW_1,congressionW_1_Hist
	congressionW_1_hist[] = 1 - congressionW_1_hist[p]
	AppendToGraph/W=$plotName congressionW_0_Hist,congressionW_1_Hist
	ModifyGraph/W=$plotName mode=6 // cityscape
	ModifyGraph/W=$plotName rgb(congressionW_0_hist)=(0,0,0)
	ModifyGraph/W=$plotName rgb(congressionW_1_hist)=(32768,32768,32768)
//	TextBox/W=$plotName/C/N=text0/F=0/A=LB/X=0.00/Y=0.00 "Congressed, " + num2str(numpnts(congressionW_0)) + "\rMisaligned, " + num2str(numpnts(congressionW_0))
	// formatting
	SetAxis/W=$plotName left 0,1
	SetAxis/W=$plotName bottom firstCon,0
	Label/W=$plotName left "Misaligned"
	Label/W=$plotName bottom "Time (min)"
End

Function MakeAndPlotFits(plotName)
	String plotName
	SetDataFolder root:
	// get a list of traces from the original window
	String xList = "", yList = ""
	String tList = TraceNameList(plotName, ";", 1)
	Variable nTraces = ItemsInList(tList)
	
	Variable i
	
	for(i = 0; i < nTraces; i += 1)
		if(!stringmatch(StringFromList(i,tList),"ave_*"))
			yList += GetWavesDataFolder(WaveRefIndexed(plotName, i, 1), 2) + ";"
			if(strlen(NameOfWave(WaveRefIndexed(plotName,i,2))) > 0)
				xList += GetWavesDataFolder(WaveRefIndexed(plotName, i, 2), 2) + ";"
			endif
		endif
	endfor
	
	plotName = ReplaceString("raw",plotName,"fit")
	KillWindow/Z $plotName
	Display/HIDE=1/N=$plotName
	
	String xName = ReplaceString("q_",plotName,"x_")
	String yName = ReplaceString("q_",plotName,"y_")
	Concatenate/O/NP=0 ylist, $yName
	Wave bigY = $yName
	Concatenate/O/NP=0 xlist, $xName
	Wave bigX = $xName
	bigX[] = (numtype(bigY[p]) == 2) ? NaN : bigX[p]
	WaveTransform zapnans bigX
	WaveTransform zapnans bigY
	// graph out the data
	AppendToGraph/W=$plotName bigY vs bigX
	Variable rV,gV,bV
	if(stringmatch(plotName,"*ch1*") == 1)
		rV=0;gV=166*257;bV=81*257
		if(stringmatch(plotName,"*con*") == 1)
			rV=0;gV=103*257;bV=51*257
		endif
	else
		rV=237*257;gV=28*257;bV=36*257
		if(stringmatch(plotName,"*con*") == 1)
			rV=152*257;gV=0*257;bV=8*257
		endif
	endif
	ModifyGraph/W=$plotName mode=2
	Label/W=$plotName left "Chromosome / Plate"
	Label/W=$plotName bottom "Time (min)"
	SetAxis/W=$plotName bottom -120,0
	CurveFit/Q line bigY /X=bigX /D /F={0.95, 5}
	ModifyGraph/W=$plotName rgb=(rV,gV,bV,32768)
	String fitName = "fit_" + yName
	String ucName = "UC_" + yName
	String lcName = "LC_" + yName
	ModifyGraph/W=$plotName mode($ucName)=7,hbFill($ucName)=4,toMode($ucName)=1
	ModifyGraph/W=$plotName lsize($ucName)=0
	ModifyGraph/W=$plotName lsize($lcName)=0
	ModifyGraph/W=$plotName rgb($fitName)=(rV,gV,bV,65535)
	ReorderTraces/W=$plotName $fitName,{$ucName,$lcName}
End