#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include <3DWaveDisplay>
#include <ImageSlider>
#include "Misseg"

////////////////////////////////////////////////////////////////////////
// Menu items
////////////////////////////////////////////////////////////////////////
Menu "Misseg"
	Submenu "Kinetochore Coloc"
		"Classify Kinetochores...", /Q, KinetochoreClassification()
		"Pick another image", /Q, DisplayUpdatedSelector()
		"Coloc Analysis...", /Q, KinetochoreColoc()
	End
End


////////////////////////////////////////////////////////////////////////
// Master functions and wrappers
////////////////////////////////////////////////////////////////////////
Function KinetochoreClassification()
	SetUpFilesToLoad()
	DisplayUpdatedSelector()
End

Function KinetochoreColoc()
	Misseg#CleanSlate()
	if(LoadCSVFiles() == -1)
		return -1
	endif
	NNCalcs()
End

////////////////////////////////////////////////////////////////////////
// Main functions
////////////////////////////////////////////////////////////////////////

Function SetUpFilesToLoad()
	SetDataFolder root:
	// first, locate the TIFF images of the D3D.dv files
	NewPath/O/Q/M="Locate TIFF images" imgDiskFolder
	if (V_flag != 0)
		DoAlert 0, "Disk folder error"
		return -1
	endif
	PathInfo/S imgDiskFolder
	Make/O/N=2/T pathWave
	pathWave[0] = S_Path
	
	// now, locate the folder containing the csv files from 3DOC in Fiji
	NewPath/O/Q/M="Locate csv files (outputs will be saved here)" csvDiskFolder
	if (V_flag != 0)
		DoAlert 0, "Disk folder error"
		return -1
	endif
	PathInfo/S csvDiskFolder
	pathWave[1] = S_Path
	
	NewDataFolder/O root:data
End

Function DisplayUpdatedSelector()
	SetDataFolder root:
	if(ItemsInList(PathList("*DiskFolder",";","")) < 2)
		DoAlert 0, "Please load the file lists first"
		return -1
	endif
	// update filelists - all images
	String imgList = IndexedFile(imgDiskFolder,-1,".tif")
	// all csv files - each time this is run, it will grow
	// since images analysed with have results saved here
	String csvList = IndexedFile(csvDiskFolder,-1,"TEXT")
	if (ItemsInList(csvList) < 3)
		DoAlert 0, "Found less than three text files"
		return -1
	endif
	// now we remove the items that have results already
	String modImgList = imgList
	String imgName, originalFileName
	Variable nFiles = ItemsInList(imgList)
	
	Variable i,j
	
	for(i = 0; i < nFiles; i += 1)
		imgName = StringFromList(i,imgList)
		originalFileName = ReplaceString(".tif",imgName,"")
		if(FindListItem(originalFileName + "_ktCat.txt", csvList) >= 0)
			modImgList = RemoveFromList(imgName,modImgList)
		endif
	endfor
	// new number of images to analyse, i.e. to go in listbox
	nFiles = ItemsInList(modImgList)
	if(nFiles == 0)
		DoAlert 0, "No images left to classify."
		return -1
	endif
	Wave/T w0 = ListToTextWave(modImgList,";")
	Make/O/N=(nFiles)/T imgFileList,imgNumList
	imgFileList[] = w0[p]
	imgNumList[] = num2str(p)
	
	// build 2D wave for listbox
	Concatenate/O/NP=1 {imgNumList,imgFileList}, selectTextWave
	Make/O/B/N=(nFiles) selectedWave = 0
	Wave selectTextWave = root:selectTextWave
	
	Wave/T w1 = ListToTextWave(csvList,";")
	Make/O/N=(ItemsInList(csvList))/T csvFileList
	csvFileList[] = w1[p]

	BuildSelectorPanel()
End

Function LoadCSVFiles()
	SetDataFolder root:
	// locate the folder containing the csv files from 3DOC in Fiji
	NewPath/O/Q/M="Locate csv files (outputs will be saved here)" csvDiskFolder
	if (V_flag != 0)
		DoAlert 0, "Disk folder error"
		return -1
	endif
	PathInfo/S csvDiskFolder
	Make/O/N=2/T pathWave
	pathWave[1] = "" // do this for consistency with previous code
	
	// each image has four text files
	// 1 ktsMod.csv
	// 2 knstrn.csv
	// 3 ktCat.txt
	// 4 log.txt
	// we need the first 3 and not the 4th
	
	// first run a check on all the files
	String fileList = IndexedFile(csvDiskFolder,-1,"TEXT")
	if (ItemsInList(fileList) < 3)
		DoAlert 0, "Found less than three text files"
		return -1
	endif
	Variable nFiles = ItemsInList(fileList)
	String ktCatList = ListMatch(fileList,"*_ktCat.txt")
	Variable nCatFiles = ItemsInList(ktCatList)
		
	Print "Checking categorisation files:", nCatFiles, "found."
	Print "There should be", ItemsInList(ListMatch(fileList,"*_ktsMod.csv"))
	
	// check and then load valid triplets
	String fileName, originalFileName
	String validKtCatList = ktCatList
	
	Variable i
	
	for(i = 0; i < nCatFiles; i += 1)
		fileName = StringFromList(i,ktCatList)
		originalFileName = ReplaceString("_ktCat.txt",fileName,"")
		if(FindListItem(originalFileName + "_ktsMod.csv", fileList) == -1)
			Print originalFileName, "missing ktsMod file"
			validKtCatList = RemoveFromList(fileName, validKtCatList)
		endif
		if(FindListItem(originalFileName + "_knstrn.csv", fileList) == -1)
			Print originalFileName, "missing knstrn file"
			validKtCatList = RemoveFromList(fileName, validKtCatList)
		endif
	endfor
	
	NewDataFolder/O/S root:data
	nCatFiles = ItemsInList(validKtCatList)
	Wave/T w0 = ListToTextWave(validKtCatList,";")
	Make/O/N=(nCatFiles)/T root:ktCatNameWave
	Wave/T ktCatNameWave = root:ktCatNameWave
	ktCatNameWave[] = w0[p]
	
	String thisFile,wName
	Variable err = GetRTError(1)

	// now we will load everything in
	for(i = 0; i < nCatFiles; i += 1)
		fileName = ktCatNameWave[i]
		thisFile = fileName
		wName = "ktCat_" + num2str(i) + "_M"
		if(CheckForData("csvDiskFolder", thisFile) > 1)
			LoadWave/M/O/J/K=1/L={0,1,0,0,0}/P=csvDiskFolder/A=$wName/Q thisFile
		else
			Make/O/N=(1,4) $(wName + "0")
		endif
		
		originalFileName = ReplaceString("_ktCat.txt",fileName,"")
		thisfile = originalFileName + "_ktsMod.csv"
		wName = "ktsMod_" + num2str(i) + "_M"
		if(CheckForData("csvDiskFolder", thisFile) > 1)
			LoadWave/M/O/J/K=1/L={0,1,0,0,0}/P=csvDiskFolder/A=$wName/Q thisFile
		else
			Make/O/N=(1,25) $(wName + "0")
		endif
		
		thisfile = originalFileName + "_knstrn.csv"
		wName = "knstrn_" + num2str(i) + "_M"
		if(CheckForData("csvDiskFolder", thisFile) > 1)
			LoadWave/M/O/J/K=1/L={0,1,0,0,0}/P=csvDiskFolder/A=$wName/Q thisFile
		else
			Make/O/N=(1,25) $(wName + "0")
		endif
		CSVChecker(i)
	endfor
	
	SetDataFolder root:	
End

Function NNCalcs()
	Wave/Z/T ktCatNameWave = root:ktCatNameWave
	Variable nWave = numpnts(ktCatNameWave)
	
	Make/O/N=(nWave,7) root:nKtWave
	Wave nKtWave = root:nKtWave
	String wName
	SetDataFolder root:data
	
	Variable i
	
	for(i = 0; i < nWave; i += 1)
		wName = "ktCat_" + num2str(i) + "_M0"
		Wave w0 = $wName
		nKtWave[i][] = GetCategories(w0,q) // 7th column will be set to 0
		nKtWave[i][6] = DimSize(w0,0)
		if((nKtWave[i][0] + nKtWave[i][1] + nKtWave[i][2] + nKtWave[i][3] + nKtWave[i][4] + nKtWave[i][5]) != nKtWave[i][6])
			Print "dataset", i, "mismatch in classification"
		endif
		// check if there are no Kt to analyse
		if(nKtWave[i][0] == nKtWave[i][6])
			continue
		else
			// 10,11,12 are x,y,z in ktsMod (w1) and knstrn (w2)
			NNWrapper(i,1)
			NNWrapper(i,4)
			NNWrapper(i,5)
		endif
	endfor
	
	SetDataFolder root:
End

STATIC Function NNWrapper(ii,cat)
	Variable ii,cat
	
	Wave nKtWave = root:nKtWave
	String wName = "nnDist_" + num2str(ii) + "_" + num2str(cat)
	Make/O/N=(nKtWave[ii][cat]) $wName
	Wave nnW = $wName
	
	Wave ktCatW = $("ktCat_" + num2str(ii) + "_M0")
	Wave ktsModW = $("ktsMod_" + num2str(ii) + "_M0")
	Wave knstrnW = $("knstrn_" + num2str(ii) + "_M0")
	if(DimSize(knstrnW,0) < 1)
		nnW[] = -1
		return -1
	endif
	
	Variable nRows = DimSize(ktCatW,0)
	Make/O/N=(1,3)/FREE ktMat	// 1 row matrix to hold x y z coords of reference kt
	
	Variable i, j = 0
	
	for(i = 0; i < nRows; i += 1)
		if(ktCatW[i][3] == cat)
			ktMat[0][] = ktsModW[i][10 + q]
			nnW[j] = Misseg#FindMinMaxOfTwo2DWaves(knstrnW,ktMat,0)
			j += 1
		endif
	endfor
End

////////////////////////////////////////////////////////////////////////
// Image functions
////////////////////////////////////////////////////////////////////////

Function LoadImageAndCSV(thisImg)
	String thisImg
	
	//we will store the currentFileName in a textwave in root
	Make/O/N=(1)/T currentFileName
	Wave/T currentFileName = root:currentFileName
	SetDataFolder root:data:
	// setup paths again? seems I don't need to
	String originalFileName = ReplaceString(".tif",thisImg,"")
	// load image and call it "image" in data folder
	currentFileName[0] = originalFileName
	ImageLoad/P=imgDiskFolder/T=tiff/O/S=0/C=-1/LR3D/Q/N=image thisImg
	// load csv file
	String thisCSV = originalFileName + "_ktsMod.csv"
	LoadWave/W/A/J/O/K=1/G/L={0,0,0,0,0}/P=csvDiskFolder thisCSV
	// this doesn't check if thisCSV exists - could cause error
	GenerateExamineWindow()
End

Function GenerateExamineWindow()
//	SetDataFolder root:data:
	WAVE/Z Image
	KillWindow/Z examine
	NewImage/N=examine Image
	ControlBar/W=examine 50
	ModifyGraph/W=examine width={Aspect,1}
	SetWindow examine hook(myhook) = hook
	Slider sld, limits={0,DimSize(Image,3)-1,1}, size={400,20},vert=0,proc=ActionProcName
	WAVE/Z XW, YW, ZW
	Variable nKts = numpnts(YW)
	Make/O/N=(nKts) liveXW, liveYW, ktCat
	Make/O/N=(nKts)/T ktNum = num2str(p + 1)
	liveXW[] = (floor(ZW[p]) == 0+1) ? XW[p] : NaN
	liveYW[] = (floor(ZW[p]) == 0+1) ? YW[p] : NaN
	ktCat[] = 1
//	Make/O/N=(6,4) colorWave={{0.5,0.99,0,0,0.352941,1},{0.5,0.99,0,0,0.376471,0.443137},{0.5,0.99,0,0,1,0.0941176},{1,1,1,1,1,1}}
	Make/O/N=(4,6) colorWave={{0.5,0.5,0.5,1},{0.99,0.99,0.99,1},{0,0,0,1},{0,0,0,1},{0,0.9,1,1},{1,0.5,0,1}}
	colorWave[][] = floor(65535*colorWave[p][q])
	MatrixTranspose colorWave
	AppendToGraph/T/W=examine liveYW/TN=kts vs liveXW
	ModifyGraph/W=examine mode(kts)=3,marker(kts)=8,msize(kts)=5,mrkThick(kts)=2
	ModifyGraph/W=examine rgb(kts)=(65535,65535,0,32768)
	ModifyGraph/W=examine zColor(kts)={ktCat,0,5,ctableRGB,0,colorWave}
	AppendToGraph/T/W=examine liveYW/TN=ktlabels vs liveXW
	ModifyGraph/W=examine mode(ktlabels)=3,textMarker(ktlabels)={ktNum,"default",0,0,2,-2.50,2.50}
	ModifyGraph/W=examine rgb(ktlabels)=(65535,65535,65535,49151)
	ktCategoryPanel(nKts)
	return 0
End

Function updateFrame(liveSlice)
	Variable liveSlice
	WAVE/Z liveXW, liveYW, XW, YW, ZW
	liveXW[] = (floor(ZW[p]) == liveSlice + 1) ? XW[p] : NaN
	liveYW[] = (floor(ZW[p]) == liveSlice + 1) ? YW[p] : NaN
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
		ModifyImage/W=examine Image plane=(sa.curval)
		updateFrame(sa.curval)
	endif
End


////////////////////////////////////////////////////////////////////////
// Panel functions
////////////////////////////////////////////////////////////////////////

Function BuildSelectorPanel()
//	DoWindow/F fileSelector	// bring panel to front if it exists
//	if(V_Flag != 0)
//		return 0							// panel already exists
//	endif
	KillWindow/Z fileSelector
	NewPanel/W=(280,56,580,532)/N=fileSelector
	DefaultGUIFont/W=fileSelector/Mac popup={"_IgorSmall",0,0}
	ListBox lb1,pos={34,13},size={233,400},frame=2,listWave=root:selectTextWave
	ListBox lb1,selWave=root:selectedWave, mode= 3,editStyle= 1
	ListBox lb1,widths= {25,400}
	
	Button Classify,pos={34,430},size={130,20},proc=ClassifyButtonProc,title="Classify"
	return 0
End

Function ClassifyButtonProc(ctrlName) : ButtonControl
	String ctrlName
 	
 	String imgName
 	WAVE/Z selectedWave = root:selectedWave
 	WAVE/Z/T selectTextWave = root:selectTextWave
 	FindValue/I=1 selectedWave
 	if(V_Value == -1)
 		return -1
 	else
 		imgName = selectTextWave[V_Value][1]
 	endif
	
	strswitch(ctrlName)
		case "Classify" :
			KillWindow/Z fileSelector
			LoadImageAndCSV(imgName)
			return 0
	endswitch
End

///	@param	nKts	number of kinetochores
Function ktCategoryPanel(nKts)
	Variable nKts
	// 30 rows, so how many columns
	Variable nRow = 30
	Variable nCol = ceil(nKts / nRow)
	
	WAVE/Z colorWave, ktCat
	KillWindow/Z ktCategory
	NewPanel/N=ktCategory/K=1/W=(0,0,20+140*nCol,150+20*nRow)/HOST=examine/EXT=0
	// labelling of columns
	DrawText/W=ktCategory 10,30,"Kinetochore category"
	// Four buttons
	Button Skip,pos={10,60+20*nRow},size={130,20},proc=CompleteButtonProc,title="Exit - no save"
	Button Reset,pos={10,82+20*nRow},size={130,20},proc=CompleteButtonProc,title="Reset kinetochores"
	Button Pass,pos={10,104+20*nRow},size={130,20},proc=CompleteButtonProc,title="Don't analyse"
	Button Complete,pos={10,126+20*nRow},fColor=(45489,65535,45489),size={130,20},proc=CompleteButtonProc,title="Complete - save"
	// add some help
	Button Skip,help={"Exits the image without saving.\rDo this if you want to take a break."}
	Button Reset,help={"Reset values of all kinetochores to 1 so that you can start over."}
	Button Pass,help={"If the image is not a good cell (not at right stage, kinetochores bad etc.)\rclick this so that the image is not analysed."}
	Button Complete,help={"When you have finished the classification, click here to save.\rImage will not be available for reclassification."}
	// insert rows
	String boxName0, valname0
	String list = "Error;Align;Poles;BG;Misalign;Misalign-ER;"
	String popupStr = "\""+list+"\""
	Variable ktVal
	Variable i,j
	
	for(i = 0; i < nKts; i += 1)
		j = floor(i / nRow)
		boxName0 = "box0_" + num2str(i)
		// row label
		DrawText/W=ktCategory 10+140*j,68+(mod(i,nRow)*20),num2str(i+1)
		// make the popups - thanks to chozo for this code
		String presetTo = StringFromList(ktCat[i],list)
		Variable mode = 1 + WhichListItem(presetTo,list) // about the same as 1+ktCat[i]
		PopUpMenu $boxName0,pos={30+140*j,53+(mod(i,nRow)*20)},size={90,14}, bodywidth=90,value= #popupStr, popvalue = presetTo, mode=mode, proc=ktPopupHelper
		ktVal = ktCat[i]
		valName0 = "val0_" + num2str(i)
		ValDisplay $valName0, pos={124+140*j,53+(mod(i,nRow)*20)}, size={18,18},bodyWidth=18,value=0,mode=1,barmisc={0,0}
		ValDisplay $valName0, limits={-1,1,0}, frame=0, zeroColor=(colorWave[ktCat[i]][0],colorWave[ktCat[i]][1],colorWave[ktCat[i]][2],colorWave[ktCat[i]][3])
	endfor
	
	return 0
End

Function ktPopupHelper(s) : PopupMenuControl
	STRUCT WMPopupAction &s

	if (s.eventCode == 2)   // mouse up
		String name = s.ctrlName
		Variable boxnum = str2num(name[5,inf])  // assumes the format to be box0_num
		Variable sel = s.popNum-1	// 1-6 (popup) vs. 0-5 (wave)
		
		WAVE/Z colorWave, ktCat
		ktCat[boxnum] = sel
		
		ValDisplay $("val0_" + num2str(boxnum)) ,win=ktCategory ,zeroColor=(colorWave[sel][0],colorWave[sel][1],colorWave[sel][2],colorWave[sel][3])
	endif
	return 0
End

Function CompleteButtonProc(ctrlName) : ButtonControl
	String ctrlName
 	
 	WAVE/Z XW, YW, ZW, ktCat
 	WAVE/T/Z currentFileName = root:currentFileName
 	String txtName = currentFileName[0] + "_ktCat.txt"
	
	strswitch(ctrlName)
		case "Skip" :
			KillWindow/Z examine
			KillWindow/Z ktCategory
			return 0
		case "Reset" :
			ktCat[] = 1
			ktCategoryPanel(numpnts(ktCat))
			break	
		case "Pass" :
			ktCat[] = 0
			Save/J/P=csvDiskFolder/W XW,YW,ZW,ktCat as txtName
			KillWindow/Z examine
			KillWindow/Z ktCategory
			return 1	
		case "Complete" :
			Save/J/P=csvDiskFolder/W XW,YW,ZW,ktCat as txtName
			KillWindow/Z examine
			KillWindow/Z ktCategory
			return 1
	endswitch
End


////////////////////////////////////////////////////////////////////////
// Utility functions
////////////////////////////////////////////////////////////////////////

STATIC Function CheckForData(pathName, filePath)
	// header (column names) are 243 characters
	// check that first data row is not blank
	String pathName	 // Name of symbolic path or ""
	String filePath		 // Name of file or partial path relative to symbolic path.

	Variable refNum

	Open/R/P=$pathName refNum as filePath
	
	String buffer, text
	Variable line = 0
	
	do
		FReadLine refNum, buffer
		if (strlen(buffer) == 0)
				return line
		endif
		line += 1
	while(line < 2)
	
	return line	   // return 2 if there was at least 1 data row
End

STATIC Function GetCategories(w,j)
	Wave w
	Variable j
	
	Duplicate/FREE/RMD=[][3] w,temp
	temp[] = (temp[p] == j) ? 1 : 0
	return sum(temp)
End

STATIC Function CSVChecker(ii)
	Variable ii
	
	Wave ktCatW = $("ktCat_" + num2str(ii) + "_M0")
	Wave ktsModW = $("ktsMod_" + num2str(ii) + "_M0")
	Wave knstrnW = $("knstrn_" + num2str(ii) + "_M0")
//	Print ii, ":", DimSize(ktCatW,0), DimSize(ktsModW,0),DimSize(knstrnW,0)
	if(DimSize(ktCatW,0) != DimSize(ktsModW,0))
		Print ii, ":", DimSize(ktCatW,0), DimSize(ktsModW,0), DimSize(knstrnW,0)
	endif
	
	return 0
End