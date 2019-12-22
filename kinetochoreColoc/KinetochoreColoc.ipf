#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include <3DWaveDisplay>
#include <ImageSlider>

////////////////////////////////////////////////////////////////////////
// Menu items
////////////////////////////////////////////////////////////////////////
Menu "Macros"
	"Classify Kinetochores...", /Q, KinetochoreClassification()
End


////////////////////////////////////////////////////////////////////////
// Master functions and wrappers
////////////////////////////////////////////////////////////////////////
Function KinetochoreClassification()
	LoadFilesAndClassify()
End

////////////////////////////////////////////////////////////////////////
// Main functions
////////////////////////////////////////////////////////////////////////
Function LoadFilesAndClassify()
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
	String imgList = IndexedFile(expDiskFolder,-1,".tif")
	Wave/T w0 = ListToTextWave(imgList,";")
	Variable nFiles = ItemsInList(imgList)
	Make/O/N=(nFiles)/T imgFileList
	imgFileList[] = w0[p]
	// now, locate the folder containing the csv files from 3DOC in Fiji
	NewPath/O/Q/M="Locate csv files (outputs will be saved here)" csvDiskFolder
	if (V_flag != 0)
		DoAlert 0, "Disk folder error"
		return -1
	endif
	PathInfo/S csvDiskFolder
	pathWave[1] = S_Path
	String csvList = IndexedFile(expDiskFolder,-1,"TEXT")
	Wave/T w1 = ListToTextWave(imgList,";")
	Make/O/N=(ItemsInList(csvList))/T csvFileList
	csvFileList[] = w1[p]
	
	String thisImg, originalFileName
	Variable i
	
	NewDataFolder/O/S root:data
	
	for (i = 0; i < nFiles; i += 1)
		thisImg = StringFromList(i, imgList)
		originalFileName = ReplaceString(".tif",ThisFile,"")
		// load image and call it image in data folder
		// load csv file
//		LoadWave/A/J/D/O/K=1/V={" "," $",0,0}/L={0,0,0,1,0}/P=expDiskFolder ThisFile
		GenerateExamineWindow()
		// how do we do flow control here because we need to wait for the result before continuing the loop
		SetDataFolder root:data:
	endfor
End

Function GenerateExamineWindow()
	WAVE/Z Image
	KillWindow/Z examine
	NewImage/N=examine Image
	ControlBar/W=examine 50
	ModifyGraph/W=examine width={Aspect,1}
	SetWindow examine hook(myhook) = hook
	Slider sld, limits={0,DimSize(Image,3)-1,1}, size={400,20},vert=0,proc=ActionProcName
	WAVE/Z XW, YW, ZW
	WAVE/Z/T ktNum
	Variable nKts = numpnts(YW)
	Make/O/N=(nKts) liveXW, liveYW, ktCat
	liveXW[] = (floor(ZW[p]) == 0+1) ? XW[p] : NaN
	liveYW[] = (floor(ZW[p]) == 0+1) ? YW[p] : NaN
	ktCat[] = 1
	Make/O/N=(6,4) colorWave={{0.5,0.99,0,0,0.352941,1},{0.5,0.99,0,0,0.376471,0.443137},{0.5,0.99,0,0,1,0.0941176},{1,1,1,1,1,1}}
	colorWave[][] = floor(65535*colorWave[p][q])
	AppendToGraph/T/W=examine liveYW/TN=kts vs liveXW
	ModifyGraph/W=examine mode(kts)=3,marker(kts)=8,msize(kts)=5,mrkThick(kts)=2
	ModifyGraph/W=examine rgb(kts)=(65535,65535,0,32768)
	ModifyGraph/W=examine zColor(kts)={ktCat,0,5,ctableRGB,0,colorWave}
	AppendToGraph/T/W=examine liveYW/TN=ktlabels vs liveXW
	ModifyGraph/W=examine mode(ktlabels)=3,textMarker(ktlabels)={ktNum,"default",0,0,2,-2.50,2.50}
	ModifyGraph/W=examine rgb(ktlabels)=(65535,65535,65535,49151)
	ktCategoryPanel(nKts)
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
		//The real code also contains some hooks for working with the displayed data
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
///////////////////////////////////////////////////////////////////////
///	@param	nKts	number of kinetochores
Function ktCategoryPanel(nKts)
	Variable nKts
	// 30 rows, so how many columns
	Variable nRow = 30
	Variable nCol = ceil(nKts / nRow)
	
	WAVE/Z colorWave, ktCat
	KillWindow/Z ktCategory
	NewPanel/N=ktCategory/K=1/W=(0,0,20+140*nCol,130+20*nRow)/HOST=examine/EXT=0
	// labelling of columns
	DrawText/W=ktCategory 10,30,"Kinetochore category"
	// Two buttons
	Button Reset,pos={10,65+20*nRow},size={80,20},proc=CompleteButtonProc,title="Reset"
	Button Complete,pos={10,90+20*nRow},size={120,20},proc=CompleteButtonProc,title="Complete"
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
End

Function ktPopupHelper(s) : PopupMenuControl
    STRUCT WMPopupAction &s
    
    if (s.eventCode == 2)   // mouse up
        String name = s.ctrlName
        Variable boxnum = str2num(name[5,inf])  // assumes the format to be box0_num
        Variable sel        = s.popNum-1        // 1-6 (popup) vs. 0-5 (wave)
        
        WAVE/Z colorWave, ktCat
        ktCat[boxnum] = sel
        
        ValDisplay $("val0_" + num2str(boxnum)) ,win=ktCategory ,zeroColor=(colorWave[sel][0],colorWave[sel][1],colorWave[sel][2],colorWave[sel][3])
    endif
    return 0
End

Function CompleteButtonProc(ctrlName) : ButtonControl
	String ctrlName
 	
 	WAVE/Z XW, YW, ZW, ktCat
 	WAVE/T/Z PathWave = root:PathWave
 	WAVE/T/Z currentFileName = root:currentFileName
 	NewPath/O/Q txtDiskFolder, PathWave[1]
 	String txtName = currentFileName[2] + "_ktCat.txt"
	
	strswitch(ctrlName)	
		case "Complete" :
			Save/J/P=txtDiskFolder/W XW,XW,YW,ktCat as txtName
			KillWindow/Z examine
			KillWindow/Z ktCategory
		case "Reset" :
			ktCat[] = 1
			ktCategoryPanel(numpnts(ktCat))
	endswitch
End
