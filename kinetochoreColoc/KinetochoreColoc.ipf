#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include <3DWaveDisplay>
#include <ImageSlider>

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
	Make/O/N=(6,4) colorWave={{0.5,1,0,0,1,0},{0.5,1,0,0,1,1},{0.5,1,0,0,0,1},{0.75,0.75,0,0,0.75,0.75}}
	colorWave[][] = floor(65535*colorWave[p][q])
	AppendToGraph/T/W=examine liveYW/TN=kts vs liveXW
	ModifyGraph/W=examine mode(kts)=3,marker(kts)=8,msize(kts)=5,mrkThick(kts)=2
	ModifyGraph/W=examine rgb(kts)=(65535,65535,0,32768)
	ModifyGraph/W=examine zColor(kts)={ktCat,0,5,ctableRGB,0,colorWave}
	AppendToGraph/T/W=examine liveYW/TN=ktlabels vs liveXW
	ModifyGraph/W=examine mode(ktlabels)=3,textMarker(ktlabels)={ktNum,"default",0,0,2,-2.50,2.50}
	ModifyGraph/W=examine rgb(ktlabels)=(65535,65535,65535,32768)
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
	// do it button
	Button DoIt,pos={10,90+20*nRow},size={120,20},proc=CompleteButtonProc,title="Complete"
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
		// condition label
//		SetVariable $boxName0,pos={30+120*j,53+(mod(i,nRow)*20)},size={50,14},value= ktCat[i], title=" "
//		PopUpMenu $boxName0,pos={30+120*j,53+(mod(i,nRow)*20)},size={50,14},value= "Error;Aligned;MisalignedNotEnsheathed;MisalignedEnsheathed;", popvalue = StringFromList(ktCat[i],popupStr)
		String presetTo = StringFromList(ktCat[i],list)
		Variable mode = 1 + WhichListItem(presetTo,list) // about the same as 1+ktCat[i]
		PopUpMenu $boxName0,pos={30+140*j,53+(mod(i,nRow)*20)},size={90,14}, bodywidth=90,value= #popupStr, popvalue = presetTo, mode=mode, proc=ktPopupHelper
		ktVal = ktCat[i]
		valName0 = "val0_" + num2str(i)
		ValDisplay $valName0, pos={124+140*j,53+(mod(i,nRow)*20)}, size={18,18},bodyWidth=18,value=0,mode=1,barmisc={0,0}
		ValDisplay $valName0, limits={-1,1,0}, frame=0, zeroColor=(colorWave[ktCat[i]][0],colorWave[ktCat[i]][1],colorWave[ktCat[i]][2],colorWave[ktCat[i]][3])
//		SetDrawEnv fillfgc=(colorWave[ktVal][0],colorWave[ktVal][1],colorWave[ktVal][2],colorWave[ktVal][3])
//		DrawOval/W=ktCategory 80+120*j,50+(mod(i,nRow)*20),98+120*j,68+(mod(i,nRow)*20)
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

//// define buttons
//Function ButtonProc(ctrlName) : ButtonControl
//	String ctrlName
//
//	Wave/T PathWave1,PathWave2
//	Variable refnum, wNum, ii
//	String expr, wNumStr, iiStr, stringForTextWave
//
//	if(StringMatch(ctrlName,"file*") == 1)
//		expr="file([[:digit:]]+)\\w([[:digit:]]+)"
//		SplitString/E=(expr) ctrlName, wNumStr, iiStr
//		// get File Path
//		Open/D/R/F="*.xls*"/M="Select Excel Workbook" refNum
//		stringForTextWave = S_filename
//	else
//		expr="dir([[:digit:]]+)\\w([[:digit:]]+)"
//		SplitString/E=(expr) ctrlName, wNumStr, iiStr
//		// set outputfolder
//		NewPath/O/Q DirOfCSVs
//		PathInfo DirOfCSVs
//		stringForTextWave = S_Path
//	endif
//
//	if (strlen(stringForTextWave) == 0) // user pressed cancel
//		return -1
//	endif
//	wNum = str2num(wNumStr)
//	ii = str2num(iiStr)
//	if (wNum == 1)
//		PathWave1[ii] = stringForTextWave
//	else
//		PathWave2[ii] = stringForTextWave
//	endif
//End
//
//Function DoItButtonProc(ctrlName) : ButtonControl
//	String ctrlName
// 	
// 	WAVE/T CondWave
//	WAVE/T PathWave1
//	Variable okvar = 0
//	
//	strswitch(ctrlName)	
//		case "DoIt" :
//			// check CondWave
//			okvar = WaveChecker(CondWave)
//			if (okvar == -1)
//				DoAlert 0, "Error: Not all conditions have a name."
//				break
//			endif
//			okvar = NameChecker(CondWave)
//			if (okvar == -1)
//				DoAlert 0, "Error: Two conditions have the same name."
//				break
//			endif
//			okvar = WaveChecker(PathWave1)
//			if (okvar == -1)
//				DoAlert 0, "Error: Not all conditions have a file to load."
//				break
//			else
//				Migrate()
//			endif
//	endswitch	
//End
