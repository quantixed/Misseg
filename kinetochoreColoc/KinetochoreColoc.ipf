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
	WAVE/Z BX, BY, BZ
	WAVE/Z/T ktNum
	Duplicate/O BX, liveBX
	Duplicate/O BY, liveBY
	Duplicate/O BY, ktCat
	liveBX[] = (BZ[p] == 0+1) ? BX[p] : NaN
	liveBY[] = (BZ[p] == 0+1) ? BY[p] : NaN
	ktCat[] = 1
	Make/O/N=(6,4) colorWave={{0.5,1,0,0,1,1},{0.5,1,0,0,1,0},{0.5,1,0,0,0,1},{0.5,0.5,0,0,0.5,0.5}}
	colorWave[][] = floor(65535*colorWave[p][q])
	AppendToGraph/T/W=examine liveBY/TN=kts vs liveBX
	ModifyGraph/W=examine mode(kts)=3,marker(kts)=8,msize(kts)=5
	ModifyGraph/W=examine rgb(kts)=(65535,65535,0,32768)
	ModifyGraph/W=examine zColor(kts)={ktCat,0,5,ctableRGB,0,colorWave}
	AppendToGraph/T/W=examine liveBY/TN=ktlabels vs liveBX
	ModifyGraph/W=examine mode(ktlabels)=3,textMarker(ktlabels)={ktNum,"default",0,0,2,-2.50,2.50}
	ModifyGraph/W=examine rgb(ktlabels)=(65535,65535,65535,32768)
End


Function updateFrame(liveSlice)
	Variable liveSlice
	WAVE/Z liveBX, liveBY, BX, BY, BZ
	liveBX[] = (BZ[p] == liveSlice + 1) ? BX[p] : NaN
	liveBY[] = (BZ[p] == liveSlice + 1) ? BY[p] : NaN
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
