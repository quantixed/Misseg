#pragma TextEncoding = "MacRoman"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// Front-end dialog to allow user to pick XML file and also to pick image - optional
// Give the option to define which counter point is which
// Give the option to change the scaling of the image
Menu "Macros"
	"Read XML...", /Q, WorkflowForXMLAnalysis()
End

Function WorkflowForXMLAnalysis()
	CleanSlate()
	SetupParamWaves()
	ReadXML("","")
	
End

// Needs XMLUtils XOP
Function ReadXML(pathName,fileName)
	String pathName		// Name of Igor symbolic path or "" to get dialog
	String fileName		// Name of file to load or "" to get dialog
	
	if ((strlen(pathName)==0) || (strlen(fileName)==0))
		// Display dialog looking for file.
		Variable refNum
		String filters = "XML Files (*.xml):.xml;"
		filters += "All Files:.*;"
		Open/D/R/P=$pathName/F=filters/M="Select XML file" refNum as fileName
		fileName = S_fileName			// S_fileName is set by Open/D
		if (strlen(fileName) == 0)		// User cancelled?
			return -2
		endif
	endif
	
	Variable fileID
	fileID = XMLopenfile(Pathname + fileName)
	xmlelemlist(fileID)
	WAVE/Z/T W_ElementList
	String XPathStr0,XpathStr1,resultMat
	
	Variable nMarkerType = 8
	Variable i,j,counter
	
	for(i = 0; i < nMarkerType; i += 1)
		resultMat = "mat_" + num2str(i+1)
		Make/O/N=(3) $resultMat
		Wave m0 = $resultMat
		XPathStr0 = "/CellCounter_Marker_File/Marker_Data/Marker_Type[" + num2str(i+1) + "]/Marker[nnn]"
		j = 1
		counter = 0
		do
			XPathStr1 = ReplaceString("nnn",XPathStr0,num2str(j))
			XMLwaveFmXpath(fileID,XPathStr1,"","")
			WAVE/Z/T M_XMLcontent
			if(numpnts(M_XMLcontent) != 0)
				if(j == 1)
					m0[][0] = str2num(M_XMLcontent[p+1])
				else
					Make/O/FREE/N=(3) tempW
					tempW[] = str2num(M_XMLcontent[p+1])
					Concatenate {tempW}, m0
				endif
			else
				counter += 1
			endif
			j += 1
		while (counter < 2)
		MatrixTranspose m0
	endfor
	xmlclosefile(fileID,0)
end

// Function to get process a point cloud (3-column 2D wav) with disc of points
///	@param	m0	matrix of points
Function PointDisc(m0)
	String m0 //prefix of waves
	Variable nPoints = dimsize(m0,0)
	
//	String xS = root + "_x"
//	String yS = root + "_y"
//	String zS = root + "_z"
//	String xList=wavelist(xs,";","")
//	String yList=wavelist(ys,";","")
//	String zList=wavelist(zs,";","")
//	String xWave,yWave,zWave,wspName,wsnName,wList,newsnName,newspName
	Make/O/N=(1,3) centroidW
	Variable dd,i,j,k,l,dist,kts
	
	for (i = 0; i < ItemsInList(nPoints); i += 1)
		MatrixOp/O xW = col(m0,0)
		MatrixOp/O yW = col(m0,1)
		MatrixOp/O zW = col(m0,2)
		// find centroid
		centroid[0][0] = mean(xW)
		centroid[0][1] = mean(yW)
		centroid[0][2] = mean(zW)
		//now do the PCA
		PCA/ALL/SEVC/SRMT/SCMT xW,yW,zW
		WAVE/Z M_C	//declare M_C wave from PCA /SCMT flag
		// because 0 = ax+by+cz+d where 
		dd = -((M_C[2][0] * centroid[0][0]) + (M_C[2][1] * centroid[0][1]) + (M_C[2][2] * centroid[0][2]))
		for(l = 0;l < nPoints;l +=1)	// use l for loop var
			//find distance to plane
			dist = ((M_C[2][0] * xw[l]) + (M_C[2][1] * yw[l]) + (M_C[2][2] * zw[l])) + d
			if (dist>=0)
				wsn[l][0]=nan
				wsn[l][1]=nan
				wsn[l][2]=nan
			else
				wsp[l][0]=nan
				wsp[l][1]=nan
				wsp[l][2]=nan
			endif
		endfor
		matrixop/o c1=col(wsn,0)
		matrixop/o c2=col(wsn,1)
		matrixop/o c3=col(wsn,2)
		WaveTransform zapnans c1
		WaveTransform zapnans c2
		WaveTransform zapnans c3
		Concatenate /o {c1,c2,c3}, $wsnName
		// now use wsn for nn calc
		kts = dimsize(wsn,0)
		Make /O /N=(kts) TempWave
		Make /O /N=(kts) DistWave
		newsnName = ReplaceString("_x",xWave,"_dist_sn")
		for (j = 0; j < kts; j+=1)
			for (k = 0; k < kts; k+=1)
				if (k==j)
					TempWave[k]=nan
				else
					TempWave[k]=sqrt(((wsn[j][0]-wsn[k][0])^2)+((wsn[j][1]-wsn[k][1])^2)+((wsn[j][2]-wsn[k][2])^2))
				endif
			endfor
			DistWave[j]=WaveMin(TempWave)
		endfor
		Duplicate /o Distwave, $newsnName
		KillWaves TempWave, DistWave
		matrixop/o c1=col(wsp,0)	//zapnans - wave needs to be broken and then zapped and reassembled
		matrixop/o c2=col(wsp,1)
		matrixop/o c3=col(wsp,2)
		WaveTransform zapnans c1
		WaveTransform zapnans c2
		WaveTransform zapnans c3
		Concatenate /o {c1,c2,c3}, $wspName
		// now use wsp for nn calc
		kts = dimsize(wsp,0)
		Make /O /N=(kts) TempWave
		Make /O /N=(kts) DistWave
		newspName = ReplaceString("_x",xWave,"_dist_sp")
		for (j = 0; j < kts; j+=1)
			for (k = 0; k < kts; k+=1)
				if (k==j)
					TempWave[k]=nan
				else
					TempWave[k]=sqrt(((wsp[j][0]-wsp[k][0])^2)+((wsp[j][1]-wsp[k][1])^2)+((wsp[j][2]-wsp[k][2])^2))
				endif
			endfor
			DistWave[j]=WaveMin(TempWave)
		endfor
		Duplicate /o Distwave, $newspName
		KillWaves TempWave, DistWave
	endfor
End

//// To calculate:
//  Pole to pole distance
//  Distiance from p1 to c and p2 to c
//  Fit a line through p1-c-p2
//  
//  Misaligned chromosomes: distance from p1,c,p2
//  
//  Fluorescence intensities in all channels(?) at each point.

STATIC Function SetupParamWaves()
	Make/O/N=(5)/T ContourDescW = {"Kinetochores Plate","Spindle Poles","BackGround","Misaligned Non-ensheathed","Misaligned Ensheathed"}
	Make/O/N=(5)/T ContourIDW = {"1","2","3","4","5"}
	Make/O/N=(1,3) ScalingW = {0.06449999660253525},{0.06449999660253525},{0.20000000298023224}
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