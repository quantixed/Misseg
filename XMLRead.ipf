#pragma TextEncoding = "MacRoman"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// Front-end dialog to allow user to pick XML file and also to pick image - optional
// Give the option to define which counter point is which
// Give the option to change the scaling of the image
Menu "Macros"
	"Misseg One File...", /Q, WorkflowForXMLAnalysis()
	"Misseg whole directory...", /Q, WorkflowForXMLAnalysisDir()
End

Function WorkflowForXMLAnalysis()
	CleanSlate()
	ReadXML("","")
	DealWithDuplicates()
	ScaleAllWaves()
	RotateAndSitUp()
	DistanceCalculations()
	MakeWavesForGraphAndPlot(0)
	MakeTheLayouts()
End

Function WorkflowForXMLAnalysisDir()
	CleanSlate()
	WorkOnDirectory()
	CollectAllMeasurements()
	MakeTheLayouts()
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

// Deal with duplicate points
STATIC Function DealWithDuplicates()
	WAVE/Z Mat_1, Mat_4, Mat_5
	Variable nPoints0,nPoints1 // in pointset
	Variable counter=0
	Variable i,j
	// Mat_4 has precedence over Mat_1
	nPoints0 = dimsize(Mat_4,0)
	nPoints1 = dimsize(Mat_1,0)
	for(i = 0; i < nPoints0; i += 1)
		for(j = 0; j < nPoints1; j += 1)
			if(Mat_1[j][0] >= Mat_4[i][0] - 2 && Mat_1[j][0] <= Mat_4[i][0] + 2)
				if(Mat_1[j][1] >= Mat_4[i][1] - 2 && Mat_1[j][1] <= Mat_4[i][1] + 2)
					Mat_1[j][] = NaN
					counter += 1
				endif
			endif
		endfor
	endfor
	if(counter > 0)
		print "Deleted", counter, "points from Mat_1 - clash with Mat_4"
	endif
	// Mat_5 has precedence over Mat_1
	nPoints0 = dimsize(Mat_5,0)
	counter = 0
	for(i = 0; i < nPoints0; i += 1)
		for(j = 0; j < nPoints1; j += 1)
			if(Mat_1[j][0] >= Mat_5[i][0] - 2 && Mat_1[j][0] <= Mat_5[i][0] + 2)
				if(Mat_1[j][1] >= Mat_5[i][1] - 2 && Mat_1[j][1] <= Mat_5[i][1] + 2)
					Mat_1[j][] = NaN
					counter += 1
				endif
			endif
		endfor
	endfor
	if(counter > 0)
		print "Deleted", counter, "points from Mat_1 - clash with Mat_5"
	endif
	// Mat_1 clashes with itself
	nPoints0 = dimsize(Mat_1,0)
	counter = 0
	for(i = 0; i < nPoints0; i += 1)
		if(NumType(Mat_1[i][0]) == 2)
			continue
		endif
		for(j = 0; j < nPoints0; j += 1)
			if(i == j)
				continue
			endif
			if(NumType(Mat_1[j][0]) == 2)
				continue
			endif
			if(Mat_1[j][0] >= Mat_1[i][0] - 2 && Mat_1[j][0] <= Mat_1[i][0] + 2)
				if(Mat_1[j][1] >= Mat_1[i][1] - 2 && Mat_1[j][1] <= Mat_1[i][1] + 2)
					if(Mat_1[j][2] >= Mat_1[i][2] - 2 && Mat_1[j][2] <= Mat_1[i][2] + 2)
						Mat_1[j][] = NaN
						counter += 1
					endif
				endif
			endif
		endfor
	endfor
	if(counter > 0)
		print "Deleted", counter, "points from Mat_1 - clash with Mat_1"
	endif
End

// Scale all waves but also delete any without points
STATIC Function ScaleAllWaves()
	String mList = WaveList("mat_*",";","")
	Variable nWaves = ItemsInList(mList)
	String mName
	Wave ScalingW = SetupParamWaves()
	
	Variable nCh = 4
	Variable i
	
	for(i = 0; i < nWaves; i += 1)
		mName = StringFromList(i,mList)
		Wave m0 = $mName
		if(dimsize(m0,0) == 1)
			KillWaves/Z m0
		else
			// Z-channel is Z and T combined so
			m0[][2] = floor((m0[p][q]-1) / nCh)
			// scale wave
			m0[][] *= ScalingW[0][q]
			// get rid of NaNs
			ZapNansFrom2DWave(m0)
		endif
	endfor
End

STATIC Function ZapNansFrom2DWave(m0)
	Wave m0
	Variable nCols = dimsize(m0,1)
	SplitWave/O/N=tmpCol m0
	String wName,originalName
	Variable i
	for(i = 0; i < nCols; i += 1)
		wName = "tmpCol" + num2str(i)
		Wave w0 = $wName
		WaveTransform ZapNans w0
	endfor
	OriginalName = NameOfWave(m0)
	Concatenate/O/KILL S_waveNames, $originalName
end

Function RotateAndSitUp()
	WAVE/Z Mat_2
	if(!WaveExists(Mat_2))
		DoAlert/T="Problem" 0, "No Mat_2 wave"
		return -1
	endif
	// find spindle midpoint, c. Stor in cWave
	Make/O/N=(1,3) cWave
	cWave[0][0] = (Mat_2[0][0] + Mat_2[1][0]) / 2
	cWave[0][1] = (Mat_2[0][1] + Mat_2[1][1]) / 2
	cWave[0][2] = (Mat_2[0][2] + Mat_2[1][2]) / 2
	// find theta and phi for spindle poles
	Variable wx = Mat_2[0][0] - cWave[0][0]
	Variable wy = Mat_2[0][1] - cWave[0][1]
	Variable wz = Mat_2[0][2] - cWave[0][2]
	// inclination/polar, theta. azimuth, phi
	Variable theta = acos(wz / (sqrt( (wx^2) + (wy^2) + (wz^2) ) ) )
	Variable phi = atan2(wy,wx)
	// Make rotation matrix
	Make/O zRotationMatrix={{cos(phi),sin(phi),0},{-sin(phi),cos(phi),0},{0,0,1}}
	Make/O yRotationMatrix={{cos(theta),0,-sin(theta)},{0,1,0},{sin(theta),0,cos(theta)}}
	// Now offset rotate all pointsets
	String wList = WaveList("mat_*",";","")
	Variable nWaves = ItemsInList(wList)
	String wName,newName
	
	Variable i
	
	for (i = 0; i < nWaves; i += 1)
		wName = StringFromList(i,wList)
		newName = "r_" + wName // name of rotated pointset
		Duplicate/O $wName,$newName
		Wave w0 = $newName
		w0[][] -= cWave[0][q] // offset every point
		MatrixMultiply w0, zRotationMatrix // rotate first axis
		Wave M_Product
		MatrixMultiply M_Product, yRotationMatrix // rotate second axis
		Duplicate/O M_Product, $newName
	endfor
	
	// determine length c (point c to point p1) 
	Variable cc = sqrt(wx^2 + wy^2 + wz^2)
	// determine length a (furthest Euclidean distance to point c)
	WAVE/Z r_mat_1
	if(!WaveExists(r_Mat_1))
		DoAlert/T="Problem" 0, "No r_Mat_1 wave"
		return -1
	endif
	Make/O/N=(1,3) oWave={{0},{0},{0}}
	Variable aa = FindMinMaxOfTwo2DWaves(r_Mat_1,oWave,1)
	// store a and c in acWave
	Make/O/N=(2) acWave = {aa,cc}
End

STATIC Function FindMinMaxOfTwo2DWaves(m0,m1,MinOrMax)
	Wave m0,m1
	Variable MinOrMax // min = 0, max = 1
	Variable nRowsA = dimsize(m0,0)
	Variable nRowsB = dimsize(m1,0)
	
	MatrixOp/O/FREE ax = col(m0,0)
	MatrixOp/O/FREE ay = col(m0,1)
	MatrixOp/O/FREE az = col(m0,2)
	MatrixOp/O/FREE bx = col(m1,0)
	MatrixOp/O/FREE by = col(m1,1)
	MatrixOp/O/FREE bz = col(m1,2)
	MatrixOp/O/FREE matAx = colRepeat(ax,nRowsB)
	MatrixOp/O/FREE matAy = colRepeat(ay,nRowsB)
	MatrixOp/O/FREE matAz = colRepeat(az,nRowsB)
	MatrixOp/O/FREE matBx = rowRepeat(bx,nRowsA)
	MatrixOp/O/FREE matBy = rowRepeat(by,nRowsA)
	MatrixOp/O/FREE matBz = rowRepeat(bz,nRowsA)
	MatrixOp/O/FREE distanceX = matAx - matBx
	MatrixOp/O/FREE distanceY = matAy - matBy
	MatrixOp/O/FREE distanceZ = matAz - matBz
	MatrixOp/O/FREE matDist = sqrt(distanceX * distanceX + distanceY * distanceY + distanceZ * distanceZ)
	WaveStats/Q/M=1 matDist
	if(MinOrMax == 1)
		return V_max
	elseif(MinOrMax == 0)
		return V_min
	else
		return -1
	endif
End

Function DistanceCalculations()
	WAVE/Z r_Mat_1,r_Mat_2,r_Mat_4,r_Mat_5
	WAVE/Z cWave,acWave,oWave
	String wList = "r_Mat_1;r_Mat_2;r_Mat_4;r_Mat_5;"
	Variable nWaves = ItemsInList(wList)
	String wName, newName
	Variable nRows
	
	Variable i
	
	for (i = 0; i < nWaves; i += 1)
		wName = StringFromList(i,wList)
		Wave m0 = $wName
		if(!WaveExists(m0))
			continue
		endif
		nRows = DimSize(m0,0)
		if(nRows == 0)
			continue
		endif
		newName = ReplaceString("r_",wName,"dist_")
		Make/O/N=(nRows,6) $newName
		Wave m1 = $newName
		// calculate distances
		// 0 = distance to plane i.e. y = 0
		m1[][0] = m0[p][2]
		// 1 = Euclidean distance to C
		MatrixOp/O/FREE mTemp = sqrt((sumrows(m0 * m0)))
		m1[][1] = mTemp[p]
		// 2 = Euclidean distance to P1
		Duplicate/O/FREE m0, mTemp0
		mTemp0[][] = m0[p][q] - r_Mat_2[0][q]
		MatrixOp/O/FREE mTemp1 = sumrows(mTemp0 * mTemp0)
		MatrixOp/O/FREE mTemp = sqrt(mTemp1)
		m1[][2] = mTemp[p]
		// 3 = Euclidean distance to P2
		mTemp0[][] = m0[p][q] - r_Mat_2[1][q]
		MatrixOp/O/FREE mTemp1 = sumrows(mTemp0 * mTemp0)
		MatrixOp/O/FREE mTemp = sqrt(mTemp1)
		m1[][3] = mTemp[p]
		// 4 = result of f(x) = x^2 over a^2... (spheroid calc)
		m1[][4] = ((m0[p][0]^2 + m0[p][1]^2) / acWave[0]^2) + (m0[p][2]^2 / acWave[1]^2)
		// 5 = normalised distance to plane i.e. distance to y = 0 normalised by c-p distance
		m1[][5] = abs(m0[p][2]) / acWave[1]
	endfor
	ClassifyMisaligned()
End

// Find average normalised distance so it can be collated later
// Classify the misalignment as a proportion and collated later
// perhaps for 1, 4 and 5 we should get mean/median, n, number, outside (spheroid).

STATIC Function Dist3D(w0,row0,w1,row1)
	Wave w0
	Variable row0
	Wave w1
	Variable row1
	MatrixOp/FREE tempmatP = row(w0,row0)
	MatrixOp/FREE tempmatQ = row(w1,row1)
	MatrixOP/FREE tempMat = tempmatP - tempmatQ
	Variable tempVar = norm(tempMat)
	return tempVar
End

Function ClassifyMisaligned()
	WAVE/Z r_Mat_1,r_Mat_2,r_Mat_4,r_Mat_5
	// get pointset of poles plus kientochores at plate
	Concatenate/O/NP=0 {r_Mat_1,r_Mat_2}, r_Mat_Spindle
	Triangulate3D/VOL r_Mat_Spindle
	Variable spindleVol = V_Value
	
	String wList = "r_Mat_4;r_Mat_5;"
	Variable nWaves = ItemsInList(wList)
	String wName, newName
	Variable nRows
	
	Variable i,j
	
	for (i = 0; i < nWaves; i += 1)
		wName = StringFromList(i,wList)
		Wave m0 = $wName
		if(!WaveExists(m0))
			continue
		endif
		nRows = DimSize(m0,0)
		if(nRows == 0)
			continue
		endif
		newName = ReplaceString("r_",wName,"hullCheck_")
		Make/O/N=(nRows) $newName
		Wave m1 = $newName
		for(j = 0; j < nRows; j += 1)
			MatrixOp/O/FREE rowW = row(m0,j)
			Concatenate/O/NP=0 {r_Mat_Spindle,rowW}, spindlePlus
			Triangulate3D/VOL spindlePlus
			m1[j] = spindleVol - V_value
		endfor
	endfor
	KillWaves/Z r_Mat_Spindle,spindlePlus
End

STATIC Function MakeWavesForGraphAndPlot(ii)
	Variable ii // iteration number
	String wList = "dist_Mat_1;dist_Mat_4;dist_Mat_5;"
	Variable nWaves = ItemsInList(wList)
	String wName, newName
	Variable nRows
	String plotName = "plot_" + num2str(ii)
	KillWindow/Z $plotName
	Display/N=$plotName
	
	Variable i
	
	for (i = 0; i < nWaves; i += 1)
		wName = StringFromList(i,wList)
		Wave m0 = $wName
		if(!WaveExists(m0))
			continue
		endif
		nRows = DimSize(m0,0)
		if(nRows == 0)
			continue
		endif
		newName = ReplaceString("dist_",wName,"p_") // p for plot
		Make/O/N=(nRows,2) $newName
		Wave m1 = $newName
		// x value is i + noise
		m1[][0] = i + gnoise(0.1)
		m1[][1] = m0[p][5]
		AppendToGraph/W=$plotName m1[][1] vs m1[][0]
	endfor
	Make/O/N=3 xPos=p
	Make/O/N=3/T xLabel={"Aligned","Misaligned","Misaligned\rEnsheathed"}
	ModifyGraph/W=$plotName userticks(bottom)={xPos,xLabel}
	TidyUpGraph(plotName)
End

////////////////////////////////////////////////////////////////////////
// Work on many files from a directory
////////////////////////////////////////////////////////////////////////


Function WorkOnDirectory()
	NewDataFolder/O/S root:data
	String expDiskFolderName, expDataFolderName
	String FileList, ThisFile
	Variable nWaves, i
	
	NewPath/O/Q/M="Please find disk folder" ExpDiskFolder
	if (V_flag!=0)
		DoAlert 0, "Disk folder error"
		Return -1
	endif
	PathInfo /S ExpDiskFolder
	ExpDiskFolderName = S_path
	FileList = IndexedFile(expDiskFolder,-1,".xml")
	Variable nFiles = ItemsInList(FileList)
	Make/O/N=(nFiles)/T root:fileNameWave
	Wave/T fileNameWave = root:fileNameWave
	
	for(i = 0; i < nFiles; i += 1)
		ThisFile = StringFromList(i,FileList)
		fileNameWave[i] = ThisFile
		expDataFolderName = "dataset_" + num2str(i)
		NewDataFolder/O/S $expDataFolderName
		ReadXML(ExpDiskFolderName,ThisFile)
		DealWithDuplicates()
		ScaleAllWaves()
		RotateAndSitUp()
		DistanceCalculations()
		MakeWavesForGraphAndPlot(i)
		SetDataFolder root:data:
	endfor
	SetDataFolder root:
End

STATIC Function CollectAllMeasurements()
	SetDataFolder root:data:	// relies on earlier load
	DFREF dfr = GetDataFolderDFR()
	String folderName
	Variable numDataFolders = CountObjectsDFR(dfr, 4)
	String wList = ""
	
	Variable i,j
	// assemble a string of semi-colon separated targets in the data folder
	for(i = 0; i < numDataFolders; i += 1)
		folderName = GetIndexedObjNameDFR(dfr, 4, i)
		wList += "root:data:" + folderName + ":thisWave;"
	endfor
	
	// we need to concatenate these waves into root (p_all_*)
	String targetWaveList = "p_Mat_1;p_Mat_4;p_Mat_5;"
	Variable nTargets = ItemsInList(targetWaveList)
	String targetName, tList, conName
	
	SetDataFolder root:
	String fullName,modtList
	
	for(i = 0; i < nTargets; i += 1)
		targetName = StringFromList(i,targetWaveList)
		tList = ReplaceString("thisWave",wList,targetName)
		modtList = tList
		// because some waves might not exist
		for(j = 0; j < numDataFolders; j +=1)
			fullName = StringFromList(j, tList)
			Wave testW = $fullName
			if(!WaveExists(testW))
				modtList = RemoveFromList(fullName,modtList)
			endif
		endfor
		// if there were no waves of that type in any of the data folders
		if(ItemsInList(modtList) > 0)
			conName = ReplaceString("p_",targetName,"p_All_")
			Concatenate/O/NP=0 modtList, $conName
		endif
	endfor
	
//	// redefine targetList for summation of each folder into row of sum_* wave
//	targetWaveList = "Img_VsArea;Img_VsPerimeter;Img_MtPeriTotal;"
//	nTargets = ItemsInList(targetWaveList)
//	String wName,tName
//	for(i = 0; i < nTargets; i += 1)
//		targetName = StringFromList(i,targetWaveList)
//		wName = "root:" + ReplaceString("Img_",targetName,"Sum_")
//		Make/O/N=(numDataFolders) $wName
//		Wave sumW0 = $wName
//		tList = ReplaceString("thisWave",wList,targetName)
//		for(j = 0; j < numDataFolders; j += 1)
//			tName = StringFromList(j,tList)
//			Wave tW0 = $tName
//			sumW0[j] = sum(tW0)
//		endfor
//	endfor
//	
//	// now count the number of Vs and Mt waves in each folder
//	Make/O/N=(numDataFolders) Count_Vs,Count_Mt
//	targetName = "Img_VsArea"
//	tList = ReplaceString("thisWave",wList,targetName)
//	for(i = 0; i < numDataFolders; i += 1)
//		tName = StringFromList(i,tList)
//		Wave tW0 = $tName
//		Count_Vs[i] = numpnts(tW0)
//	endfor
//	targetName = "Img_MtPeriTotal"
//	tList = ReplaceString("thisWave",wList,targetName)
//	for(i = 0; i < numDataFolders; i += 1)
//		tName = StringFromList(i,tList)
//		Wave tW0 = $tName
//		Count_Mt[i] = numpnts(tW0)
//	endfor
//	WAVE/Z Sum_MtPeriTotal,Sum_VsPerimeter
//	MatrixOp/O Ratio_CountVsPerCountMito = Count_Vs / Count_Mt
//	MatrixOp/O Ratio_CountVsPerSumMitoPerim = Count_Vs / Sum_MtPeriTotal
//	MatrixOp/O Ratio_SumVsPerimPerSumMitoPerim = Sum_VsPerimeter / Count_Mt
//	
//	// now we will pull the Img_ClusterSizes into a matrix
//	targetName = "Img_ClusterSizes"
//	tList = ReplaceString("thisWave",wList,targetName)
//	conName = ReplaceString("Img_",targetName,"All_")
//	Concatenate/O/NP=1 tList, tempMat
//	// Img_ClusterSizes lists the number of clusters containing n vesicles, so we'll find the average
//	MatrixOp/O/NTHR=0 $conName = sumRows(tempMat) / numCols(tempMat)
End

STATIC Function MakeTheLayouts()
	DoWindow/K allPlotLayout
	NewLayout/N=allPlotLayout
	String modList = WinList("plot_*",";","WIN:1")
	Variable nWindows = ItemsInList(modList)
	Variable PlotsPerPage = 15 // 5 x 3
	String plotName
	String exString = "Tile/A=(" + num2str(ceil(PlotsPerPage/3)) + ",3)"
	
	Variable i,pgNum=1
	
	for(i = 0; i < nWindows; i += 1)
		plotName = StringFromList(i,modList)
		AppendLayoutObject/W=allPlotLayout/PAGE=(pgnum) graph $plotName
		if(mod((i + 1),PlotsPerPage) == 0 || i == (nWindows -1)) // if page is full or it's the last plot
			LayoutPageAction/W=allPlotLayout size(-1)=(595, 842), margins(-1)=(18, 18, 18, 18)
			ModifyLayout/W=allPlotLayout units=0
			ModifyLayout/W=allPlotLayout frame=0,trans=1
			Execute /Q exString
			if (i != nWindows -1)
				LayoutPageAction/W=allPlotLayout appendpage
				pgNum += 1
				LayoutPageAction/W=allPlotLayout page=(pgNum)
			endif
		endif
	endfor
	SavePICT/PGR=(1,-1)/E=-2/W=(0,0,0,0) as "plots.pdf"
End

////////////////////////////////////////////////////////////////////////
// Utility functions
////////////////////////////////////////////////////////////////////////

STATIC Function TidyUpGraph(plotName)
	String plotName
	SetAxis/A/N=1/E=1 left
	Label left "Normalized distance to plate"
	ModifyGraph mode=3,marker=19,mrkThick=0,rgb=(0,0,0,32768)
	SetAxis bottom -0.5,2.5
End


STATIC Function/WAVE SetupParamWaves()
	Make/O/N=(5)/T ContourDescW = {"Kinetochores Plate","Spindle Poles","BackGround","Misaligned Non-ensheathed","Misaligned Ensheathed"}
	Make/O/N=(5)/T ContourIDW = {"1","2","3","4","5"}
	Make/O/N=3/D ScalingW = {0.06449999660253525,0.06449999660253525,0.20000000298023224}
	MatrixTranspose ScalingW
	return scalingW
End

STATIC Function CleanSlate()
	SetDataFolder root:
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
		if(Stringmatch(name,"*Packages*") != 1)
			KillDataFolder $name		
		endif
	endfor
End