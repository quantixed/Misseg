#pragma TextEncoding = "UTF-8"
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
	MakeTheLayouts("dist",5,3)
	MakeTheLayouts("spher",5,3)
End

Function WorkflowForXMLAnalysisDir()
	CleanSlate()
	WorkOnDirectory()
	CollectAllMeasurements()
	MakeTheLayouts("dist",5,3)
	MakeTheLayouts("spher",5,3)
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
	
	// determine length c (point c to point p1) 
	Variable cc = sqrt(wx^2 + wy^2 + wz^2)
	
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
		newName = "rn_" + wName // name of normalised rotated pointset
		Duplicate/O M_Product, $newName
		Wave w1 = $newName
		w1 /= cc
	endfor
	
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
		Make/O/N=(nRows,10) $newName
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
		// Now calculate Spherical Coordinates
		// 6 = r
		m1[][6] = sqrt(m0[p][0]^2 + m0[p][1]^2 + m0[p][2]^2)
		// 7 = theta (polar angle)
		m1[][7] = acos(abs(m0[p][2]) / m1[p][6]) // abs will give +ve x,y and 0,pi/2 theta
		// 8 = phi (azimuthal)
		m1[][8] = atan2(abs(m0[p][1]),abs(m0[p][0])) // abs gives 0,pi/2 theta
		// 9 = r (normalized)
		m1[][9] = m1[p][6] / abs(acWave[1])
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
	
	// THE DISTANCE PLOT - distance from point to plane
	String wList = "dist_Mat_1;dist_Mat_4;dist_Mat_5;"
	Variable nWaves = ItemsInList(wList)
	String wName, newName
	Variable nRows
	String plotName = "dist_" + num2str(ii)
	KillWindow/Z $plotName
	Display/N=$plotName/HIDE=1
	
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
		m1[][1] = m0[p][9]
		AppendToGraph/W=$plotName m1[][1] vs m1[][0]
	endfor
	Make/O/N=3 xPos=p
	Make/O/N=3/T xLabel={"Aligned","Misaligned","Misaligned\rEnsheathed"}
	ModifyGraph/W=$plotName userticks(bottom)={xPos,xLabel}
	SetAxis/W=$plotName/A/N=1/E=1 left
	Label/W=$plotName left "Normalized distance to plate"
	ModifyGraph/W=$plotName mode=3,marker=19,mrkThick=0,rgb=(0,0,0,32768)
	SetAxis/W=$plotName bottom -0.5,2.5
	
	// NOW MAKE SPHERICAL COORD GRAPHS
	wList = "dist_Mat_1;dist_Mat_2;dist_Mat_4;dist_Mat_5;"
	nWaves = ItemsInList(wList)
	plotName = "spher_" + num2str(ii)
	String tName // for tracename
	// make colorWave
	Make/O/N=(4,4) cW={{0,0,0,32768},{65535,0,0,32768},{2,39321,1,32768},{0,0,65535,32768}}
	MatrixTranspose cW
	KillWindow/Z $plotName
	Display/N=$plotName/HIDE=1
	
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
		tName = "c7_" + wName
		AppendToGraph/W=$plotName/L=lo m0[][7]/TN=$tName vs m0[][9]
		ModifyGraph/W=$plotName rgb($tName)=(cW[i][0],cW[i][1],cW[i][2],cW[i][3])
		tName = "c8_" + wName
		AppendToGraph/W=$plotName/L=hi m0[][8]/TN=$tName vs m0[][9]
		ModifyGraph/W=$plotName rgb($tName)=(cW[i][0],cW[i][1],cW[i][2],cW[i][3])
	endfor
	ModifyGraph/W=$plotName mode=3,marker=19,mrkThick=0
	ModifyGraph/W=$plotName standoff=0,axisEnab(lo)={0,0.45},axisEnab(hi)={0.55,1},freePos(lo)=0,freePos(hi)=0
	SetAxis/W=$plotName/A/N=1/E=1 bottom
	SetAxis/W=$plotName lo 0,pi/2
	SetAxis/W=$plotName hi 0,pi/2
	Label/W=$plotName bottom "r"
	Label/W=$plotName lo "θ"
	Label/W=$plotName hi "φ"
	Make/O/N=3 yPos={0,pi/4,pi/2}
	Make/O/N=3/T yLabel={"0","π/4","π/2"}
	ModifyGraph/W=$plotName userticks(lo)={yPos,yLabel}
	ModifyGraph/W=$plotName userticks(hi)={yPos,yLabel}
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
	String targetWaveList = "dist_Mat_1;dist_Mat_2;dist_Mat_4;dist_Mat_5;"
	targetWavelist += "rn_Mat_1;rn_Mat_2;rn_Mat_4;rn_Mat_5;"
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
			conName = "all_" + targetName
			Concatenate/O/NP=0 modtList, $conName
		endif
	endfor
End

STATIC Function MakeTheLayouts(prefix,nRow,nCol)
	String prefix
	Variable nRow, nCol
	
	String layoutName = "all"+prefix+"Layout"
	DoWindow/K $layoutName
	NewLayout/N=$layoutName
	String modList = WinList(prefix+"_*",";","WIN:1")
	Variable nWindows = ItemsInList(modList)
	Variable PlotsPerPage = nRow * nCol
	String plotName
	String exString = "Tile/A=(" + num2str(ceil(PlotsPerPage/nCol)) + ","+num2str(nCol)+")"
	
	Variable i,pgNum=1
	
	for(i = 0; i < nWindows; i += 1)
		plotName = StringFromList(i,modList)
		AppendLayoutObject/W=$layoutName/PAGE=(pgnum) graph $plotName
		if(mod((i + 1),PlotsPerPage) == 0 || i == (nWindows -1)) // if page is full or it's the last plot
			LayoutPageAction/W=$layoutName size(-1)=(595, 842), margins(-1)=(18, 18, 18, 18)
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
	String fileName = layoutName + ".pdf"
	SavePICT/PGR=(1,-1)/E=-2/W=(0,0,0,0) as fileName
End

Function MakeTheGizmos()
	SetDataFolder root:
	String wList = WaveList("all_rn_*",";","")
	Variable nWaves = ItemsInList(wList)
	String wName, newName
	Variable nRowsOctant
	Variable i
	
	for(i = 0; i < nWaves; i += 1)
		wName = StringFromList(i,wList)
		newName = ReplaceString("all_rn_",wName,"all_rnOct_") // for the upper octant
		Duplicate/O $wName, $newName
		Wave w0 = $newName
		w0[][] = abs(w0[p][q])
		newName = ReplaceString("all_rn_",wName,"all_rnCube_") // for the whole cube
		nRowsOctant = dimsize(w0,0)
		Make/O/N=(nRowsOctant*8,3) $newName
		Wave w1 = $newName
		w1[][] = w0[mod(p,nRowsOctant)][q]
		MakeCube(w1)
	endfor
	
	String gizList = "dataset;octant;cube;",gizName
	Variable nGiz = ItemsInList(gizList)
	String modList = wList + ReplaceString("all_rn_",wList,"all_rnOct_") + ReplaceString("all_rn_",wList,"all_rnCube_")
	
	for(i = 0; i < nGiz; i += 1)
		gizName = StringFromList(i,gizList)
		KillWindow/Z $gizName
		NewGizmo/N=$gizName
		Wave gW1 = $StringFromList(i*4+0,modList)
		Wave gW2 = $StringFromList(i*4+1,modList)
		Wave gW4 = $StringFromList(i*4+2,modList)
		Wave gW5 = $StringFromList(i*4+3,modList)
		AppendToGizmo/N=$gizName/D Scatter=gW1,name=scatter1
		AppendToGizmo/N=$gizName/D Scatter=gW2,name=scatter2
		AppendToGizmo/N=$gizName/D Scatter=gW4,name=scatter4
		AppendToGizmo/N=$gizName/D Scatter=gW5,name=scatter5
 		ModifyGizmo/N=$gizName ModifyObject=scatter1,objectType=scatter,property={ size,0.2}
 		ModifyGizmo/N=$gizName ModifyObject=scatter2,objectType=scatter,property={ size,0.2}
 		ModifyGizmo/N=$gizName ModifyObject=scatter4,objectType=scatter,property={ size,0.2}
 		ModifyGizmo/N=$gizName ModifyObject=scatter5,objectType=scatter,property={ size,0.2}
		ModifyGizmo/N=$gizName ModifyObject=scatter1,objectType=scatter,property={ color,0,0,0,0.5}
		ModifyGizmo/N=$gizName ModifyObject=scatter2,objectType=scatter,property={ color,1,0,0,1}
		ModifyGizmo/N=$gizName ModifyObject=scatter4,objectType=scatter,property={ color,3.0518e-05,0.6,1.5259e-05,1}
		ModifyGizmo/N=$gizName ModifyObject=scatter5,objectType=scatter,property={ color,0,0,1,1}
	endfor
End


////////////////////////////////////////////////////////////////////////
// Utility functions
////////////////////////////////////////////////////////////////////////

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

STATIC Function MakeCube(m0)
	Wave m0
	Variable nRows = dimSize(m0,0)/8
	// make x
	Make/O/N=(nRows)/FREE ca=1
	Make/O/N=(nRows)/FREE cb=-1
	Concatenate/O/NP=0/FREE {ca,ca,ca,ca,cb,cb,cb,cb},c0
	//make y
	Concatenate/O/NP=0/FREE {ca,cb,ca,cb,ca,cb,ca,cb},c1
	//make z
	Concatenate/O/NP=0/FREE {ca,ca,cb,cb,ca,ca,cb,cb},c2
	//make mat
	Concatenate/O/NP=1/FREE {c0,c1,c2},octMat
	m0[][] = m0[p][q] * octMat[p][q]
End