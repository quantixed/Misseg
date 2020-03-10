#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma ModuleName = Misseg

//Notes on Misseg:
//This ipf grew out of XML Analysis.
//1. The original aim was to analyse and plot the spatial distribution of pointsets
//   The pointsets were described in XML files using CellCounter in Fiji.
//   They correspond to the position of kinetochores in hyperstacks.
//   The pointsets are brought into Igor and then displayed.
//   The points are classified in groups in CellCounter.
//   This base workflow is WorkflowForXMLAnalysisDir()
//2. The next iteration was to read the intensity of each channel in a sphere and analyse that:
//   This workflow is called WorkflowForImageAnalysisDir()
//3. After this we added the ability to clip images centred on each point and do all kinds of stuff
//   For eample rotating and averaging. 2 & 3 are not currently used but are kept in this ipf.
//4. A version of Misseg where the exclusion zone is used to rescale all the coordinates
//   So that their position reflects their depth into the ER that surrounds the spindle.
//   This version requires segmented images to be fed in.
//   A "hybrid" version was created that does this analysis and then uses the original method
//   to display. Shows a box plot of the exclusion zone ratio.


////////////////////////////////////////////////////////////////////////
// Menu items
////////////////////////////////////////////////////////////////////////
Menu "Misseg"
	Submenu "Misseg functions"
		"Misseg (Spatial)...", /Q, WorkflowForXMLAnalysisDir()
		"Misseg (Exclusion Zone)...", /Q, WorkflowForXMLAnalysisWithER()
		"Image Analysis (Spatial)...", /Q, WorkflowForImageAnalysisDir()
		"Image Clipping Extension...", /Q, WorkflowForImageClipping()
	End
End

////////////////////////////////////////////////////////////////////////
// Master functions and wrappers
////////////////////////////////////////////////////////////////////////

Function WorkflowForXMLAnalysisDir()
	CleanSlate()
	StartingPanelforXML()
End

Function WorkflowForImageAnalysisDir()
	CleanSlate()
	StartingPanelForIA()
End

Function WorkflowForImageClipping()
	CleanSlate()
	StartingPanelForIC()
End

Function WorkflowForXMLAnalysisWithER()
	CleanSlate()
	StartingPanelforER()
End

// After clicking do it on the panel we come here
STATIC Function XMLWrapperFunc()
	KillWindow/Z SetUp
	WorkOnDirectory()
	CollectAllMeasurements()
	MakeTheLayouts("dist",5,3)
	MakeTheLayouts("spher",5,3)
	MakeTheGizmos()
End

// After clicking do it on the panel we come here
STATIC Function IAWrapperFunc()
	KillWindow/Z SetUp
	WorkOnDirectoryIA()
	MakeTheLayouts("intens",5,3)
	SummariseIntensityMeasurements(4)
	MakeTheLayouts("mean",6,2)
End

// After clicking do it on the panel we come here
STATIC Function ICWrapperFunc()
	KillWindow/Z SetUp
	WorkOnDirectoryIC()
	WAVE/Z gVarWave = root:gVarWave
	// size the images on the layout to clipsize
	MakeTheLayouts("p_rgb",floor((gVarWave[1]/2) * 1.5),ceil(gVarWave[1]/2),filtVar = 1)
	MakeTheLayouts("p_rgb",floor((gVarWave[1]/2) * 1.5),ceil(gVarWave[1]/2),filtVar = 2)
	MakeTheLayouts("p_rgb",floor((gVarWave[1]/2) * 1.5),ceil(gVarWave[1]/2),filtVar = 4)
	MakeTheLayouts("p_rgb",floor((gVarWave[1]/2) * 1.5),ceil(gVarWave[1]/2),filtVar = 5)
	CollateImagesToAverageByDistance(20,10)
	MakeTheLayouts("p_rot",floor((gVarWave[1]/2) * 1.5),ceil(gVarWave[1]/2),filtVar = 1)
	MakeTheLayouts("p_rot",floor((gVarWave[1]/2) * 1.5),ceil(gVarWave[1]/2),filtVar = 2)
	MakeTheLayouts("p_rot",floor((gVarWave[1]/2) * 1.5),ceil(gVarWave[1]/2),filtVar = 4)
	MakeTheLayouts("p_rot",floor((gVarWave[1]/2) * 1.5),ceil(gVarWave[1]/2),filtVar = 5)
	KillTheLayout("p_rgb",1)
	KillTheLayout("p_rot",1)
End

// After clicking Do It (or Hybrid) on the panel we come here
STATIC Function ERWrapperFunc(VARIABLE optVar)
	KillWindow/Z SetUp
	WorkOnDirectoryER(optVar)
	CollectAllMeasurements()
	MakeTheGizmos()
	PlotRatios()
End

////////////////////////////////////////////////////////////////////////
// Main functions
////////////////////////////////////////////////////////////////////////
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
	// store the name of original image in a textwave called fNameWave
	XPathStr0 = "/CellCounter_Marker_File/Image_Properties/Image_Filename"
	XMLwaveFmXpath(fileID,XPathStr0,"","")
	Make/O/N=(1)/T fNameWave
	// cannot assume first row of M_XMLcontent holds filename
	fNameWave = M_XMLcontent[DimSize(M_XMLcontent,0)-1]
	// finally close the XML file
	xmlclosefile(fileID,0)
end

// Deal with duplicate points
STATIC Function DealWithDuplicates()
	WAVE/Z Mat_1, Mat_2, Mat_4, Mat_5
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
		print "Deleted", counter, "points from Mat_1 - clash with Mat_4 (", GetDataFolder(0) ,")"
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
		print "Deleted", counter, "points from Mat_1 - clash with Mat_5 (", GetDataFolder(0) ,")"
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
		print "Deleted", counter, "points from Mat_1 - clash with Mat_1 (", GetDataFolder(0) ,")"
	endif
	// too many points in Mat_2
	nPoints0 = dimsize(Mat_2,0)
	if(nPoints0 > 2)
		// assume that 1st two rows are genuine (the code works on this basis)
		DeletePoints 2,(nPoints0 - 2), Mat_2
		print "Deleted", nPoints0 - 2, "points from Mat_1 - too many points (", GetDataFolder(0) ,")"
	endif
End

// Scale all waves but also delete any without points
STATIC Function ScaleAllWaves(correctVar,nCh)
	Variable correctVar // bit 0 = correct Z, bit 1 = scale the values
	Variable nCh
	correctVar = trunc(correctVar)
	String mList = WaveList("mat_*",";","")
	Variable nWaves = ItemsInList(mList)
	String mName
	Wave/Z ScalingW = root:ScalingW
	if(!WaveExists(ScalingW))
		Print "No scaling wave available."
		return -1
	endif
	
	Variable i
	
	for(i = 0; i < nWaves; i += 1)
		mName = StringFromList(i,mList)
		Wave m0 = $mName
		if(dimsize(m0,0) == 1)
			if(m0[0][2] == -1 || sum(m0) == 0)
				KillWaves/Z m0
			endif
		else
			if((correctVar & 2^0) != 0)
				// Z-channel is Z and C combined so
				m0[][2] = ceil((m0[p][q]) / nCh) - 1 // -1 to make 0-based
			endif
			if((correctVar & 2^1) != 0)
				// scale wave
				m0[][] *= ScalingW[0][q]
			endif
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

// Find the angle for rotation of clip relative to closest spindle pole
// We will make a new wave named ont_ (for orientation) that contains the
// angle for rotation, the closest pole and the distance to it 
STATIC Function OrientationOfObjects()
	String mList = WaveList("mat_*",";","")
	Variable nWaves = ItemsInList(mList)
	String mName, newName
	Variable nRows,closestPole
	WAVE/Z mat_2
	Variable dist0,dist1
	
	Variable i,j
	
	for(i = 0; i < nWaves; i += 1)
		mName = StringFromList(i,mList)
		Wave m0 = $mName
		nRows = DimSize(m0,0)
		newName = ReplaceString("mat_",mName,"ont_")
		Make/O/N=(nRows,3) $newName
		Wave w1 = $newName
		for(j = 0; j < nRows; j += 1)
			// find the closest pole
			dist0 = sqrt((m0[j][0] - mat_2[0][0])^2 + (m0[j][1] - mat_2[0][1])^2)
			dist1 = sqrt((m0[j][0] - mat_2[1][0])^2 + (m0[j][1] - mat_2[1][1])^2)
			// store shortest distance in column 0
			w1[j][0] = min(dist0,dist1)
			if(dist0 > dist1)
				closestPole = 1
			else
				closestPole = 0
			endif
			// store clostest pole
			w1[j][1] = closestPole
			// work out the angle to the nearest pole for each row of the object and store that
			w1[j][2] = 90 - (atan2(mat_2[closestPole][1] - m0[j][1],mat_2[closestPole][0] - m0[j][0]) * (180/Pi))
			// special case of mat_2
			if(stringmatch(mName, "mat_2") == 1)
				if(j == 0)
					closestPole = 1
				elseif(j == 1)
					closestPole = 0
				endif
				w1[j][2] = 90 - (atan2(mat_2[closestPole][1] - m0[j][1],mat_2[closestPole][0] - m0[j][0]) * (180/Pi))
			endif
		endfor
	endfor
End

Function RotateAndSitUp(CorrectVar)
	Variable CorrectVar
	WAVE/Z Mat_2
	if(!WaveExists(Mat_2))
		DoAlert/T="Problem" 0, "No Mat_2 wave"
		return -1
	endif
	WAVE/Z cWave, scalingW = root:scalingW
	// at this point we will have a cWave if FindCentre has been done (Exclusion Zone workflow)
	// that cWave (if it exists) will need scaling.
	// if it doesn't exist we need to make it, based on the standard method
	// if we are using the hybrid approach we need to deal with that scenario
	if(CorrectVar == 0)
		KillWaves/Z cWave
	endif
	if(!WaveExists(cWave))
		// find spindle midpoint, c. Store in cWave
		Make/O/N=(1,3) cWave
		cWave[0][0] = (Mat_2[0][0] + Mat_2[1][0]) / 2
		cWave[0][1] = (Mat_2[0][1] + Mat_2[1][1]) / 2
		cWave[0][2] = (Mat_2[0][2] + Mat_2[1][2]) / 2
	else
		// cWave needs scaling
		cWave[][] *= ScalingW[0][q]
	endif
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
		// Make normalized set for all but mat_3
		if(stringmatch(wName,"mat_3") == 0)
			newName = "rn_" + wName // name of normalised rotated pointset
			Duplicate/O M_Product, $newName
			Wave w1 = $newName
			if((correctVar & 2^0) != 0) // if bit 0 is set, norm to ER and not poles
				Wave/Z r1 = $(ReplaceString("mat_",wName,"rat_"))
				if(!WaveExists(r1))
					continue
				endif
				w1[][] /= sqrt(w1[p][0]^2 + w1[p][1]^2 + w1[p][2]^2) / r1[p] // very inefficient
			else
				w1 /= cc
			endif
		endif
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

Function FindTheCentre()
	WAVE/Z Mat_1,Mat_2
	if(!WaveExists(Mat_1))
		DoAlert/T="Problem" 0, "No Mat_1 wave"
		return -1
	endif
	// find centroid of plate kinetochores, c. Store in cWave
	Make/O/N=(1,3) cWave
	WaveStats/Q/M=1/RMD=[][0] Mat_1
	cWave[0][0] = V_Avg
	WaveStats/Q/M=1/RMD=[][1] Mat_1
	cWave[0][1] = V_Avg
	WaveStats/Q/M=1/RMD=[][2] Mat_1
	cWave[0][2] = V_Avg
	// determine length c (point c to point p1) 
	Variable wx = Mat_2[0][0] - cWave[0][0]
	Variable wy = Mat_2[0][1] - cWave[0][1]
	Variable wz = Mat_2[0][2] - cWave[0][2]
	Variable cc = sqrt(wx^2 + wy^2 + wz^2)
	return cc
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
	// get pointset of poles plus kinetochores at plate
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

STATIC Function LoadImageAndAnalyse(ImageDiskFolderName,ImageFileName,nCh,rr)
	String ImageDiskFolderName		// Name of Igor symbolic path or "" to get dialog
	String ImageFileName		// Name of file to load or "" to get dialog
	Variable nCh,rr // number of channels and radius size in pixels from dialog
	// now convert .dv to .tif
	if(StringMatch(ImageFileName,"*.tif") == 0)
		ImageFileName = RemoveEnding(ImageFileName,".dv") + ".tif"
	endif
	
	NewPath/O/Q ImageDiskFolderPath, ImageDiskFolderName
	ImageLoad/T=tiff/N=OriginalImage/O/P=ImageDiskFolderPath/S=0/LR3D/C=-1 ImageFileName
	Wave/Z OriginalImage
	if(!WaveExists(OriginalImage))
		DoAlert/T="Problem" 0, "Image did not load"
		return -1
	endif
	// Rearrange the image. Upon load, channels are arrayed consecutively
	// will be x y z c
	if(RejigImage(OriginalImage,nCh,1) != 1)
		DoAlert/T="Problem" 0, "Image rearrangement did not execute."
		return -1
	else
		Wave ImgMat
	endif
	// declare the waves from a previous load
	WAVE/Z mat_1,mat_2,mat_3,mat_4,mat_5
	if(!WaveExists(mat_3))
		DoAlert/T="Problem" 0, "No mat_3 (background counter locations)"
		return -1
	else
		WAVE bgW = MakeBgWave(imgMat,mat_3,rr)
	endif
	// subtract background values from signed 16-bit ImgMat
	ImgMat[][][][] -= bgW[0][s]
	MeasureIntensities(mat_1,imgMat,rr)
	MeasureIntensities(mat_2,imgMat,rr)
	MeasureIntensities(mat_4,imgMat,rr)
	MeasureIntensities(mat_5,imgMat,rr)
	TidyUp()
	return 1
End

STATIC Function LoadImageAndClip(ImageDiskFolderName,ImageFileName,nCh,clipSize,ii)
	String ImageDiskFolderName		// Name of Igor symbolic path or "" to get dialog
	String ImageFileName		// Name of file to load or "" to get dialog
	Variable nCh,clipSize // number of channels and size of clip in pixels from dialog
	Variable ii // index of which dataset we're working on
	// now convert .dv to .tif
	if(StringMatch(ImageFileName,"*.tif") == 0)
		ImageFileName = RemoveEnding(ImageFileName,".dv") + ".tif"
	endif
	
	NewPath/O/Q ImageDiskFolderPath, ImageDiskFolderName
	ImageLoad/T=tiff/N=OriginalImage/O/P=ImageDiskFolderPath/S=0/LR3D/C=-1 ImageFileName
	Wave/Z OriginalImage
	if(!WaveExists(OriginalImage))
		DoAlert/T="Problem" 0, "Image did not load"
		return -1
	endif
	// Rearrange the image. Upon load, channels are arrayed consecutively
	// will be x y z c
	if(RejigImage(OriginalImage,nCh,0) != 1)
		DoAlert/T="Problem" 0, "Image rearrangement did not execute."
		return -1
	else
		Wave imgMat16
	endif
	// We need to reference gVarWave to get rCh,gCh,bCh
	Wave gVarWave = root:gVarWave
	// Now find the min and max values (trimmed)
	MakeMinMaxWave(imgMat16,gVarWave[5],gVarWave[6],gVarWave[7])
	// declare the waves from a previous load
	WAVE/Z mat_1,mat_2,mat_4,mat_5
	TakeClips(mat_1,imgMat16,clipSize,ii)
	TakeClips(mat_2,imgMat16,clipSize,ii)
	TakeClips(mat_4,imgMat16,clipSize,ii)
	TakeClips(mat_5,imgMat16,clipSize,ii)
	TidyUp()
	return 1
End

STATIC Function LoadThresholdedImage(ImageDiskFolderName,ImageFileName)
	String ImageDiskFolderName		// Name of Igor symbolic path or "" to get dialog
	String ImageFileName		// Name of file to load or "" to get dialog
	// now convert .dv to .tif
	if(StringMatch(ImageFileName,"*.tif") == 0)
		ImageFileName = RemoveEnding(ImageFileName,".dv") + ".tif"
	endif
	
	NewPath/O/Q ImageDiskFolderPath, ImageDiskFolderName
	ImageLoad/T=tiff/N=OriginalImage/O/P=ImageDiskFolderPath/S=0/LR3D/C=-1 ImageFileName
	Wave/Z OriginalImage
	if(!WaveExists(OriginalImage))
		DoAlert/T="Problem" 0, "Image did not load"
		return -1
	endif
	// we need to find the centre ~ this is also done later in RotateAndSitUp
	Variable cc = FindTheCentre() // this is the half spindle length, in pixels.
	// list of waves to process
	String wList = "mat_1;mat_2;mat_4;mat_5;"
	String mName
	Variable nWaves = ItemsInList(wList)
	Variable i
	for(i = 0; i < nWaves; i += 1)
		mName = StringfromList(i, wList)
		Wave/Z m0 = $mName
		if(WaveExists(m0) == 1)
			Resample3DLine(m0,cc)
			InterpolateAndFindThreshold(m0,OriginalImage)
		endif
	endfor
	return 1
End

STATIC Function RejigImage(ImageMat,nCh,USVar)
	WAVE ImageMat
	Variable nCh
	Variable USVar // unsigned is 0, signed is 1
	Variable nFrames = DimSize(ImageMat,2)
	if(mod(nFrames,nCh) > 0)
		return -1
	endif
	// 16-bit signed or unsigned
	if(USVar == 0)
		Make/O/N=(dimSize(imageMat,0),dimSize(imageMat,1),nFrames/nCh,nCh)/U/W ImgMat16
		ImgMat16[][][][] = ImageMat[p][q][(r * nCh) + mod(s,nCh)]
	elseif(USVar == 1)
		Make/O/N=(dimSize(imageMat,0),dimSize(imageMat,1),nFrames/nCh,nCh)/W ImgMat
		ImgMat[][][][] = ImageMat[p][q][(r * nCh) + mod(s,nCh)]
	endif
	return 1
End

Function/WAVE MakeBgWave(ImageMat,bgMat,rr)
	Wave ImageMat,bgMat
	Variable rr

	Variable nPoints = DimSize(bgmat,0)
	Variable nChannels = DimSize(ImageMat,3)
	Make/O/N=(nPoints,nChannels)/FREE bgValueW
	Variable xPos,yPos,zPos,xPosMin,xPosMax,yPosMin,yPosMax
	
	Variable i,j
	
	for(i = 0; i < nPoints; i += 1)
		xPos = round(bgMat[i][0])
		yPos = round(bgMat[i][1])
		zPos = round(bgMat[i][2])
		// first check if bg area bumps into edges on x
		if(xPos - rr >= 0)
			xPosMin = xPos - rr
			xPosMax = xPos + rr
			if(xPos + rr >= DimSize(ImageMat,0))
				xPosMax = DimSize(ImageMat,0) - 1
				xPosMin = xPosMax - rr * 2
			endif
		else
			xPosMin = 0
			xPosMax = rr * 2
		endif
		// now check if bg area bumps into edges on y
		if(yPos - rr >= 0)
			yPosMin = yPos - rr
			yPosMax = yPos + rr
			if(yPos + 2 > DimSize(ImageMat,1))
				yPosMax = DimSize(ImageMat,1) - 1
				yPosMin = yPosMax - rr* 2
			endif
		else
			yPosMin = 0
			yPosMax = rr * 2
		endif
		// now measure the bg values, in one z-slice only
		for(j = 0; j < nChannels; j += 1)
			if(rr == 0)
				bgValueW[i][j] = ImageMat[xPos][yPos][zPos][j]
			else
				WaveStats/Q/M=1/RMD=[xPosMin,xPosMax][yPosMin,yPosMax][zPos][j] ImageMat
				bgValueW[i][j] = V_Avg
			endif
		endfor
	endfor
	// 1 row wave with average bg value for each channel in each column
	MatrixOp/O bgW = averageCols(bgValueW)
	return bgW
End

Function MeasureIntensities(objToMeasure,ImageMat,rr)
	Wave/Z objToMeasure,ImageMat
	Variable rr // radius in pixels
	if(!WaveExists(objToMeasure))
		return -1
	endif
	Variable nPoints = DimSize(objToMeasure,0)
	if(nPoints == 0)
		return -1
	endif
	Variable nChannels = dimsize(ImageMat,3)
	// the objects to measure are stored in mat_n, store intensities in int_n
	String wName = ReplaceString("mat_",NameOfWave(objToMeasure),"int_")
	Make/O/N=(nPoints,nChannels) $wName
	Wave w0 = $wName
	// for each row in the objects to measure we'll excise a hypercube of image data
	// Image is arrange x y z c
	Variable xPos,yPos,zPos
	Variable xPosMin,xPosMax,yPosMin,yPosMax,zPosMin,zPosMax
	Wave/z ScalingW = root:ScalingW
	Variable zr = ceil(rr * (ScalingW[0][0] / ScalingW[0][2])) // rr in terms of z scaling
	
	Variable i,j
	
	for(i = 0; i < nPoints; i += 1)
		// integers stored in XML are 0-based for X and Y. Z has been converted
		xPos = round(objToMeasure[i][0])
		yPos = round(objToMeasure[i][1])
		zPos = round(objToMeasure[i][2])
		// first check if cube area bumps into edges on x
		if(xPos - rr >= 0)
			xPosMin = xPos - rr
			xPosMax = xPos + rr
			if(xPos + rr >= DimSize(ImageMat,0))
				xPosMax = DimSize(ImageMat,0) - 1
				xPosMin = xPosMax - (rr * 2)
			endif
		else
			xPosMin = 0
			xPosMax = rr * 2
		endif
		// now check if cube area bumps into edges on y
		if(yPos - rr >= 0)
			yPosMin = yPos - rr
			yPosMax = yPos + rr
			if(yPos + rr > DimSize(ImageMat,1))
				yPosMax = DimSize(ImageMat,1) - 1
				yPosMin = yPosMax - (rr * 2)
			endif
		else
			yPosMin = 0
			yPosMax = rr * 2
		endif
		// now check if cube area bumps into edges on z
		if(zPos - zr >= 0)
			zPosMin = zPos - zr
			zPosMax = zPos + zr
			if(zPos + zr > DimSize(ImageMat,2))
				zPosMax = DimSize(ImageMat,2) - 1
				zPosMin = zPosMax - (zr * 2)
			endif
		else
			zPosMin = 0
			zPosMax = zr * 2
			// for large rr it's possible we run out of z-slices, warn about this.
			if(zPosMax >= DimSize(ImageMat,2))
				Print "zDim exceeded for", GetDataFolder(0)
			endif
		endif
		// save some computing time is rr = 0
		if(rr == 0)
			w0[i][] = ImageMat[xPos][yPos][zPos][q]
		else
			// excise the hypercube -> cubeImg
			Duplicate/O/FREE/RMD=[xPosMin,xPosMax][yPosMin,yPosMax][zPosMin,zPosMax][] ImageMat, cubeImg
			// mask according to spherical ROI and measure mean voxel density for each channel
			Wave theReturnWave = SendToSphere(cubeImg)
			// mean voxel density is stored in each column, store these in int_n for each row of Objects to Measure.
			w0[i][] = theReturnWave[0][q]
			KillWaves/Z theReturnWave
		endif
	endfor
End

Function/WAVE SendToSphere(cubeImg)
	Wave cubeImg
	// deduce variable rr from the size of cubeImg
	Variable rr = (DimSize(cubeImg,0) - 1) / 2
	// Because of scaling we will use this instead to make sphereMask
	Variable xmid = (DimSize(cubeImg,0) - 1) / 2
	Variable ymid = (DimSize(cubeImg,1) - 1) / 2
	Variable zmid = (DimSize(cubeImg,2) - 1) / 2
	Wave/Z sphereMask = root:sphereMask
	if(!WaveExists(sphereMask))
		MakeSphereMask(xMid,yMid,zMid)
	endif
	Wave/Z sphereMask = root:sphereMask
	// multiply cubeImg by mask
	cubeImg *= sphereMask
	Variable nChannels = DimSize(cubeImg,3)
	// how many voxels for the mean intensity?
	Variable roiSizeInVoxels = sum(sphereMask) / nChannels
	Make/O/N=(1,nChannels) theROIAverageWave
	
	Variable i
	
	for(i = 0; i < nChannels; i += 1)
		WaveStats/Q/M=1/RMD=[][][][i] cubeImg
		theROIAverageWave[][i] = V_sum / roiSizeInVoxels
	endfor
	return theROIAverageWave
End

STATIC Function MakeSphereMask(xMid,yMid,zMid)
	Variable xMid,yMid,zMid
	Variable rr = xMid
	Wave gVarWave = root:gVarWave
	Variable nCh = gVarWave[0]
	Make/O/N=(xMid * 2 + 1,yMid * 2 + 1,zMid * 2 + 1,nCh)/B root:sphereMask = 0
	Wave sphereMask = root:sphereMask
	if(rr == 0)
		sphereMask = 1
		return 1
	endif
	Wave/Z ScalingW = root:ScalingW
	if(!WaveExists(ScalingW))
		DoAlert/T="Problem" 0, "No scaling wave"
		return -1
	endif
	// centre of cube is rr,rr,rr
	// scalingW has pixel dimensions in columns
	Variable rs = rr * ScalingW[0][0]
	sphereMask[][][][] = (sqrt( (ScalingW[0][0] * (p-xmid))^2 + (ScalingW[0][1] * (q-ymid))^2 + (ScalingW[0][2] * (r-zmid))^2) <= rs) ? 1 : 0
	return 1
End

Function TidyUp()
	WAVE/Z M_xmlcontent,OriginalImage,ScalingW,W_ElementList,W_xmlcontentnodes,ImgMat,ImgMat16
	KillWaves/Z M_xmlcontent,OriginalImage,ScalingW,W_ElementList,W_xmlcontentnodes,ImgMat,ImgMat16
End

STATIC Function MakeIntWavesForGraphAndPlot(ii)
	Variable ii // iteration number
	
	// THE INTENSITY PLOT - distance from point to plane
	String wList = "int_1;int_2;int_4;int_5;"
	Variable nWaves = ItemsInList(wList)
	String wName, newName
	Variable nRows, nChannels
	String plotName = "intens_" + num2str(ii)
	KillWindow/Z $plotName
	Display/N=$plotName/HIDE=1
	Wave/T objWave = root:objWave
	Wave/T ChWave = root:ChWave
	
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
		nChannels = DimSize(m0,1)
		for(j = 0; j < nChannels; j += 1)
			newName = ReplaceString("int_",wName,"p_") + "_" + num2str(j) // p for plot
			Make/O/N=(nRows,2) $newName
			Wave m1 = $newName
			// x value distribute by loc group by channel
			m1[][0] = (i) + (j * nChannels) + gnoise(0.1)
			m1[][1] = m0[p][j]
			AppendToGraph/W=$plotName m1[][1] vs m1[][0]
		endfor
	endfor
	Make/O/N=16 xPos = p
	Make/O/N=16/T xLabel
	xLabel[0,;4]=ObjWave[0]
	xLabel[1,;4]=ObjWave[1]
	xLabel[2,;4]=ObjWave[3]
	xLabel[3,;4]=ObjWave[4]
	ModifyGraph/W=$plotName userticks(bottom)={xPos,xLabel}
	SetAxis/W=$plotName/A/N=1/E=1 left
	Label/W=$plotName left "Fluorescence Intensity"
	ModifyGraph/W=$plotName mode=3,marker=19,mrkThick=0,rgb=(0,0,0,32768)
	SetAxis/W=$plotName bottom -0.5,15.5
	ModifyGraph/W=$plotName tkLblRot(bottom)=90
	// add channel labels
	for(i = 0; i < nChannels; i += 1)
		SetDrawEnv/W=$plotName xcoord= bottom,textxjust= 1,textyjust= 2, fsize= 10
		DrawText/W=$plotName ((nChannels - 1) / 2 + (nChannels * i)),0,ChWave[i]
	endfor
End

STATIC Function MakeMinMaxWave(ImageMat,rCh,gCh,bCh)
	Wave ImageMat
	Variable rCh,gCh,bCh // index of chunk for each channel. -1 means blank
	String rgbList = num2str(rCh) + ";" + num2str(gCh) + ";" + num2str(bCh) + ";"
	Variable nZ = DimSize(ImageMat,2)
	Make/O/N=(nZ,2,3) minMaxW // rows for z, columns min and max, layers r g b
	Variable chSelect
	Variable minLevel, maxLevel, range
	Variable lowerCentile = 0.0035 / 2
	Variable upperCentile = 1 - lowerCentile
	
	Variable i,j
	
	for(i = 0; i < nZ; i += 1)
		for(j = 0; j < 3; j += 1)
			chSelect = str2num(StringFromList(j,rgbList))
			if(chSelect == -1)
				continue
			endif
			Duplicate/O/FREE/RMD=[][][i][chSelect] ImageMat, tempW
			Redimension/N=(-1,-1) tempW
			ImageHistogram/I tempW
			WAVE/Z W_ImageHist
			Integrate W_ImageHist/D=cumHist
			minLevel = cumHist[0]
			maxLevel = cumHist[dimsize(cumHist,0) - 1]
			range = maxLevel - minLevel
			FindLevel/Q cumHist, minLevel + range * lowerCentile
			minMaxW[i][0][j] = V_levelX
			FindLevel/Q cumHist, minLevel + range * upperCentile
			minMaxW[i][1][j] = V_levelX
			//MatrixOp/O thresholdedImage = clip(Image,lowerThreshold,upperThreshold)
			KillWaves/Z W_ImageHist, cumHist
		endfor
	endfor
End

Function TakeClips(objToMeasure,ImageMat,clipSize,ii)
	Wave/Z objToMeasure,ImageMat
	Variable clipSize // clipSize in pixels
	Variable ii // index to tell us what dataset we're working on
	if(!WaveExists(objToMeasure))
		return -1
	endif
	Variable nPoints = DimSize(objToMeasure,0)
	if(nPoints == 0)
		return -1
	endif
	// We need to reference gVarWave to get rCh,gCh,bCh
	Wave gVarWave = root:gVarWave
	Variable nChannels = dimsize(ImageMat,3)
	// for each row in the objects to measure we'll excise a layer of image data
	// Image is arrange x y z c
	Variable xPos,yPos,zPos
	Variable xPosMin,xPosMax,yPosMin,yPosMax
	Variable rr = (clipSize - 1) / 2
	String clipName,plotName
	
	Variable i
	
	for(i = 0; i < nPoints; i += 1)
		// integers stored in XML are 0-based for X and Y. Z has been converted
		xPos = round(objToMeasure[i][0])
		yPos = round(objToMeasure[i][1])
		zPos = round(objToMeasure[i][2])
		// first check if clip area bumps into edges on x, pass if so
		xPosMin = xPos - rr
		xPosMax = xPos + rr
		if(xPosMin < 0 || xPosMax >= DimSize(ImageMat,0))
			continue
		endif
		// now check if cube area bumps into edges on y
		yPosMin = yPos - rr
		yPosMax = yPos + rr
		if(yPosMin < 0 || yPosMax >= DimSize(ImageMat,1))
			continue
		endif
		if(zPos < 0 || zPos >= DimSize(ImageMat,2))
			continue
		endif
		// excise the clip
		clipName = ReplaceString("mat_",NameOfWave(objToMeasure),"clp_") + "_" + num2str(i)
		Duplicate/O/RMD=[xPosMin,xPosMax][yPosMin,yPosMax][zPos][] ImageMat, $clipName
		Wave clip0 = $clipName
		// Make the appropriate RGB images here and make a hidden window - we'll put them in a layout later
		Wave clip1 = MakeTheImage(clip0, zPos,gVarWave[5],gVarWave[6],gVarWave[7])
		KillWaves/Z clip0
		// plotName is p_rgb_0_1_2 for dataset 0, obj 1, row 2
		plotName = ReplaceString("clp_",clipName,"p_rgb_"+num2str(ii)+"_")
		KillWindow/Z $plotName
		NewImage/N=$plotName/HIDE=1/S=0 clip1
		ModifyGraph/W=$plotName width={Plan,1,top,left}
		// Now we'll make a rotated version
		Wave clip2 = MakeRotatedImage(clip1)
		// plotName is p_rgb_0_1_2 for dataset 0, obj 1, row 2
		plotName = ReplaceString("p_rgb_",plotName,"p_rot_")
		KillWindow/Z $plotName
		NewImage/N=$plotName/HIDE=1/S=0 clip2
		ModifyGraph/W=$plotName width={Plan,1,top,left}
	endfor
End

STATIC Function/WAVE MakeTheImage(clip0,zPos,rCh,gCh,bCh)
	Wave clip0
	Variable zPos
	Variable rCh,gCh,bCh // index of chunk for each channel. -1 means blank
	
	WAVE/Z minMaxW
	if(!WaveExists(minMaxW))
		DoAlert/T="Problem" 0, "Missing the min and max wave"
	endif
	String clipName = NameOfWave(clip0)
	Variable nCh = DimSize(clip0,3)
	Variable xSize = DimSize(clip0,0)
	Variable ySize = DimSize(clip0,1)	// should be the same as xSize
	Variable totalPx = xSize * ySize

	String newName = ReplaceString("clp_",clipName,"clp_rgb_")
	Make/O/N=(xSize,ySize,3)/W/U $newName
	Wave rgbClip = $newName
	// for testing
//	rgbClip[][][0] = clip0[p][q][0][rCh]
//	rgbClip[][][1] = clip0[p][q][0][gCh]
//	rgbClip[][][2] = clip0[p][q][0][bCh]
	if(rCh == -1)
		rgbClip[][][0] = 0
	else
		rgbClip[][][0] = limit(65535 * (clip0[p][q][0][rCh] - minMaxW[zPos][0][0]) / (minMaxW[zPos][1][0] - minMaxW[zPos][0][0]),0,65535)
	endif
	if(gCh == -1)
		rgbClip[][][1] = 0
	else
		rgbClip[][][1] = limit(65535 * (clip0[p][q][0][gCh] - minMaxW[zPos][0][1]) / (minMaxW[zPos][1][1] - minMaxW[zPos][0][1]),0,65535)
	endif
	if(bCh == -1)
		rgbClip[][][2] = 0
	else
		rgbClip[][][2] = limit(65535 * (clip0[p][q][0][bCh] - minMaxW[zPos][0][2]) / (minMaxW[zPos][1][2] - minMaxW[zPos][0][2]),0,65535)
	endif
	return rgbClip
End

STATIC Function/WAVE MakeRotatedImage(clip1)
	Wave clip1
	String clipName = NameOfWave(clip1)
	// clipname is this so let's parse it clip_rgb_1_2
	String expr = "clp\\wrgb\\w([[:digit:]]+)\\w([[:digit:]]+)"
	String matNo,rowNo
	SplitString/E=(expr) clipName, matNo, rowNo
	String wName = "ont_" + matNo
	Wave theAngleW = $wName

	String newName = ReplaceString("clp_rgb_",clipName,"clp_rot_")
	Duplicate/O clip1, $newName
	Wave rotClip = $newName
	ImageRotate/A=(theAngleW[str2num(rowNo)][2])/O/Q/RGBA=(0,0,0) rotClip
	
	Variable xSizeOrig = DimSize(clip1,0)
	Variable ySizeOrig = DimSize(clip1,1)	// should be the same as xSize
	Variable diagPx = ceil(sqrt(xSizeOrig^2 + ySizeOrig^2))
	if(mod(diagPx,2) == 0)
		diagPx += 1
	endif
	// final image will be diagpx x diag px, what is its current size after rotation
	Variable xSizeNew = DimSize(rotClip,0)
	Variable ySizeNew = DimSize(rotClip,1)	// should be the same as xSize
	Variable xExpand = floor((diagPx - xSizeNew) / 2)
	Variable yExpand = floor((diagPx - ySizeNew) / 2)
	// pad bottom right
	ImageTransform/N={(xExpand),(yExpand)}/O padimage rotClip
	// rotate by 180 and pad again
	ImageRotate/O/F rotClip
	xSizeNew = DimSize(rotClip,0)
	ySizeNew = DimSize(rotClip,1)
	xExpand = (diagPx - xSizeNew)
	yExpand = (diagPx - ySizeNew)
	ImageTransform/N={(xExpand),(yExpand)}/O padimage rotClip
	// rotate back to original orientation
	ImageRotate/O/F rotClip
	return rotClip
End

// Work on many files from a directory
Function WorkOnDirectory()
	Wave gVarWave = root:gVarWave
	Wave/T pathWave = root:pathWave
	Variable nCh = gVarWave[0]
	Make/O/N=3 ScalingW = {gVarWave[2],gVarWave[3],gVarWave[4]}
	MatrixTranspose ScalingW
	
	NewDataFolder/O/S root:data
	String expDiskFolderName, expDataFolderName
	String FileList, ThisFile
	Variable i
	// XML Dir now pre-selected
	String pathToFolder = pathWave[0]
	NewPath/O/Q ExpDiskFolder, pathToFolder
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
		ScaleAllWaves(3,nCh)
		RotateAndSitUp(0)
		DistanceCalculations()
		MakeWavesForGraphAndPlot(i)
		MakeIndividualGizmo()
		SetDataFolder root:data:
	endfor
	SetDataFolder root:
End

// Work on many files from a directory
Function WorkOnDirectoryIA()
	Wave gVarWave = root:gVarWave
	Wave/T pathWave = root:pathWave
	Variable nCh = gVarWave[0]
	Variable rr = gVarWave[1]
	Make/O/N=3 ScalingW = {gVarWave[2],gVarWave[3],gVarWave[4]}
	MatrixTranspose ScalingW
		
	NewDataFolder/O/S root:data
	String expDiskFolderName, expDataFolderName
	String FileList, ThisFile
	Variable nWaves, i
	// XML Dir now pre-selected
	String pathToFolder = pathWave[0]
	NewPath/O/Q ExpDiskFolder, pathToFolder
	PathInfo /S ExpDiskFolder
	ExpDiskFolderName = S_path
	FileList = IndexedFile(expDiskFolder,-1,".xml")
	Variable nFiles = ItemsInList(FileList)
	Make/O/N=(nFiles)/T root:fileNameWave
	Wave/T fileNameWave = root:fileNameWave
	// TIFF Dir pre-selected
	pathToFolder = pathWave[1]
	NewPath/O/Q ImageDiskFolder, pathToFolder
	PathInfo /S ImageDiskFolder
	String ImageDiskFolderName = S_path
	String ImageFileName
	
	for(i = 0; i < nFiles; i += 1)
		ThisFile = StringFromList(i,FileList)
		fileNameWave[i] = ThisFile
		expDataFolderName = "dataset_" + num2str(i)
		NewDataFolder/O/S $expDataFolderName
		ReadXML(ExpDiskFolderName,ThisFile)
		DealWithDuplicates()
		WAVE/Z/T fNameWave
		ImageFileName = fNameWave[0]
		ScaleAllWaves(1,nCh)
		if(LoadImageAndAnalyse(ImageDiskFolderName,ImageFileName,nCh,rr) != -1)
			MakeIntWavesForGraphAndPlot(i)
		endif
		SetDataFolder root:data:
	endfor
	SetDataFolder root:
End

// Work on many files from a directory
Function WorkOnDirectoryIC()
	Wave gVarWave = root:gVarWave
	Wave/T pathWave = root:pathWave
	Variable nCh = gVarWave[0]
	Variable clipSize = gVarWave[1]
	Make/O/N=3 ScalingW = {gVarWave[2],gVarWave[3],gVarWave[4]}
	MatrixTranspose ScalingW
	
	NewDataFolder/O/S root:data
	String expDiskFolderName, expDataFolderName
	String FileList, ThisFile
	Variable nWaves, i
	// XML Dir now pre-selected
	String pathToFolder = pathWave[0]
	NewPath/O/Q ExpDiskFolder, pathToFolder
	PathInfo /S ExpDiskFolder
	ExpDiskFolderName = S_path
	FileList = IndexedFile(expDiskFolder,-1,".xml")
	Variable nFiles = ItemsInList(FileList)
	Make/O/N=(nFiles)/T root:fileNameWave
	Wave/T fileNameWave = root:fileNameWave
	// TIFF Dir pre-selected
	pathToFolder = pathWave[1]
	NewPath/O/Q ImageDiskFolder, pathToFolder
	PathInfo /S ImageDiskFolder
	String ImageDiskFolderName = S_path
	String ImageFileName
	
	for(i = 0; i < nFiles; i += 1)
		ThisFile = StringFromList(i,FileList)
		fileNameWave[i] = ThisFile
		expDataFolderName = "dataset_" + num2str(i)
		NewDataFolder/O/S $expDataFolderName
		ReadXML(ExpDiskFolderName,ThisFile)
		DealWithDuplicates()
		WAVE/Z/T fNameWave
		ImageFileName = fNameWave[0]
		ScaleAllWaves(1,nCh)
		OrientationOfObjects()
		LoadImageAndClip(ImageDiskFolderName,ImageFileName,nCh,clipSize,i)
		SetDataFolder root:data:
	endfor
	SetDataFolder root:
End

// Work on many files from a directory
Function WorkOnDirectoryER(VARIABLE OptVar)
	Wave gVarWave = root:gVarWave
	Wave/T pathWave = root:pathWave
	// not using row 1 of gVarWave
	Variable nCh = gVarWave[0]
	Make/O/N=3 ScalingW = {gVarWave[2],gVarWave[3],gVarWave[4]}
	MatrixTranspose ScalingW
		
	NewDataFolder/O/S root:data
	String expDiskFolderName, expDataFolderName
	String FileList, ThisFile
	Variable nWaves, i
	// XML Dir now pre-selected
	String pathToFolder = pathWave[0]
	NewPath/O/Q ExpDiskFolder, pathToFolder
	PathInfo /S ExpDiskFolder
	ExpDiskFolderName = S_path
	FileList = IndexedFile(expDiskFolder,-1,".xml")
	Variable nFiles = ItemsInList(FileList)
	Make/O/N=(nFiles)/T root:fileNameWave
	Wave/T fileNameWave = root:fileNameWave
	// TIFF Dir pre-selected
	pathToFolder = pathWave[1]
	NewPath/O/Q ImageDiskFolder, pathToFolder
	PathInfo /S ImageDiskFolder
	String ImageDiskFolderName = S_path
	String ImageFileName
	
	for(i = 0; i < nFiles; i += 1)
		ThisFile = StringFromList(i,FileList)
		fileNameWave[i] = ThisFile
		expDataFolderName = "dataset_" + num2str(i)
		NewDataFolder/O/S $expDataFolderName
		ReadXML(ExpDiskFolderName,ThisFile)
		DealWithDuplicates()
		WAVE/Z/T fNameWave
		ImageFileName = fNameWave[0]
		ScaleAllWaves(1,nCh) // this doesn't scale but it corrects Z
		LoadThresholdedImage(ImageDiskFolderName,ImageFileName)
		ScaleAllWaves(2,nCh) // this will scale without correcting Z
		RotateAndSitUp(optVar) // 1 is scale to ER (excl zone) 0 is scale to Poles (hybrid)
		DistanceCalculations()
		MakeIndividualGizmo()
		TidyUp()
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
	// these ratio waves will be here if ER or hybrid is run
	targetWavelist += "rat_1;rat_2;rat_4;rat_5;"
	Variable nTargets = ItemsInList(targetWaveList)
	String targetName, tList, conName
	
	SetDataFolder root:
	String fullName, modtList
	
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

Function SummariseIntensityMeasurements(nCh)
	Variable nCh
	SetDataFolder root:data:	// relies on earlier load
	DFREF dfr = GetDataFolderDFR()
	String folderName
	Variable numDataFolders = CountObjectsDFR(dfr, 4)
	String wList = "", sList = ""
	
	Variable i,j
	// assemble a string of semi-colon separated targets in the data folder
	for(i = 0; i < numDataFolders; i += 1)
		folderName = GetIndexedObjNameDFR(dfr, 4, i)
		wList += "root:data:" + folderName + ":thisWave;"
	endfor
	
	// these are the waves we want to summarise
	String targetWaveList = "Int_1;Int_2;Int_4;Int_5;"
	Variable nTargets = ItemsInList(targetWaveList)
	String targetName, tList, wName
	
	SetDataFolder root:
	String fullName
	
	for(i = 0; i < nTargets; i += 1)
		targetName = StringFromList(i,targetWaveList)
		tList = ReplaceString("thisWave",wList,targetName)
		wName = "root:summary_" + targetName
		Make/O/N=(numDataFolders,nCh) $wName
		Wave summaryW = $wName
		sList += wName + ";"
		for(j = 0; j < numDataFolders; j +=1)
			fullName = StringFromList(j, tList)
			Wave testW = $fullName
			if(!WaveExists(testW))
				summaryW[j][] = NaN
			else
				MatrixOp/O/FREE meanW = averageCols(testW)
				summaryW[j][] = meanW[0][q]
			endif
		endfor
	endfor
	// this gives us summary waves for each of the objects, each column is a separate channel
	Concatenate/O/KILL/NP=2 sList, summaryMatrix
	// now we have cells in rows, channels in columns and objects in layers
	ImageTransform/G=1 transposeVol summaryMatrix
	WAVE/Z M_VolumeTranspose
	SplitWave/SDIM=2/O/N=channel M_VolumeTranspose
	KillWaves/Z M_VolumeTranspose
	Graph2DWaves(nCh)
End

STATIC Function Graph2DWaves(nCh)
	Variable nCh
	SetDataFolder root:
	String wName, plotName, newName
	Variable nRows
	Make/O/N=(nCh) xPos=p
	Wave/T objWave = root:objWave
	Make/O/N=(nCh)/T xLabel={objWave[0],objWave[1],objWave[3],objWave[4]}
	Wave/T ChWave = root:ChWave
	
	Variable i
	
	for(i = 0; i < nCh; i += 1)
		wName = "channel" + num2str(i)
		plotName = "mean_" + num2str(i)
		newName = "p_mean_" + num2str(i)
		Wave w0 = $wName
		nRows = DimSize(w0,0)
		Make/O/N=(nRows,nCh) xW
		xW[][] = q + gnoise(0.1)
		Redimension/N=(nRows * nCh) w0, xW
		Concatenate/O/KILL {xW,w0}, $newName
		Wave w1 = $newName
		KillWindow/Z $plotName
		Display/N=$plotName/HIDE=1 w1[][1] vs w1[][0]
		ModifyGraph/W=$plotName userticks(bottom)={xPos,xLabel}
		SetAxis/W=$plotName/A/N=1/E=1 left
		Label/W=$plotName left "Average intensity"
		ModifyGraph/W=$plotName mode=3,marker=19,mrkThick=0,rgb=(0,0,0,32768)
		SetAxis/W=$plotName bottom -0.5,(nCh-0.5)
		TextBox/W=$plotName/C/N=text0/F=0/X=0.00/Y=0.00 ChWave[i]
	endfor
End

STATIC Function PlotRatios()
	Variable Obj = 4
	SetDataFolder root:
	String wName, newName
	Variable nRows
	Wave/T objWave = root:objWave
	Make/O/N=(Obj)/T xLabel={objWave[0],objWave[1],objWave[3],objWave[4]}
	String ObjString = "1;2;4;5;"
	String plotName = "ratioBoxPlot"
	WAVE/Z colorWave = MakeColorWave()

	Variable i
	
	for(i = 0; i < Obj; i += 1)
		wName = "all_rat_" + StringFromList(i,ObjString)
		Wave w0 = $wName
		WaveTransform zapnans w0
		w0[] = log(w0[p]) / log(2)
		nRows = Max(nRows,numpnts(w0))
	endfor
	
	KillWindow/Z $plotName
	Display/N=$plotName
	
	for(i = 0; i < Obj; i += 1)
		wName = "all_rat_" + StringFromList(i,ObjString)
		Wave w0 = $wName
		newName = "all_bp_" + StringFromList(i,ObjString)
		Make/O/N=(nRows,Obj) $newName = NaN
		Wave m0 = $newName
		m0[0,numpnts(w0) - 1][i] = w0[p]
		AppendBoxPlot/W=$plotName m0 vs xLabel
		ModifyBoxPlot/W=$plotName trace=$newName,markers={19,-1,19},markerSizes={2,2,2},markerThick={0,0,0}
		ModifyBoxPlot/W=$plotName trace=$newName,whiskerMethod=4
		ModifyGraph/W=$plotName rgb($newName)=(colorWave[i][0],colorWave[i][1],colorWave[i][2],colorWave[i][3])
	endfor
	
	ModifyGraph/W=$plotName toMode=-1
	Label/W=$plotName left "Ratio (Log\B2\M)"
	SetAxis/A/N=1/W=$plotName left
	ModifyGraph/W=$plotName zero(left)=4
End


STATIC Function MakeTheLayouts(prefix,nRow,nCol,[iter, filtVar])
	String prefix
	Variable nRow, nCol
	Variable iter	// this is if we are doing multiple iterations of the same layout
	Variable filtVar // this is the object we want to filter for
	if(ParamIsDefault(filtVar) == 0)
		String filtStr = prefix + "_*_" + num2str(filtVar) + "_*"	// this is if we want to filter for this string from the prefix
	endif
	
	String layoutName = "all"+prefix+"Layout"
	DoWindow/K $layoutName
	NewLayout/N=$layoutName
	String allList = WinList(prefix+"_*",";","WIN:1")
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
	nWindows = ItemsInList(modList)
	Variable PlotsPerPage = nRow * nCol
	String exString = "Tile/A=(" + num2str(ceil(PlotsPerPage/nCol)) + ","+num2str(nCol)+")"
	
	Variable pgNum=1
	
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
	String fileName
	if(!ParamIsDefault(iter))
		fileName = layoutName + num2str(iter) + ".pdf"
	else
		fileName = layoutName + ".pdf"
	endif
	if(ParamIsDefault(filtVar) == 0)
		fileName = ReplaceString(".pdf",fileName, "_" + num2str(filtVar) + ".pdf")
	endif
	String folderStr
	// specific to this ipf
	WAVE/Z/T PathWave
	NewPath/C/O/Q/Z outputFolder PathWave[1]
	SavePICT/O/WIN=$layoutName/P=outputFolder/PGR=(1,-1)/E=-2/W=(0,0,0,0) as fileName
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
	Variable nGiz = ItemsInList(gizList), bigNum
	String modList = wList + ReplaceString("all_rn_",wList,"all_rnOct_") + ReplaceString("all_rn_",wList,"all_rnCube_")
	String pictName
	
	for(i = 0; i < nGiz; i += 1)
		gizName = StringFromList(i,gizList)
		KillWindow/Z $gizName
		NewGizmo/N=$gizName
		// add lights and properties first
		AppendToGizmo/N=$gizName attribute specular={1,1,1,1,1032},name=specular0
		AppendToGizmo/N=$gizName attribute shininess={20,40},name=shininess0
		AppendToGizmo/N=$gizName light=Directional,name=light0
		ModifyGizmo/N=$gizName modifyObject=light0,objectType=light,property={ position,0.4822,-0.5171,0.7071,0.0000}
		ModifyGizmo/N=$gizName modifyObject=light0,objectType=light,property={ direction,0.4822,-0.5171,0.7071}
		ModifyGizmo/N=$gizName setDisplayList=0, object=light0
		ModifyGizmo/N=$gizName setDisplayList=1, attribute=shininess0
		ModifyGizmo/N=$gizName setDisplayList=2, attribute=specular0
		// select scatters
		Wave gW1 = $StringFromList(i*4+0,modList)
		Wave gW2 = $StringFromList(i*4+1,modList)
		Wave gW4 = $StringFromList(i*4+2,modList)
		Wave gW5 = $StringFromList(i*4+3,modList)
		// add in reverse order
		AppendToGizmo/N=$gizName/D Scatter=gW5,name=scatter5
		AppendToGizmo/N=$gizName/D Scatter=gW4,name=scatter4
		AppendToGizmo/N=$gizName/D Scatter=gW2,name=scatter2
		AppendToGizmo/N=$gizName/D Scatter=gW1,name=scatter1
		ModifyGizmo/N=$gizName ModifyObject=scatter1,objectType=scatter,property={ size,0.05}
		ModifyGizmo/N=$gizName ModifyObject=scatter2,objectType=scatter,property={ size,0.4}
		ModifyGizmo/N=$gizName ModifyObject=scatter4,objectType=scatter,property={ size,0.2}
		ModifyGizmo/N=$gizName ModifyObject=scatter5,objectType=scatter,property={ size,0.2}
      ModifyGizmo/N=$gizName ModifyObject=scatter1,objectType=scatter,property={ color,0.5,0.5,0.5,1}
      ModifyGizmo/N=$gizName ModifyObject=scatter2,objectType=scatter,property={ color,0,0,0,1}
      ModifyGizmo/N=$gizName ModifyObject=scatter4,objectType=scatter,property={ color,90/255,96/255,255/255,1}
      ModifyGizmo/N=$gizName ModifyObject=scatter5,objectType=scatter,property={ color,255/255,113/255,24/255,1}
		bigNum = max(wavemax(gW1),wavemax(gW2),wavemax(gW4),wavemax(gW5))
		bigNum = max(bigNum,abs(min(wavemin(gW1),wavemin(gW2),wavemin(gW4),wavemin(gW5))))+0.1
		if(numtype(bigNum) > 0)
			bigNum = 10
		endif
		if(i == 1)
			ModifyGizmo/N=$gizName setOuterBox={-0.1,bigNum,-0.1,bigNum,-0.1,bigNum}
		else
			ModifyGizmo/N=$gizName setOuterBox={-bigNum,bigNum,-bigNum,bigNum,-bigNum,bigNum}
		endif
		ModifyGizmo/N=$gizName scalingOption=0
		// do some rotations and save - see notes below
		DoRotationsAndSave(gizName)
	endfor
	// make and save "empty gizmo showing the axis cue"
	gizName = "axisCue"
	KillWindow/Z $gizName
	NewGizmo/N=$gizName
	ModifyGizmo/N=$gizName showAxisCue=1
	DoRotationsAndSave(gizName)
End

// Notes about rotations
// In the RotateAndSitUp function, poles were aligned in y and then rotated so that they became Z
// So, x = X, y = Z, z = Y
// Our natural view of the spindle is xy so we need to see XZ
// Note that earlier versions of this code we looked at YZ
// The notation I will use is +X+Z this means positive X on x (points right), positive Z on y (points up), i.e. Y away
STATIC Function DoRotationsAndSave(gizName)
	String gizName
	String pictName
	// +X+Z
	ModifyGizmo/N=$gizName setQuaternion={sin(pi/4),0,0,cos(pi/4)}
	pictName = gizName + "_+X+Z.png"
	SavePICT/O/P=TiffPath/E=-5/B=300 as pictName
	// +X+Y
	ModifyGizmo/N=$gizName setQuaternion={0,0,0,0}
	pictName = gizName + "_+X+Y.png"
	SavePICT/O/P=TiffPath/E=-5/B=300 as pictName
	// +Y+Z
	ModifyGizmo/N=$gizName setQuaternion={0.5,0.5,0.5,0.5}
	pictName = gizName + "_+Y+Z.png"
	SavePICT/O/P=TiffPath/E=-5/B=300 as pictName
	// upper quadrant 45 degree view
	Make/O/N=4/FREE q1,q2,qr
	MakeQuaternion(1,0,0,pi/2,q1) // 90 CCW about X
	MakeQuaternion(0,1,0,0.75 * pi,q2) // 135 CCW about Y
	MultiplyQuaternions(q1,q2,qr)
	q1[] = qr[p]
	MakeQuaternion(1,0,0,1.75 * pi,q2) // 315 CCW about X
	MultiplyQuaternions(q1,q2,qr)
	ModifyGizmo/N=$gizName setQuaternion={qr[0],qr[1],qr[2],qr[3]}
	pictName = gizName + "_UQ.png"
	SavePICT/O/P=TiffPath/E=-5/B=300 as pictName
	// set back to +X+Z view
	ModifyGizmo/N=$gizName setQuaternion={sin(pi/4),0,0,cos(pi/4)}
End

Function MakeIndividualGizmo()
	String wList = WaveList("rn_*",";","")
	Variable nWaves = ItemsInList(wList)
	String wName, newName
	Variable nRowsOctant
	Variable i
	
	for(i = 0; i < nWaves; i += 1)
		wName = StringFromList(i,wList)
		Wave w0 = $wName
		newName = ReplaceString("rn_",wName,"rnCube_") // for the whole cube
		nRowsOctant = dimsize(w0,0)
		Make/O/N=(nRowsOctant*8,3) $newName
		Wave w1 = $newName
		w1[][] = w0[mod(p,nRowsOctant)][q]
		MakeCube(w1)
	endfor
	
	String gizName = "g_" + GetDataFolder(0)
	wList = ReplaceString("rn_",wList,"rnCube_")
	String pictName
	
	KillWindow/Z $gizName
	NewGizmo/N=$gizName
	Wave/Z gW1 = $StringFromList(0,wList) // first two should always be present
	Wave/Z gW2 = $StringFromList(1,wList)
	Wave/Z gW4 = $StringFromList(2,wList)
	Wave/Z gW5 = $StringFromList(3,wList)
	// add in reverse order
	if(WaveExists(gW5) == 1)
		AppendToGizmo/N=$gizName/D Scatter=gW5,name=scatter5
		ModifyGizmo/N=$gizName ModifyObject=scatter5,objectType=scatter,property={ size,0.2}
		ModifyGizmo/N=$gizName ModifyObject=scatter5,objectType=scatter,property={ color,255/255,113/255,24/255,1}
	endif
	if(WaveExists(gW4) == 1)
		AppendToGizmo/N=$gizName/D Scatter=gW4,name=scatter4
		ModifyGizmo/N=$gizName ModifyObject=scatter4,objectType=scatter,property={ size,0.2}
		ModifyGizmo/N=$gizName ModifyObject=scatter4,objectType=scatter,property={ color,90/255,96/255,255/255,1}
	endif
	AppendToGizmo/N=$gizName/D Scatter=gW2,name=scatter2
	ModifyGizmo/N=$gizName ModifyObject=scatter2,objectType=scatter,property={ size,0.4}
	ModifyGizmo/N=$gizName ModifyObject=scatter2,objectType=scatter,property={ color,0,0,0,1}
	AppendToGizmo/N=$gizName/D Scatter=gW1,name=scatter1
	ModifyGizmo/N=$gizName ModifyObject=scatter1,objectType=scatter,property={ size,0.05}
	ModifyGizmo/N=$gizName ModifyObject=scatter1,objectType=scatter,property={ color,0.5,0.5,0.5,1}
	
//	ModifyGizmo/N=$gizName insertDisplayList=0, attribute=blendFunc0
//	AppendToGizmo/N=$gizName attribute blendFunction={770,771},name=blendFunc0
//	ModifyGizmo/N=$gizName insertDisplayList=0, opName=enableBlend, operation=enable, data=3042
	ModifyGizmo/N=$gizName scalingOption=0
	ModifyGizmo/N=$gizName setQuaternion={0.5,0.5,0.5,0.5}
	RescaleGizmo(2,gizName)
	pictName = gizName + ".png"
	SavePICT/O/P=TiffPath/E=-5/B=300 as pictName
	KillWindow/Z $gizName
	KillWaves/Z gW1,gW2,gW4,gW5
End

Function RescaleGizmo(boxSize,gizName)
	Variable boxSize
	String gizName
	ModifyGizmo/N=$gizName setOuterBox={-boxSize,boxSize,-boxSize,boxSize,-boxSize,boxSize}
	AppendToGizmo/D/Z/N=$gizName sphere={1/boxSize,25,25},name=sphere0
	ModifyGizmo/N=$gizName modifyObject=sphere0, objectType=sphere, property={radius,1/boxSize}
	ModifyGizmo/N=$gizName modifyObject=sphere0,objectType=Sphere,property={drawStyle,100011}
	ModifyGizmo/N=$gizName modifyObject=sphere0,objectType=Sphere,property={color,0.600000,0.000015,0.480003,1.000000}
End

Function CollateImagesToAverage(nSample,nIter)
	Variable nSample // to determine sample size
	Variable nIter
	
	SetDataFolder root:data:	// relies on earlier load
	DFREF dfr = GetDataFolderDFR()
	String folderName
	Variable numDataFolders = CountObjectsDFR(dfr, 4)
	if(numDataFolders == 0)
		return -1
	endif
	String wList = ""
	String subList
	
	Variable i,ii
	// assemble a string of semi-colon separated targets in all data folders
	for(i = 0; i < numDataFolders; i += 1)
		folderName = "root:data:" + GetIndexedObjNameDFR(dfr, 4, i)
		SetDataFolder folderName
		subList = WaveList("clp_rot_*",";","")
		subList = ReplaceString("clp_",subList,folderName + ":" + "clp_")
		wList += subList
		SetDataFolder root:data:
	endfor
	SetDataFolder root:
	String targets = "1;2;4;5;"
	Variable nTargets = ItemsInList(targets)
	String theTarget,catList,wName,plotname,imgList,newName
	
	// check that we have enough to do the average of nSample and adjust if not
	for(i = 0; i < nTargets; i += 1)
		theTarget = StringFromList(i,targets)
		catList = ListMatch(wList, "*rot_" + theTarget + "_*")
		nSample = min(ItemsInList(catList),nSample)
	endfor
	
	for(ii = 0; ii < nIter; ii += 1)
		for(i = 0; i < nTargets; i += 1)
			theTarget = StringFromList(i,targets)
			catList = ListMatch(wList, "*rot_" + theTarget + "_*")
			Make/O/N=(ItemsInList(catList))/FREE allRowsW = p
			StatsSample/N=(nSample)/Q allRowsW
			WAVE/Z W_Sampled
			Make/O/N=(nSample)/T/FREE allRotW = StringFromList(W_Sampled[p],catList)
			// reassign catList to the sampled imgs
			wfprintf catList, "%s;", allRotW
			wName = "all_ave_" + theTarget
			Concatenate/O catList, $wName
			Wave w0 = $wName
			ImageTransform averageRGBimages w0
			WAVE/Z M_AverageRGBIMage
			KillWaves/Z w0
			Duplicate/O M_AverageRGBImage, $wName
			Wave w0 = $wName
			plotName = ReplaceString("all_ave_",wName,"p_ave_")
			KillWindow/Z $plotName
			NewImage/N=$plotName/HIDE=1/S=0 M_AverageRGBIMage
			ModifyGraph/W=$plotName width={Plan,1,top,left}
			RGB2Montage(w0,4)
		endfor
		MakeTheLayouts("p_mtg",6,1,iter = ii)
		KillTheLayout("p_mtg",0)
	endfor
End

Function CollateImagesToAverageByDistance(nSample,nIter)
	Variable nSample // to determine sample size
	Variable nIter

	SetDataFolder root:data:	// relies on earlier load
	DFREF dfr = GetDataFolderDFR()
	String folderName
	Variable numDataFolders = CountObjectsDFR(dfr, 4)
	if(numDataFolders == 0)
		return -1
	endif
	String wList = ""
	String subList
	
	Variable i,ii
	// assemble a string of semi-colon separated targets in all data folders
	for(i = 0; i < numDataFolders; i += 1)
		folderName = "root:data:" + GetIndexedObjNameDFR(dfr, 4, i)
		SetDataFolder folderName
		subList = WaveList("clp_rot_*",";","")
		subList = ReplaceString("clp_",subList,folderName + ":" + "clp_")
		wList += subList
		SetDataFolder root:data:
	endfor
	SetDataFolder root:
	String targets = "1;2;4;5;"
	Variable nTargets = ItemsInList(targets)
	String theTarget,catList,wName,plotname,imgList,newName
	// first check whether nSample exceeds the upper half of mat_5 objects by distance to nearest pole
	catList = ListMatch(wList, "*rot_5_*")
	nSample = min(floor(ItemsInList(catList) / 2),nSample)
	
	// now check that we have enough to do the average of nSample and adjust if not
	for(i = 0; i < nTargets; i += 1)
		theTarget = StringFromList(i,targets)
		catList = ListMatch(wList, "*rot_" + theTarget + "_*")
		nSample = min(ItemsInList(catList),nSample)
	endfor
	
	for(ii = 0; ii < nIter; ii += 1)
		for(i = 0; i < nTargets; i += 1)
			theTarget = StringFromList(i,targets)
			catList = ListMatch(wList, "*rot_" + theTarget + "_*")
			if(cmpstr(theTarget,"5") == 0)
				Wave allRowsW = GiveMeTheValidWaveIndex(catList)
			else
				Make/O/N=(ItemsInList(catList))/FREE allRowsW = p
			endif
			StatsSample/N=(nSample)/Q allRowsW
			WAVE/Z W_Sampled
			Make/O/N=(nSample)/T/FREE allRotW = StringFromList(W_Sampled[p],catList)
			// reassign catList to the sampled imgs
			wfprintf catList, "%s;", allRotW
			wName = "all_ave_" + theTarget
			Concatenate/O catList, $wName
			Wave w0 = $wName
			ImageTransform averageRGBimages w0
			WAVE/Z M_AverageRGBIMage
			KillWaves/Z w0
			Duplicate/O M_AverageRGBImage, $wName
			Wave w0 = $wName
			plotName = ReplaceString("all_ave_",wName,"p_ave_")
			KillWindow/Z $plotName
			NewImage/N=$plotName/HIDE=1/S=0 M_AverageRGBIMage
			ModifyGraph/W=$plotName width={Plan,1,top,left}
			RGB2Montage(w0,4)
		endfor
		MakeTheLayouts("p_mtg",6,1,iter = ii)
		KillTheLayout("p_mtg",0)
	endfor
End

STATIC Function/WAVE GiveMeTheValidWaveIndex(catList)
	String catList
	Variable nImages = ItemsInList(catList)
	Variable i
	
	String imgName, ontName // although it says Name, these are the full address
	String expr = "root\\Wdata\\Wdataset\\w([[:digit:]]+)\\Wclp\\wrot\\w([[:digit:]]+)\\w([[:digit:]]+)"
	String datasetNo,matNo,rowNo
	Make/O/N=(nImages)/FREE distW, indexW=p
	
	for(i = 0; i < nImages; i += 1)
		imgName = StringFromList(i,catList)
		SplitString/E=(expr) imgName, datasetNo, matNo, rowNo
		ontName = "root:data:dataset_" + datasetNo + ":ont_" + matNo
		Wave ontW = $ontName
		distW[i] = ontW[str2num(rowNo)][0]
	endfor
	/// now we have a wave with all the distances in it in order with catList
	Sort/R distW, distW, indexW
	Make/O/N=(floor(nImages / 2)) root:furthest5W/WAVE=indexWave
	indexWave[] = indexW[p] 
	return indexWave
End

////////////////////////////////////////////////////////////////////////
// Utility functions
////////////////////////////////////////////////////////////////////////

STATIC Function CleanSlate()
	SetDataFolder root:
	String fullList = WinList("*", ";","WIN:65543")
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

// Use this to kill layouts and contents made with MakeTheLayouts
STATIC Function KillTheLayout(prefix,killWavesToo)
	String prefix
	Variable killWavesToo // 0 means no, 1 (or anything else) means yes)
	String layoutName = "all"+prefix+"Layout"
	DoWindow/K $layoutName
	
	SetDataFolder root:
	String fullList = WinList(prefix + "_*", ";","WIN:1")
	Variable allItems = ItemsInList(fullList)
	String windowName, imgName, imgWaveName
	Variable i
 	
 	if(killWavesToo == 0)
		for(i = 0; i < allItems; i += 1)
			windowName = StringFromList(i, fullList)
			KillWindow/Z $windowName	
		endfor
 	else
		for(i = 0; i < allItems; i += 1)
			windowName = StringFromList(i, fullList)
			imgName = StringFromList(0, ImageNameList(windowName,";"))
			Wave imgWave = ImageNameToWaveRef(windowName,imgName)
			KillWindow/Z $windowName	
			KillWaves/Z imgWave
		endfor
	endif
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

STATIC Function RGB2Montage(masterImage,grout)
	Wave masterImage
	Variable grout
	Variable nRows = 1
	Variable nColumns = 4
	
	String mtgName = "mtg_" + NameOfWave(masterImage)
	
	if(!WaveExists(masterImage))
		Print "Image does not exist"
		return 0
	endif
	
	// This is for RGB images only
	Variable cSize = dimsize(masterImage,2)
	Variable xSize = dimsize(masterImage,1)
	Variable ySize = dimsize(masterImage,0)
	Variable x1 = (xSize * nColumns) + (grout * (nColumns-1))
	Variable y1 = (ySize * nRows) + (grout * (nRows-1))
	if(cSize != 3)
		Print NameOfWave(masterImage) + " is not an RGB image."
		return -1
	else
		Make/W/U/O/N=(x1,y1,3) $mtgName
	endif
	Wave newMontage = $mtgName
	// make it white
	newMontage[][] = 65535
	
	Variable xPos=0,yPos=0
	Variable i,j
	
	for(i = 0; i < nColumns; i += 1)
		xPos = mod(i,nColumns)
		yPos = floor(i/nColumns)
		x1 = (xSize * xPos) + (grout * xPos)
		y1 = (ySize * yPos) + (grout * yPos)
		// at this position we paste in the channels
		for(j = 0; j < 3; j += 1)
			if(i == 3)
				Duplicate/O/FREE/RMD=[][][j] masterImage, subImage
			else
				Duplicate/O/FREE/RMD=[][][i] masterImage, subImage
			endif
			ImageTransform/P=(j)/INSI=subImage/INSX=(x1)/INSY=(y1) InsertImage newMontage
		endfor
	endfor
	String plotName = "p_mtg_" + NameOfWave(masterImage)
	KillWindow/Z $plotName
	NewImage/N=$plotName/S=0 newMontage
	ModifyGraph/W=$plotName width={Plan,1,top,left}
End

STATIC Function Resample3DLine(matA,acDist)
	Wave matA // 3-column wave of x,y,z end points
	Variable acDist // distance from c to A (i.e. centre of spindle to spindle pole)
	WAVE/Z cWave
	// Find length between the two points
	Duplicate/O/FREE matA, mTemp0
	mTemp0[][] = cWave[0][q] - matA[p][q]
	MatrixOp/O/FREE mTemp1 = sumrows(mTemp0 * mTemp0)
	MatrixOp/O/FREE mTemp2 = sqrt(mTemp1)
	// Make a matrix to hold the resampled line with xyz coords
	Variable nLayer = ceil(acDist * 3) // number of pixels for 3 half-spindle lengths
	String mName = NameOfWave(matA) + "_rs"
	// 3D wave with rows as different objects, columns as xyz and layers as interpolated positions
	Variable nRow = DimSize(MatA,0)
	Make/O/N=(nRow,3,nLayer) $mName
	Wave matB = $mName
	// fill in values in the wave
	matB[][][] = ((1 - (r / mTemp2[p])) * cWave[0][q]) + ((r / mTemp2[p]) * matA[p][q])
End


STATIC Function InterpolateAndFindThreshold(m0,imgMat)
	Wave m0, imgMat	// m0 is the original coordinate set.
	String mName = NameOfWave(m0) + "_rs"
	Wave m1 = $mName // 3D wave with the xyz coords of all interpolated lines
	Variable nObjects = DimSize(m1,0)
	// Make a wave to hold the threshold xyz coords	
	Make/O/N=(nObjects,3) $(ReplaceString("_rs",mName,"_th"))
	Wave threshCoordsMat = $(ReplaceString("_rs",mName,"_th"))
	// make a copy of the 3D wave and transpose so that it's more simple to extract layers for objects
	MatrixOp/O m2 = transposeVol(m1,3)
	
	Variable i,rowNo
	
	for(i = 0; i < nObjects; i += 1)
		MatrixOP/O/FREE tempMat = m2[][][i] // take the correct coords for the object
		Wave w0 = CallInterp3d(ImgMat, tempMat)	// could not get Interp3DPath to work
//		FindLevel/Q/P w0, 128
		Findlevel/B=2/EDGE=1/Q/P w0,255
		if(V_Flag == 1)
			Print mName, num2str(i), "threshold not found"
			threshCoordsMat[i][] = NaN
		elseif(numtype(V_levelX) == 2)
			threshCoordsMat[i][] = NaN
		else
			rowNo = ceil(V_levelX)
			threshCoordsMat[i][] = tempMat[V_levelX][q]
		//	Print mName, num2str(i), "threshold found"
		endif
	endfor
	// Make a wave to hold the ratio of distance from c to obj : c to thresh
	Make/O/N=(nObjects) $(ReplaceString("mat_",NameOfWave(m0),"rat_"))
	Wave ratioW = $(ReplaceString("mat_",NameOfWave(m0),"rat_"))
	WAVE/Z cWave
	// Find length between the objects in m0
	Duplicate/O/FREE m0, mTemp0
	mTemp0[][] = cWave[0][q] - m0[p][q]
	MatrixOp/O/FREE mTemp01 = sumrows(mTemp0 * mTemp0)
	MatrixOp/O/FREE mTemp02 = sqrt(mTemp01)
	// Find length between the objects in m1
	Duplicate/O/FREE threshCoordsMat, mTemp1
	mTemp1[][] = cWave[0][q] - threshCoordsMat[p][q]
	MatrixOp/O/FREE mTemp11 = sumrows(mTemp1 * mTemp1)
	MatrixOp/O/FREE mTemp12 = sqrt(mTemp11)
	ratioW[] = mTemp02 / mTemp12
	
	KillWaves/Z m2, resultW, m1, threshCoordsMat
End

STATIC Function/WAVE CallInterp3D(matA,matB)
	Wave matA,matB
	Variable nCall = DimSize(matB,0)
	Make/O/N=(nCall) resultW
	Variable i
	for(i = 0; i < nCall; i += 1)
		resultW[i] = Interp3D(matA, matB[i][0],matB[i][1],matB[i][2])
	endfor
	return resultW
End

// Note that this was written for XML where the following were picked out (by index)
//1 Kinetochores Plate
//2 Spindle Poles
//3 BackGround
//4 Misaligned Non-ensheathed"
//5 Misaligned Ensheathed

////////////////////////////////////////////////////////////////////////
// Panel functions
////////////////////////////////////////////////////////////////////////

Function StartingPanelForXML()
	// make global text wave to store paths, object and channel info
	Make/T/O/N=2 PathWave
	Make/T/O/N=5 ObjWave={"Aligned","Poles","Background","Misaligned","Ensheathed"}
	Make/T/O/N=4 ChWave={"DNA","ER","Tubulin","CENP-C"}
	// make global numeric wave for other variables
	Make/O/N=8 gVarWave={4,8,0.06449999660253525,0.06449999660253525,0.20000000298023224}
	// note that this panel will not deal with less/more than 4 channels and 5 objects
	DoWindow/K SetUp
	NewPanel/N=SetUp/K=1/W=(81,73,774,298)
	SetDrawLayer UserBack
	SetDrawEnv linefgc= (65535,65535,65535),fillfgc= (49151,53155,65535)
	DrawRect 9,116,176,218
	Button SelectDir1,pos={12,10},size={140,20},proc=ButtonProcXML,title="Select XML Dir"
	Button SelectDir2,pos={12,41},size={140,20},proc=ButtonProcXML,title="Select Output Dir"
	SetVariable Dir1,pos={188,13},size={480,14},value= PathWave[0],title="XML Directory"
	SetVariable Dir2,pos={188,44},size={480,14},value= PathWave[1],title="Output Directory"
	SetVariable Obj0,pos={214,80},size={194,14},value= ObjWave[0],title="Object 1"
	SetVariable Obj1,pos={214,110},size={194,14},value= ObjWave[1],title="Object 2"
	SetVariable Obj2,pos={214,140},size={194,14},value= ObjWave[2],title="Object 3"
	SetVariable Obj3,pos={214,170},size={194,14},value= ObjWave[3],title="Object 4"
	SetVariable Obj4,pos={214,200},size={194,14},value= ObjWave[4],title="Object 5"
	
	SetVariable ChSetVar,pos={12,70},size={166,14},title="How many channels?"
	SetVariable ChSetVar,format="%g",value= gVarWave[0]
//	SetVariable RrSetVar,pos={12,90},size={166,14},title="Analysis radius (px)"
//	SetVariable RrSetVar,format="%g",value= gVarWave[1]
	SetVariable xVar,pos={27,135},size={126,14},title="x size (nm)"
	SetVariable xVar,format="%g",value= gVarWave[2]
	SetVariable yVar,pos={27,160},size={126,14},title="y size (nm)"
	SetVariable yVar,format="%g",value= gVarWave[3]
	SetVariable zVar,pos={27,185},size={126,14},title="z size (nm)"
	SetVariable zVar,format="%g",value= gVarWave[4]

	SetVariable Ch0,pos={468,80},size={194,14},value= ChWave[0],title="Channel 1"
	SetVariable Ch1,pos={468,110},size={194,14},value= ChWave[1],title="Channel 2"
	SetVariable Ch2,pos={468,140},size={194,14},value= ChWave[2],title="Channel 3"
	SetVariable Ch3,pos={468,170},size={194,14},value= ChWave[3],title="Channel 4"
	Button DoIt,pos={564,194},size={100,20},proc=ButtonProcXML,title="Do It"
End
 
// define buttons
Function ButtonProcXML(ctrlName) : ButtonControl
	String ctrlName
 
		Wave/T PathWave, ObjWave, ChWave
		Wave gVarWave
		Variable refnum,okVar
 
		strswitch(ctrlName)
 
			case "SelectDir1"	:
				// set XML directory
				NewPath/Q/O/M="Locate folder with XML files" XMLPath
				PathInfo XMLPath
				PathWave[0] = S_Path
				PathWave[1] = S_Path
				break
 
			case "SelectDir2"	:
				// set TIFF directory
				NewPath/Q/O/M="Locate output folder with images" TiffPath
				PathInfo TiffPath
				PathWave[1] = S_Path
				break
 
			case "DoIt" :
				// check CondWave
				okvar = WaveChecker(PathWave)
				if (okvar == -1)
					Print "Error: Not all directories are selected."
					break
				endif
				okvar = NameChecker(objWave)
				if (okvar == -1)
					Print "Error: Two objects have the same name."
					break
				endif
				okvar = NameChecker(ChWave)
				if (okvar == -1)
					Print "Error: Two channels have the same name."
					break
				else
					XMLWrapperFunc()
					return -1
				endif
		EndSwitch
End

Function StartingPanelForIA()
	// make global text wave to store paths, object and channel info
	Make/T/O/N=2 PathWave
	Make/T/O/N=5 ObjWave={"Aligned","Poles","Background","Misaligned","Ensheathed"}
	Make/T/O/N=4 ChWave={"DNA","Tubulin","ER","CENP-C"}
	// make global numeric wave for other variables
	Make/O/N=8 gVarWave={4,8,0.06449999660253525,0.06449999660253525,0.20000000298023224}
	// note that this panel will not deal with less/more than 4 channels and 5 objects
	DoWindow/K SetUp
	NewPanel/N=SetUp/K=1/W=(81,73,774,298)
	SetDrawLayer UserBack
	SetDrawEnv linefgc= (65535,65535,65535),fillfgc= (49151,53155,65535)
	DrawRect 9,116,176,218
	Button SelectDir1,pos={12,10},size={140,20},proc=ButtonProcIA,title="Select XML Dir"
	Button SelectDir2,pos={12,41},size={140,20},proc=ButtonProcIA,title="Select TIFF Dir"
	SetVariable Dir1,pos={188,13},size={480,14},value= PathWave[0],title="XML Directory"
	SetVariable Dir2,pos={188,44},size={480,14},value= PathWave[1],title="TIFF Directory"
	SetVariable Obj0,pos={214,80},size={194,14},value= ObjWave[0],title="Object 1"
	SetVariable Obj1,pos={214,110},size={194,14},value= ObjWave[1],title="Object 2"
	SetVariable Obj2,pos={214,140},size={194,14},value= ObjWave[2],title="Object 3"
	SetVariable Obj3,pos={214,170},size={194,14},value= ObjWave[3],title="Object 4"
	SetVariable Obj4,pos={214,200},size={194,14},value= ObjWave[4],title="Object 5"
	
	SetVariable ChSetVar,pos={12,70},size={166,14},title="How many channels?"
	SetVariable ChSetVar,format="%g",value= gVarWave[0]
	SetVariable RrSetVar,pos={12,90},size={166,14},title="Analysis radius (px)"
	SetVariable RrSetVar,format="%g",value= gVarWave[1]
	SetVariable xVar,pos={27,135},size={126,14},title="x size (nm)"
	SetVariable xVar,format="%g",value= gVarWave[2]
	SetVariable yVar,pos={27,160},size={126,14},title="y size (nm)"
	SetVariable yVar,format="%g",value= gVarWave[3]
	SetVariable zVar,pos={27,185},size={126,14},title="z size (nm)"
	SetVariable zVar,format="%g",value= gVarWave[4]

	SetVariable Ch0,pos={468,80},size={194,14},value= ChWave[0],title="Channel 1"
	SetVariable Ch1,pos={468,110},size={194,14},value= ChWave[1],title="Channel 2"
	SetVariable Ch2,pos={468,140},size={194,14},value= ChWave[2],title="Channel 3"
	SetVariable Ch3,pos={468,170},size={194,14},value= ChWave[3],title="Channel 4"
	Button DoIt,pos={564,194},size={100,20},proc=ButtonProcIA,title="Do It"
End
 
// define buttons
Function ButtonProcIA(ctrlName) : ButtonControl
	String ctrlName
 
		Wave/T PathWave, ObjWave, ChWave
		Wave gVarWave
		Variable refnum,okVar
 
		strswitch(ctrlName)
 
			case "SelectDir1"	:
				// set XML directory
				NewPath/Q/O/M="Locate folder with XML files" XMLPath
				PathInfo XMLPath
				PathWave[0] = S_Path
				break
 
			case "SelectDir2"	:
				// set TIFF directory
				NewPath/Q/O/M="Locate folder with images" TiffPath
				PathInfo TiffPath
				PathWave[1] = S_Path
				break
 
			case "DoIt" :
				// check CondWave
				okvar = WaveChecker(PathWave)
				if (okvar == -1)
					Print "Error: Not all directories are selected."
					break
				endif
				okvar = NameChecker(objWave)
				if (okvar == -1)
					Print "Error: Two objects have the same name."
					break
				endif
				okvar = NameChecker(ChWave)
				if (okvar == -1)
					Print "Error: Two channels have the same name."
					break
				else
					IAWrapperFunc()
					return -1
				endif
		EndSwitch
End

Function StartingPanelForIC()
	// make global text wave to store paths, object and channel info
	Make/T/O/N=2 PathWave
	Make/T/O/N=5 ObjWave={"Aligned","Poles","Background","Misaligned","Ensheathed"}
	Make/T/O/N=4 ChWave={"DNA","Tubulin","ER","CENP-C"}
	// make global numeric wave for other variables
	Make/O/N=8 gVarWave={4,61,0.06449999660253525,0.06449999660253525,0.20000000298023224,2,1,3}
	// note that this panel will not deal with less/more than 4 channels and 5 objects
	DoWindow/K SetUp
	NewPanel/N=SetUp/K=1/W=(81,73,774,298)
	SetDrawEnv linefgc= (65535,65535,65535),fillfgc= (49151,53155,65535)
	DrawRect 9,116,176,218
	SetDrawEnv fillfgc= (65535,0,0)
	DrawOval 527,80,542,95
	SetDrawEnv fillfgc= (0,65535,0)
	DrawOval 527,110,542,125
	SetDrawEnv fillfgc= (0,0,65535)
	DrawOval 527,140,542,155
	Button SelectDir1,pos={12,10},size={140,20},proc=ButtonProcIC,title="Select XML Dir"
	Button SelectDir2,pos={12,41},size={140,20},proc=ButtonProcIC,title="Select TIFF Dir"
	SetVariable Dir1,pos={188,13},size={480,14},title="XML Directory"
	SetVariable Dir1,value= PathWave[0]
	SetVariable Dir2,pos={188,44},size={480,14},title="TIFF Directory"
	SetVariable Dir2,value= PathWave[1]
	SetVariable Obj0,pos={214,80},size={132,14},title="Object 1"
	SetVariable Obj0,value= ObjWave[0]
	SetVariable Obj1,pos={214,110},size={132,14},title="Object 2"
	SetVariable Obj1,value= ObjWave[1]
	SetVariable Obj2,pos={214,140},size={132,14},title="Object 3"
	SetVariable Obj2,value= ObjWave[2]
	SetVariable Obj3,pos={214,170},size={132,14},title="Object 4"
	SetVariable Obj3,value= ObjWave[3]
	SetVariable Obj4,pos={214,200},size={132,14},title="Object 5"
	SetVariable Obj4,value= ObjWave[4]
	SetVariable ChSetVar,pos={12,70},size={166,14},title="How many channels?"
	SetVariable ChSetVar,format="%g",value= gVarWave[0]
	SetVariable RrSetVar,pos={12,90},size={166,14},title="Clip size (px)"
	SetVariable RrSetVar,format="%g",value= gVarWave[1]
	SetVariable xVar,pos={27,135},size={126,14},title="x size (nm)"
	SetVariable xVar,format="%g",value= gVarWave[2]
	SetVariable yVar,pos={27,160},size={126,14},title="y size (nm)"
	SetVariable yVar,format="%g",value= gVarWave[3]
	SetVariable zVar,pos={27,185},size={126,14},title="z size (nm)"
	SetVariable zVar,format="%g",value= gVarWave[4]
	SetVariable Ch0,pos={371,80},size={139,14},title="Channel 1"
	SetVariable Ch0,value= ChWave[0]
	SetVariable Ch1,pos={371,110},size={139,14},title="Channel 2"
	SetVariable Ch1,value= ChWave[1]
	SetVariable Ch2,pos={371,140},size={139,14},title="Channel 3"
	SetVariable Ch2,value= ChWave[2]
	SetVariable Ch3,pos={371,170},size={139,14},title="Channel 4"
	SetVariable Ch3,value= ChWave[3]
	SetVariable rCh,pos={549,80},size={126,14},title="Red",format="%g"
	SetVariable rCh,value= gVarWave[5]
	SetVariable gCh,pos={549,110},size={126,14},title="Green"
	SetVariable gCh,format="%g",value= gVarWave[6]
	SetVariable bCh,pos={549,140},size={126,14},title="Blue",format="%g"
	SetVariable bCh,value= gVarWave[7]
	Button DoIt,pos={564,194},size={100,20},proc=ButtonProcIC,title="Do It"
End
 
// define buttons
Function ButtonProcIC(ctrlName) : ButtonControl
	String ctrlName
 
		Wave/T PathWave, ObjWave, ChWave
		Wave gVarWave
		Variable refnum,okVar
 
		strswitch(ctrlName)
 
			case "SelectDir1"	:
				// set XML directory
				NewPath/Q/O/M="Locate folder with XML files" XMLPath
				PathInfo XMLPath
				PathWave[0] = S_Path
				break
 
			case "SelectDir2"	:
				// set TIFF directory
				NewPath/Q/O/M="Locate folder with images" TiffPath
				PathInfo TiffPath
				PathWave[1] = S_Path
				break
 
			case "DoIt" :
				// check CondWave
				okvar = WaveChecker(PathWave)
				if (okvar == -1)
					Print "Error: Not all directories are selected."
					break
				endif
				okvar = NameChecker(objWave)
				if (okvar == -1)
					Print "Error: Two objects have the same name."
					break
				endif
				okvar = NameChecker(ChWave)
				if (okvar == -1)
					Print "Error: Two channels have the same name."
					break
				else
					// make sure clip size is not an even number
					if(mod(gVarWave[1],2) == 0)
						gVarWave[1] += 1
					endif
					ICWrapperFunc()
					return -1
				endif
		EndSwitch
End

Function StartingPanelForER()
	// make global text wave to store paths, object and channel info
	WAVE/Z/T PathWave
	if(!WaveExists(PathWave))
		Make/T/O/N=2 PathWave
	endif
	Make/T/O/N=5 ObjWave={"Aligned","Poles","Background","Misaligned","Ensheathed"}
	Make/T/O/N=4 ChWave={"DNA","Tubulin","ER","CENP-C"}
	// make global numeric wave for other variables
	Make/O/N=8 gVarWave={4,8,0.06449999660253525,0.06449999660253525,0.20000000298023224}
	// note that this panel will not deal with less/more than 4 channels and 5 objects
	DoWindow/K SetUp
	NewPanel/N=SetUp/K=1/W=(81,73,774,298)
	SetDrawLayer UserBack
	SetDrawEnv linefgc= (65535,65535,65535),fillfgc= (49151,53155,65535)
	DrawRect 9,116,176,218
	Button SelectDir1,pos={12,10},size={140,20},proc=ButtonProcER,title="Select XML Dir"
	Button SelectDir2,pos={12,41},size={140,20},proc=ButtonProcER,title="Select TIFF Dir"
	SetVariable Dir1,pos={188,13},size={480,14},value= PathWave[0],title="XML Directory"
	SetVariable Dir2,pos={188,44},size={480,14},value= PathWave[1],title="TIFF Directory"
	SetVariable Obj0,pos={214,80},size={194,14},value= ObjWave[0],title="Object 1"
	SetVariable Obj1,pos={214,110},size={194,14},value= ObjWave[1],title="Object 2"
	SetVariable Obj2,pos={214,140},size={194,14},value= ObjWave[2],title="Object 3"
	SetVariable Obj3,pos={214,170},size={194,14},value= ObjWave[3],title="Object 4"
	SetVariable Obj4,pos={214,200},size={194,14},value= ObjWave[4],title="Object 5"
	
	SetVariable ChSetVar,pos={12,70},size={166,14},title="How many channels?"
	SetVariable ChSetVar,format="%g",value= gVarWave[0]
//	SetVariable RrSetVar,pos={12,90},size={166,14},title="Analysis radius (px)"
//	SetVariable RrSetVar,format="%g",value= gVarWave[1]
	SetVariable xVar,pos={27,135},size={126,14},title="x size (nm)"
	SetVariable xVar,format="%g",value= gVarWave[2]
	SetVariable yVar,pos={27,160},size={126,14},title="y size (nm)"
	SetVariable yVar,format="%g",value= gVarWave[3]
	SetVariable zVar,pos={27,185},size={126,14},title="z size (nm)"
	SetVariable zVar,format="%g",value= gVarWave[4]

	SetVariable Ch0,pos={468,80},size={194,14},value= ChWave[0],title="Channel 1"
	SetVariable Ch1,pos={468,110},size={194,14},value= ChWave[1],title="Channel 2"
	SetVariable Ch2,pos={468,140},size={194,14},value= ChWave[2],title="Channel 3"
	SetVariable Ch3,pos={468,170},size={194,14},value= ChWave[3],title="Channel 4"
	Button Hybrid,pos={444,194},size={100,20},proc=ButtonProcER,title="Hybrid"
	Button DoIt,pos={564,194},size={100,20},proc=ButtonProcER,title="Do It"
End
 
// define buttons
Function ButtonProcER(ctrlName) : ButtonControl
	String ctrlName
 
		Wave/T PathWave, ObjWave, ChWave
		Wave gVarWave
		Variable refnum,okVar
 
		strswitch(ctrlName)
 
			case "SelectDir1"	:
				// set XML directory
				NewPath/Q/O/M="Locate folder with XML files" XMLPath
				PathInfo XMLPath
				PathWave[0] = S_Path
				break
 
			case "SelectDir2"	:
				// set TIFF directory
				NewPath/Q/O/M="Locate folder with images" TiffPath
				PathInfo TiffPath
				PathWave[1] = S_Path
				break
				
			case "Hybrid" :
				// check CondWave
				okvar = WaveChecker(PathWave)
				if (okvar == -1)
					Print "Error: Not all directories are selected."
					break
				endif
				okvar = NameChecker(objWave)
				if (okvar == -1)
					Print "Error: Two objects have the same name."
					break
				endif
				okvar = NameChecker(ChWave)
				if (okvar == -1)
					Print "Error: Two channels have the same name."
					break
				else
					ERWrapperFunc(0)
					return -1
				endif
				
			case "DoIt" :
				// check CondWave
				okvar = WaveChecker(PathWave)
				if (okvar == -1)
					Print "Error: Not all directories are selected."
					break
				endif
				okvar = NameChecker(objWave)
				if (okvar == -1)
					Print "Error: Two objects have the same name."
					break
				endif
				okvar = NameChecker(ChWave)
				if (okvar == -1)
					Print "Error: Two channels have the same name."
					break
				else
					ERWrapperFunc(1)
					return -1
				endif
		EndSwitch
End

STATIC function WaveChecker(TextWaveToCheck)
	Wave/T TextWaveToCheck
	Variable nRows = numpnts(TextWaveToCheck)
	Variable len
	
	Variable i
	
	for(i = 0; i < nRows; i += 1)
		len = strlen(TextWaveToCheck[i])
		if(len == 0)
			return -1
		elseif(numtype(len) == 2)
			return -1
		endif
	endfor
	return 1
End

STATIC function NameChecker(TextWaveToCheck)
	Wave/T TextWaveToCheck
	Variable nRows = numpnts(TextWaveToCheck)
	Variable len
	
	Variable i,j
	
	for(i = 0; i < nRows; i += 1)
		for(j = 0; j < nRows; j += 1)
			if(j > i)
				if(cmpstr(TextWaveToCheck[i], TextWaveToCheck[j], 0) == 0)
					return -1
				endif
			endif
		endfor
	endfor
	return 1
End

Function VizSphereMask()
	WAVE/Z sphereMask = root:sphereMask
	Duplicate/O/RMD=[][][][0]/FREE sphereMask, theTempMat
	Redimension/N=(-1,-1,-1) theTempMat
	// rows are x values, columns are y values and layers are z-values
	Variable xSize = DimSize(theTempMat,0)
	Variable ySize = DimSize(theTempMat,1)
	Variable zSize = DimSize(theTempMat,2)
	Variable nRows = xSize * ySize * zSize
	Make/O/N=(nRows,3) cubicLoc
	cubicLoc[][0] = mod(p,xSize)
	cubicLoc[][1] = mod(floor(p/ySize),xSize)
	cubicLoc[][2] = floor(p/(xSize * ySize))
	Make/O/N=(nRows) cubicVal
	cubicVal[] = theTempMat[cubicLoc[p][0]][cubicLoc[p][1]][cubicLoc[p][2]]
	MatrixOp/O colorWaveGiz = colRepeat(cubicVal,4)
	// set zeros to blank
	colorWaveGiz[][3] = (cubicVal[p] == 0) ? 0 : 1
	colorWaveGiz[][1] = 0
	MatrixOp/O sizeWaveGiz = colRepeat(cubicVal,3)
	sizeWaveGiz[][] = (cubicVal[p] == 0) ? 0.1 : 0.75
	NewGizmo
	AppendToGizmo DefaultScatter= root:cubicLoc
	ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ sizeType,1}
	ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ sizeWave,root:sizeWaveGiz}
	ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ scatterColorType,1}
	ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ colorWave,root:colorWaveGiz}
End

STATIC Function/WAVE MakeColorWave()
	// This is for the box plot
	// rows are the four items on box plot
	// originally 1,2,4,5 were grey,red,green,blue
	// now they are grey,orange,blue,red
	Make/O/N=(4,4) colorWave = 0
	colorWave[][3] = 127 // alpha is 0.5
	colorWave[0][0,2] = 127 // grey
	colorWave[1][0,2] = {{0},{0},{0}} // black
	colorWave[2][0,2] = {{90},{96},{255}} // blue
	colorWave[3][0,2] = {{255},{113},{24}} // orange
	colorWave *= 257
	return colorWave
End

// q1 and q2 are waves that correspond to {x,y,z,w} quaternions.
// New quaternion stored in qr which is the quaternion product q1*q2.
Function MultiplyQuaternions(q1,q2,qr)
	Wave q1,q2,qr
	
	Variable w2 = q2[3]
	Variable w1 = q1[3]
	qr[3] = w2 * w1 - (q2[0] * q1[0] + q2[1] * q1[1] + q2[2] * q1[2])
	Make/N=4/FREE vcross = 0
 	vcross[0] = (q1[1] * q2[2]) - (q1[2] * q2[1])
 	vcross[1] = (q1[2] * q2[0]) - (q1[0] * q2[2])
 	vcross[2] = (q1[0] * q2[1]) - (q1[1] * q2[0])
 	MatrixOP/FREE aa = w2 * q1 + w1 * q2 + vcross
	qr[0] = aa[0]
	qr[1] = aa[1]
	qr[2] = aa[2]
	Variable NN = norm(qr)
	qr /= NN
End

Function MakeQuaternion(Ax,Ay,Az,theta,w)
	Variable Ax,Ay,Az,theta
	Wave w
	Variable N = sqrt(Ax^2 + Ay^2 + Az^2)
	w[0] = Ax * sin(theta/2)/N
	w[1] = Ay * sin(theta/2)/N
	w[2] = Az * sin(theta/2)/N
	w[3] = cos(theta/2)
End