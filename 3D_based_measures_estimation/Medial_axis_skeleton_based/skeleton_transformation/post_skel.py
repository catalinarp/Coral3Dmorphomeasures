import vtk
import getJunctionEndPointIds

import numpy as np
import vtk.util.numpy_support as nps
import math as m

def postProcessSkeleton(line_data, poly_data, min_size_end=4):
	
	# assign medial thickness to vertices
	line_data = cleanBranches(line_data, min_size_end)




def cleanBranches(line_data, min_size_end = 4):
	
	br = getJunctionEndPointIds(line_data)
	end_br = br['end_branches']
	end_p = br['end_points']	
	line_data = reverseEndBranches(line_data, end_p, end_br)
	line_data = deleteShortEndBranches(line_data, end_br, min_end_br_length)
	line_data = mergeBranches(line_data)
	return line_data


def cleanSkeleton(line_data, poly_data=None, change_thick = True, \
		min_end_br_length = 4):
	"""
	Clean skeleton graph: short terminal branches are removed and connecting 
	branches are merged, end branch vertices are ordered (end point is point ID 
	0 of cell). Medial thickness can be redefined. Unconnected components are
	removed. 
	
	"""

	# remove 
	line_data = cleanBranches(line_data, min_end_br_length)
	if change_thick:
		line_data = tryNewThickness(line_data, poly_data)
	line_data = preprocessPolyData(line_data, mm_to_cm = False)
	br = getJunctionEndPointIds(line_data)

	# somehow the end branches need to be reversed again
	line_data = reverseEndBranches(line_data, br['end_points'],br['end_branches'])
	return line_data



def reverseEndBranches(poly_line,end_points, end_branches):
	for i,br in enumerate(end_branches):
		branch = poly_line.GetCell(br)
		n_points = branch.GetNumberOfPoints()
		if branch.GetPointId(0) != end_points[i]:
			poly_line.ReverseCell(br)
	return poly_line

def mergeBranches(line_data):
	complete = False
	while not complete:
		complete = True
		merge_ids = getJunctionEndPointIds(line_data)['merge_ids']
		if not len(merge_ids):
			return line_data
		to_delete = []
		lines = []
		in_cells = vtk.vtkIdList()
		compid = 0
		for pid in merge_ids:
			line_data.GetPointCells(pid, in_cells)
			cids = [in_cells.GetId(0), in_cells.GetId(1)]
			for cid in cids:
				if cid in to_delete:
					complete = False
					compid = pid	
			if complete == False and compid == pid:
				compid = 0
				continue

			np_l = [line_data.GetCell(cids[0]).GetNumberOfPoints(),\
				line_data.GetCell(cids[1]).GetNumberOfPoints()]
			new_line = vtk.vtkLine()
			new_line.GetPointIds().SetNumberOfIds(np_l[0] + np_l[1] -1)
			counter = 0
			for c in range(2):
				line = line_data.GetCell(cids[c])
				if c == 0:
					rev = (line.GetPointId(0) == pid)
				else:
					rev = (line.GetPointId(np_l[c]-1) == pid)
				if rev:
					rang = [np_l[c] -1 -i for i in range(np_l[c])]
				else:
					rang = [i for i in range(np_l[c])]
				if c == 1:
					rang.pop(0)
				for i in rang:
					new_line.GetPointIds().SetId(counter, line.GetPointId(i))
					counter += 1
				to_delete.append(cids[c])
			lines.append(new_line)

		cells = vtk.vtkCellArray()
		
		for i in range(line_data.GetNumberOfCells()):
			if i not in to_delete:
				cells.InsertNextCell(line_data.GetCell(i))
		for i in range(len(lines)):
			cells.InsertNextCell(lines[i])

		poly_new = vtk.vtkPolyData()
		poly_new.DeepCopy(line_data)
		poly_new.SetLines(cells)
		line_data = poly_new
		line_data = preprocessPolyData(line_data)
	return line_data

def deleteShortEndBranches(line_data, end_branches, min_points = 4):
	in_cells = vtk.vtkIdList()
	for br in end_branches:
		cel = line_data.GetCell(br)
		num_p = cel.GetNumberOfPoints()
		if num_p < min_points:
			line_data.DeleteCell(br)
	line_data.RemoveDeletedCells()
	return line_data


def tryNewThickness(line_data, poly_data):
	
	n_points = line_data.GetNumberOfPoints()

	point_locator = vtk.vtkPointLocator()
	point_locator.SetDataSet(poly_data)
	point_locator.BuildLocator()

	medial_thickness = np.zeros(n_points)
	for i in range(n_points):
		skel_point_loc = line_data.GetPoint(i)
		poly_id = point_locator.FindClosestPoint(skel_point_loc)
		poly_point_loc = poly_data.GetPoint(poly_id)
		dist = vtk.vtkMath.Distance2BetweenPoints(skel_point_loc, poly_point_loc)
		medial_thickness[i] = m.sqrt(dist) * 2

	line_data.GetPointData().SetScalars(nps.numpy_to_vtk(medial_thickness))
	return line_data

def cleanBranches(line_data, min_end_br_length = 4):
	
	br = getJunctionEndPointIds(line_data)
	end_br = br['end_branches']
	end_p = br['end_points']	
	line_data = reverseEndBranches(line_data, end_p, end_br)
	line_data = deleteShortEndBranches(line_data, end_br, min_end_br_length)
	line_data = mergeBranches(line_data)
	return line_data


def preprocessPolyData(poly_data):
	"""
	Removes unconnected components (artefacts) from vtkPolyData object.

	Args:
		poly_data (vtk.vtkPolyData): mesh or graph with artefacts
		mm_to_cm (bool): mesh coordinates are in mm but should be cm

	Returns:
		vtk.vtkPolyData: mesh without artifacts
	"""
	
	poly_large = getLargestClosedSurface(poly_data)
	poly_clean = cleanPolyData(poly_large)
	return poly_clean

def getLargestClosedSurface(poly_data):
	"""
	Takes the largest connected component of a vtkPolyData object.
	
	Args:
		poly_data (vtk.vtkPolyData): mesh or graph with unconnected components

	Returns:
		vtk.vtkPolyData: mesh or graph single connected component
		
	"""
	conn_filter = vtk.vtkPolyDataConnectivityFilter()
	conn_filter.SetExtractionModeToLargestRegion()
	conn_filter.SetInputData(poly_data)
	conn_filter.Update()
	return conn_filter.GetOutput()

def cleanPolyData(polydata):
	"""
	Removes double and unused points and degenerate cells.

	Args:
		poly_data (vtk.vtkPolyData): mesh or graph

	Returns:
		vtk.vtkPolyData: mesh or graph where all vertices are used
	"""
	clean_filter = vtk.vtkCleanPolyData()
	clean_filter.SetInputData(polydata)
	clean_filter.Update()
	return clean_filter.GetOutput()

