from ..helpers.load_data import readVTK
from ..helpers import local_directories as ldir
import copy 
import math as m 
# from ...helpers import vtk_local as vtk 
# import numpy as np
# import copy
# from ...helpers.shape_data.load_data import readVTK
# import math as m

def getSkeletonDistances(coral_name, skel_prop_dict, min_length = 4, \
		line_dir = ldir.LINE):
	
	"""
	Obtains measures directly related to distances between nodes in a skeleton 
	graph (branch length, branching rate, branch spacing).

	Args:
		coral_name (str): sample ID of the coral (filename without subscript)
		skel_prop_dict (Dict): dictonionary that contains TODO
		min_length (int): minimum number of vertices in branch to be considered
			for analyses
		line_dir (str): directory in which the skeleton graph is saved
	
	Returns:
		Dict: 
	"""

	end_points = copy.deepcopy(skel_prop_dict[coral_name]['end_points'])
	end_branches = copy.deepcopy(skel_prop_dict[coral_name]['end_branches'])
	
	del skel_prop_dict
	
	poly_line = readVTK(coral_name, line_dir)
	lengths = getBranchLengths(poly_line, min_length)
	br_sp1 = getBranchSpacing(poly_line, end_points, end_branches)
	br_sp2 = getEndSpacing(poly_line, end_points, end_branches)
	return {**lengths, **br_sp1, **br_sp2}

def getBranchLengths(poly_line, min_length = 4):
	"""
	Obtain branch length and brancing rate of all branches in a skeleton graph.
	Branch length is sum of all edges in a branch.


	"""
	n_lines = poly_line.GetNumberOfCells()
	lengths = np.zeros(n_lines)
	rates = np.zeros(n_lines)
	enough_points = np.zeros(n_lines, dtype = bool)
	points = poly_line.GetPoints()
	for i in range(n_lines):
		branch = poly_line.GetCell(i)
		lengths[i] = getLineLength(branch,points)
		rates[i] = getLineRate(branch,points)
		enough_points[i] = longEnough(branch, min_length)
	return {'br_length': lengths, 'br_rate': rates, 'long_enough': enough_points}


def longEnough(branch, min_points = 4):
	n_points = branch.GetNumberOfPoints()
	return (n_points >= min_points)

def getLineLength(branch, points):
	n_points = branch.GetNumberOfPoints()
	l = 0
	for i in range(1,n_points):
		pid0 = branch.GetPointId(i-1)
		pid1 = branch.GetPointId(i)
		p0 = points.GetPoint(pid0)
		p1 = points.GetPoint(pid1)
		l += m.sqrt(vtk.vtkMath.Distance2BetweenPoints(p0, p1)) 
	return l

def getLineRate(branch, points):
	"""
	"""

	n_points = branch.GetNumberOfPoints()
	pid0 = branch.GetPointId(0)
	pid1 = branch.GetPointId(n_points - 1)
	p0 = points.GetPoint(pid0)
	p1 = points.GetPoint(pid1)
	l = m.sqrt(vtk.vtkMath.Distance2BetweenPoints(p0,p1))
	return l

def getBranchLengths(poly_line, min_length = 4):
	n_lines = poly_line.GetNumberOfCells()
	lengths = np.zeros(n_lines)
	rates = np.zeros(n_lines)
	enough_points = np.zeros(n_lines, dtype = bool)
	points = poly_line.GetPoints()
	for i in range(n_lines):
		branch = poly_line.GetCell(i)
		lengths[i] = getLineLength(branch,points)
		rates[i] = getLineRate(branch,points)
		enough_points[i] = longEnough(branch, min_length)
	return {'br_length': lengths, 'br_rate': rates, 'long_enough': enough_points}

def getSpacing(ps, ep):
	vecs = ps - ep
	vecs_sq = vecs ** 2
	dist = np.sqrt(np.sum(vecs_sq, axis = 1))
	return np.min(dist)


def getSkelPointLocation(poly_data):
	n_points = poly_data.GetNumberOfPoints()
	coords = np.zeros((n_points,3))
	for i in range(n_points):
		coords[i] = poly_data.GetPoint(i)
	return coords

def getPointIdsBranch(branch):
	n_points = branch.GetNumberOfPoints()
	branch_points = np.zeros(n_points, dtype = int)
	for i in range(n_points):
		branch_points[i] = branch.GetPointId(i)
	return branch_points


def getBranchSpacing(poly_data, end_points, end_branches):
	
	n_ends = len(end_points)
	branch_spacing = np.zeros(n_ends)

	point_coords = getSkelPointLocation(poly_data)
	n_points = point_coords.shape[0]
	for i in range(n_ends):
		
		# get branch and point ids in branch
		branch = poly_data.GetCell(end_branches[i])
		pid_branch = getPointIdsBranch(branch)
		ep_coord = point_coords[end_points[i]]

		# select all the points, except those in own branch
		mask = np.ones(n_points, dtype = bool)
		mask[pid_branch] = False 
		point_coord_branch = point_coords[mask]

		# get closest point outside own branch
		branch_spacing[i] = getSpacing(point_coord_branch, ep_coord)
		
	return {'br_spacing_v1': branch_spacing, 'br_spacing_v1_loc': end_points}

def getEndSpacing(poly_data, end_points, end_branches):
	n_ends = len(end_points)
	end_point_spacing = np.zeros(n_ends)

	point_coords = getSkelPointLocation(poly_data)
	for i in range(n_ends):
		ep_coord = point_coords[end_points[i]]
		end_points_temp = np.delete(end_points, i)
		point_coords_temp = point_coords[end_points_temp]
		end_point_spacing[i] = getSpacing(point_coords_temp, ep_coord)

	return {'br_spacing_v2':end_point_spacing, 'br_spacing_v2_loc':end_points}