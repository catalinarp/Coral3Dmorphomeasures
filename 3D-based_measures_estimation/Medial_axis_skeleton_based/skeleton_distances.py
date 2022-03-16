import copy 
import math as m 
import numpy as np
import vtk
from .basic_skeleton_measures import getJunctionEndPointIds

def getSkeletonDistances(poly_line, min_length = 4):
	
	"""
	Obtains measures directly related to distances between nodes in a skeleton 
	graph (branch length, branching rate, branch spacing).

	Args:
		poly_line (vtk.vtkPolyData): skeleton graph
		min_length (int): minimum number of vertices in branch to be considered
			for analyses
	
	Returns:
		Dict: dictionary of numpy arrays with keys
			- br_length: branch length of each branch (float,ordered by cell ID)
			- br_rate: brancing rate of each branch (float, ordered by cell ID)
			- long_enough: suitable for analysis (bool, ordered by cell ID)
			- br_spacing_v1: branch spacing v1 (float)
			- br_spacing_v1_loc: vertex IDs of br_spacing_v1
			- br_spacing_v2: branch spacing v2 (float)
			- br_spacing_v2_loc: vertex IDs of br_spacing_v2
	"""

	# get end vertices and corresponding end branches
	skel_info = getJunctionEndPointIds(poly_line)
	end_points = skel_info['end_points']
	end_branches = skel_info['end_branches']

	# get measures related to branch distances 
	lengths = getBranchLengths(poly_line, min_length)
	br_sp1 = getBranchSpacing(poly_line, end_points, end_branches)
	br_sp2 = getEndSpacing(poly_line, end_points, end_branches)
	return {**lengths, **br_sp1, **br_sp2}

def getBranchLengths(poly_line, min_length = 4):
	"""
	Obtain branch length and brancing rate of all branches in a skeleton graph.
	Branch length is sum of all edge lengths in a branch. Branch rate is the 
	distance between the first and last vertex of a branch.

	Args:
		poly_line (vtk.vtkPolydata): skeleton graph
		min_length (int) : number of vertices a branch should have for analysis

	Returns:
		Dict: dictionary of numpy arrays with keys
			- br_length: branch length of each branch (float,ordered by cell ID)
			- br_rate: brancing rate of each branch (float, ordered by cell ID)
			- long_enough: suitable for analysis (bool, ordered by cell ID)
	"""
	
	# obtain number of cells (/branches)
	n_lines = poly_line.GetNumberOfCells()
	lengths = np.zeros(n_lines)
	rates = np.zeros(n_lines)
	enough_points = np.zeros(n_lines, dtype = bool)
	points = poly_line.GetPoints()
	# loop over all branches
	for i in range(n_lines):
		branch = poly_line.GetCell(i)
		
		# obtain branch length, branching rate and suitability for analysis
		lengths[i] = getLineLength(branch,points)
		rates[i] = getLineRate(branch,points)
		enough_points[i] = longEnough(branch, min_length)
	# dictionary of all metrics per branch ID
	return {'br_length': lengths, 'br_rate': rates, 'long_enough': enough_points}

def getLineLength(branch,points):
	"""
	Obtains branch length of a branch (length of all edges).

	Args:
		branch (vtk.vtkCell): branch of skeleton graph
	
	Returns:
		float: branch length
	"""
	n_points = branch.GetNumberOfPoints()
	l = 0
	
	# loop over every edge
	for i in range(1,n_points):
		pid0 = branch.GetPointId(i-1)
		pid1 = branch.GetPointId(i)
		p0 = points.GetPoint(pid0)
		p1 = points.GetPoint(pid1)
		# obtain length and sum to previous lengths
		l += m.sqrt(vtk.vtkMath.Distance2BetweenPoints(p0, p1)) 
	return l

def getLineRate(branch,points):
	"""
	Obtains branching rate of a branch (distance between first and last vertex).

	Args:
		branch (vtk.vtkCell): branch of skeleton graph
	
	Returns:
		float: branching rate
	"""

	n_points = branch.GetNumberOfPoints()
	
	# get first vertex
	pid0 = branch.GetPointId(0)
	p0 = points.GetPoint(pid0)

	# get last vertex
	pid1 = branch.GetPointId(n_points - 1)
	p1 = points.GetPoint(pid1)
	
	# calculate distance
	l = m.sqrt(vtk.vtkMath.Distance2BetweenPoints(p0,p1))
	return l

def longEnough(branch, min_points = 4):
	"""
	Returns whether the number of vertices in a branch meets the minimum amount
	of branch vertices to be suitable for analysis.

	Args:
		branch (vtk.vtkCell): branch of skeleton graph 
		min_points (int): minimum amount of vertices to be long enough
	
	Returns:
		Bool: True if long enough, False otherwise
	"""
	n_points = branch.GetNumberOfPoints()
	return (n_points >= min_points)

def getBranchSpacing(poly_line, end_points, end_branches):
	"""
	Get branch_spacing_v1: the shortest distance between a branch tip and a 
	vertex of the skeleton graph that is not part of the corresponding branch
	(see [1]).

	Args:
		poly_line (vtk.vtkPolyData): skeleton graph
		end_points (np.array,int): vertex IDs of branch tips
		end_branches (np.array,int): cell IDs of branches of end_points
	
	Returns:
		Dict: dictionary of numpy arrays with keys
			- br_spacing_v1: branch spacing v1 (float)
			- br_spacing_v1_loc: vertex IDs of br_spacing_v1


	[1]: OUR PAPER
	"""
	n_ends = len(end_points)
	branch_spacing = np.zeros(n_ends)

	# get all coordinates of vertices in the skeleton graph
	point_coords = getSkelPointLocation(poly_line)
	n_points = point_coords.shape[0]
	
	# loop over every terminal branch
	for i in range(n_ends):
		
		# get branch and point ids in branch
		branch = poly_line.GetCell(end_branches[i])
		pid_branch = getPointIdsBranch(branch)
		ep_coord = point_coords[end_points[i]]

		# select all the points, except those in own branch
		mask = np.ones(n_points, dtype = bool)
		mask[pid_branch] = False 
		point_coord_branch = point_coords[mask]

		# get closest point outside own branch
		branch_spacing[i] = getSpacing(point_coord_branch, ep_coord)
		
	return {'br_spacing_v1': branch_spacing, 'br_spacing_v1_loc': end_points}

def getEndSpacing(poly_line, end_points, end_branches):
	"""
	Get branch_spacing_v2: the shortest distance between a branch tip and any
	other branch tip (see [1]).

	Args:
		poly_line (vtk.vtkPolyData): skeleton graph
		end_points (np.array,int): vertex IDs of branch tips
		end_branches (np.array,int): cell IDs of branches of end_points
	
	Returns:
		Dict: dictionary of numpy arrays with keys
			- br_spacing_v2: branch spacing v2 (float)
			- br_spacing_v2_loc: vertex IDs of br_spacing_v1


	[1]: OUR PAPER
	"""
	n_ends = len(end_points)
	end_point_spacing = np.zeros(n_ends)

	# get coordinates of end point vertices in the graph
	point_coords = getSkelPointLocation(poly_line)
	print(point_coords)
	end_coords = copy.deepcopy(point_coords[end_points,:])
	print(end_coords)
	del point_coords
	
	# loop over end points
	for i in range(n_ends):
		print(i)
		# get coordinates of end point
		coord_tip = end_coords[i]
		print(coord_tip)
		# get the other end points
		other_coords = np.delete(end_coords, i,axis=0)
		print(other_coords)
		# find shortest distance
		end_point_spacing[i] = getSpacing(other_coords, coord_tip)

	return {'br_spacing_v2':end_point_spacing, 'br_spacing_v2_loc':end_points}


def getSkelPointLocation(poly_data):
	"""
	Get coordinates of all vertices in a polydata object.

	Args:
		poly_data (vtk.vtkPolyData): polydataobject
	Returns:
		np.array (float): array of coordinates of all vertices (ordered by ID)
	"""
	n_points = poly_data.GetNumberOfPoints()
	coords = np.zeros((n_points,3))
	for i in range(n_points):
		coords[i] = poly_data.GetPoint(i)
	return coords

def getPointIdsBranch(branch):
	"""
	Obtain all vertex IDs in a branch.

	Args:
		branch (vtk.vtkCell): branch of skeleton graph\

	Returns:
		np.array (int): vertex IDs in branch
	"""
	n_points = branch.GetNumberOfPoints()
	branch_points = np.zeros(n_points, dtype = int)
	
	# loop over all branch vertices and get ID
	for i in range(n_points):
		branch_points[i] = branch.GetPointId(i)
	return branch_points

def getSpacing(coords, point):
	"""
	Find sortest distance between a point and a list of coordinates.

	Args:
		coords (np.array, float): array of coordinates
		point (np.array, float): coordinate of reference point
	
	Returns:
		shortest distance between c
	"""

	# get distance between point and all coords
	vecs = coords - point
	vecs_sq = vecs ** 2
	dist = np.sqrt(np.sum(vecs_sq, axis = 1))

	# find shortest distance
	return np.min(dist)
