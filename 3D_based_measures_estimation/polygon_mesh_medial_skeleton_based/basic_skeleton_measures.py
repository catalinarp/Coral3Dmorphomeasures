import vtk
import numpy as np
import vtk.util.numpy_support as nps

def getJunctionEndPointIds(poly_line):
	"""
	Obtains point ids of end vertices (d = 1), junction vertices (d = 3) and 
	vertices that are end of a branch but have only d=2 (merge ids). Also the
	cell ids of the branches with end vertices are obtained.

	Args:
		poly_line (vtk.vtkPolyData): skeleton graph

	Returns:
		Dict: dictionary with following key-value pairs:
		- junctions: all point ids at branch junctions (np.array, int)
		- end_points: all point ids of end points (np.array, int)
		- end_branches: all cell ids of branches with end points (np.array, int)
		- merge_ids: all point ids of vertices where branches should be merged
		  (np.array, int)
	"""
	n_cells = poly_line.GetNumberOfCells()
	in_cells = vtk.vtkIdList()
	junction_ids = []
	end_point_ids = []
	end_branch_ids = []
	merge_ids = []

	# loop over every branch
	for i in range(n_cells):
		branch = poly_line.GetCell(i)
		n_points = branch.GetNumberOfPoints()
		
		# check the first and last vertex of branch
		for j in [0, n_points - 1]:
			point_id = branch.GetPointId(j)

			# check in how many cells these vertices are present
			poly_line.GetPointCells(point_id, in_cells)
			n_neigh = in_cells.GetNumberOfIds()
			
			# 1: end point, 2: branch should be merges, 3: junction point
			if n_neigh == 1:
				end_point_ids.append(point_id)
				end_branch_ids.append(i)
			elif n_neigh == 2:
				merge_ids.append(point_id)
			else:
				junction_ids.append(point_id)

	# get the set of junction IDs to prevent overlap
	junction_ids = np.asarray(list(set(junction_ids)), dtype =int)
	end_point_ids = np.asarray(end_point_ids, dtype=int)
	end_branch_ids = np.asarray(end_branch_ids, dtype=int)
	if merge_ids:
		merge_ids = np.asarray(merge_ids, dtype = int)
	return {
				'junctions': junction_ids, 
				'end_points': end_point_ids, 
				'end_branches':end_branch_ids, 
				'merge_ids': merge_ids
		}

######## IMPROVE

def getBasicSkelMeasures(poly_line):
	"""
	Get some basic measures of a skeleton graph.
	NOTE: the maximum and minimum thickness described here is something else 
	than da and db (see skel_thickness.py)
	
	Args:
		poly_line (vtk.vtkPolyData): skeleton graph

	Returns:
		Dict: dictrionary with following keys and values:
		- 
	"""
	junc_end = getJunctionEndPointIds(poly_line)
	thickness = getThicknessDistribution(poly_line)
	min_max_avg = getMinMaxAvgThick(poly_line, thickness['medial_thickness'])
	ep_raw = junc_end['end_points']
	eb_raw = junc_end['end_branches']
	root_id = selectRootId(min_max_avg['avg_thick'], eb_raw, ep_raw)
	junc_end['end_points'] = np.delete(ep_raw, np.where(ep_raw == root_id['root_point']))
	junc_end['end_branches'] = np.delete(eb_raw, np.where(eb_raw == root_id['root_branch']))
	return {**junc_end, **thickness, **min_max_avg, **root_id}

def selectRootId(avg_thickness, end_branch_ids,end_point_ids):
	"""
	Select root branch and vertex by selecting the end branch with the highest 
	average medial thickness. 


	Args:
		avg_thickness (np.array,float): average thickness of every branch 
		end_branch_ids (np.array, int): cell ids of all end branches
		end_point_ids (np.array, int): cell ids of all end points
	
	NOTE: the order of end_branch_ids and end_point_ids should be ordered, such
	that 'end_branch_ids[i]' is the cell id that containt point id 
	'end_point_ids[i]'.
	'
	Returns:
		Dict: dictionary with following key-value pairs
		- 
	"""
	end_thick = avg_thickness[end_branch_ids]
	root_idx = np.argmax(end_thick)
	return {'root_point':end_point_ids[root_idx], 'root_branch':end_branch_ids[root_idx]}

def getThicknessDistribution(poly_line):
	"""
	Obtain the scalars (meant to be medial thickness) of a skeleton graph as a 
	numpy array.

	Args:
		poly_line (vtk.vtkPolyData): skeleton graph

	Returns:
		Dict: dictionary with key-value:
		- medial_thickness: point scalars ordered by point id (np.array, float)
	"""
	scalars = poly_line.GetPointData().GetScalars()
	scalars_np = nps.vtk_to_numpy(scalars)
	return {'medial_thickness':scalars_np}

def getIdsBranch(branch):
	n_points = branch.GetNumberOfPoints()
	ids = np.zeros(n_points, dtype = int)
	for i in range(n_points):
		ids[i] = branch.GetPointId(i)
	return ids

def getMinMaxAvgThickIds(thickness, pids):
	thickvals = thickness[pids]
	return [np.min(thickvals),np.max(thickvals),np.mean(thickvals)]

def getMinMaxAvgThick(poly_line, thickness):
	"""
	"""
	n_lines = poly_line.GetNumberOfCells()
	min_max_avg = np.zeros((n_lines,3))
	for i in range(n_lines):
		line = poly_line.GetCell(i)
		ids = getIdsBranch(line)
		min_max_avg[i] = getMinMaxAvgThickIds(thickness, ids)
	min_max_avg = min_max_avg.T
	return {
			'min_thick': min_max_avg[0], 
			'max_thick': min_max_avg[1], 
			'avg_thick': min_max_avg[2]
		}