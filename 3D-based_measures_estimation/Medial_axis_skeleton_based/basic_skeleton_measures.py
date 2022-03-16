import vtk
import numpy as np

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
