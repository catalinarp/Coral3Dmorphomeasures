import vtk
import vtk.util.numpy_support as nps
import numpy as np 
import copy
import math as m 
from .basic_skeleton_measures import getBasicSkelMeasures

def getBranchWidthAndAngles(poly_line):
	

	skel_info = getBasicSkelMeasures(poly_line)
	end_branches = skel_info['end_branches']

	max_thick = MaxThickEnds(poly_line, end_branches)
	min_thick = getMinimumThickness(poly_line,skel_info)
	term_thick = getTerminalThickness(poly_line, skel_info)
	end_angles = getEndBranchAngle(poly_line, skel_info, min_thick['db_loc'])
	min_thick_ends = {'db': min_thick['db'][end_branches], 'db_loc': min_thick['db_loc'][end_branches]}
	return {**max_thick, **min_thick_ends, **term_thick, **end_angles}

def MaxThickEnds(poly_line, end_branches): 
	n_branches = len(end_branches)
	ids = np.zeros(n_branches, dtype = int)
	thickness = np.zeros(n_branches)
	for i,br in enumerate(end_branches):
		branch = poly_line.GetCell(br)
		n_points = branch.GetNumberOfPoints()
		pid = branch.GetPointId(n_points- 1)
		thick = poly_line.GetPointData().GetScalars().GetValue(pid)
		ids[i] = pid
		thickness[i] = thick
	return {'da': thickness, 'da_loc': ids}

def getMinimumThickness(poly_line, basic_info): 
	end_branch_ids = copy.deepcopy(basic_info['end_branches'])
	root_id = copy.deepcopy(basic_info['root_point'])
	end_point_ids = copy.deepcopy(basic_info['end_points'])

	del basic_info

	root_loc = poly_line.GetPoint(root_id)
	del root_id, end_point_ids
	n_branch = poly_line.GetNumberOfCells()
	min_thick = np.zeros(n_branch)
	min_id = np.zeros(n_branch, dtype = int)
	
	for i in range(n_branch):
		if i in end_branch_ids:
			ordered = True 
		else:
			ordered = False
		branch = poly_line.GetCell(i)
		bid, thick = getMinimumThicknessBranch(branch, poly_line, root_loc, ordered)

		min_thick[i] = thick 
		min_id[i] = bid

	return {'db': min_thick, 'db_loc':min_id}


def getMinimumThicknessBranch(branch, poly_line, root_loc, ordered = False): 
	n_points = branch.GetNumberOfPoints()
	junc_id = branch.GetPointId(n_points - 1)

	if not ordered:
		dist = np.zeros(2)
		for i,pid in enumerate([0, n_points -1]):
			jid = branch.GetPointId(pid)
			jloc = poly_line.GetPoint(jid)
			dist[i] = vtk.vtkMath.Distance2BetweenPoints(jloc, root_loc)
		min_dist = np.argmin(dist)
		if min_dist == 0:
			junc_id = branch.GetPointId(0)

	junc_loc = poly_line.GetPoint(junc_id)
	junc_d = poly_line.GetPointData().GetScalars().GetValue(junc_id)
	min_sphere = np.zeros(n_points)
	
	for i in range(branch.GetNumberOfPoints()):
		pid = branch.GetPointId(i)
		if pid == junc_id:
			min_sphere[i] = 10000
			continue
		p_loc = poly_line.GetPoint(pid)
		p_d = poly_line.GetPointData().GetScalars().GetValue(pid)

		dist = m.sqrt(vtk.vtkMath.Distance2BetweenPoints(junc_loc, p_loc))
		rad_sum = (junc_d + p_d)/2
		min_sphere[i] = np.abs(dist - rad_sum)

	br_pid = np.argmin(min_sphere)

	pid = branch.GetPointId(br_pid)

	thick = poly_line.GetPointData().GetScalars().GetValue(pid)
	return pid, thick

def getTerminalThickness(poly_line,basic_info):
	end_points = basic_info['end_points']
	thickness = basic_info['medial_thickness']
	return {'dc': thickness[end_points], 'dc_loc':end_points}

def getEndBranchAngle(poly_line, basic_info, db_locs):

	end_branches = copy.deepcopy(basic_info['end_branches'])

	del basic_info
	angles = np.zeros(len(end_branches))
	norms = np.zeros((len(end_branches),2,3))
	juncs = np.zeros(len(end_branches), dtype = int)
	for i,br in enumerate(end_branches):

		angles[i],norms[i],juncs[i] = getBranchAngle(br, db_locs, poly_line)
	return {'angle':angles, 'angle_normals': norms, 'angle_loc': juncs}


def getBranchAngle(branch_id, sphere_ids, poly_line):

	branch = poly_line.GetCell(branch_id)
	n_points = branch.GetNumberOfPoints()
	junc_id = branch.GetPointId(n_points - 1)

	junc_loc = np.asarray(poly_line.GetPoint(junc_id))
	br_loc = np.asarray(poly_line.GetPoint(sphere_ids[branch_id]))

	br_normal = br_loc - junc_loc

	in_cells = vtk.vtkIdList()
	poly_line.GetPointCells(junc_id,in_cells)
	n_cells = in_cells.GetNumberOfIds()
	angles = np.zeros(n_cells)
	norms = np.zeros((n_cells,3))

	for i in range(n_cells):

		neigh_id = in_cells.GetId(i)
		neigh_loc = np.asarray(poly_line.GetPoint(sphere_ids[neigh_id]))
		neigh_normal = neigh_loc - junc_loc
		x = -1
		while np.all((neigh_normal == 0)):
			
			neigh_branch = poly_line.GetCell(neigh_id)
			sp_id = neigh_branch.GetPointId(int(neigh_branch.GetNumberOfPoints()/2)+x)
			neigh_loc = np.asarray(poly_line.GetPoint(sp_id))
			neigh_normal = neigh_loc - junc_loc
			x +=1

		norms[i] = neigh_normal
		angles[i] = calcAngleNormals(br_normal, neigh_normal)

	norms = norms[np.where(angles > 0.001)]
	angles = angles[np.where(angles > 0.0001)]
	norm_set = np.asarray([br_normal, norms[np.argmin(angles)]])

	return np.min(angles),norm_set,junc_id


def calcAngleNormals(normal1, normal2):

	normal_array = np.zeros((2,3))
	normal_array[0] = normal1
	normal_array[1] = normal2
	normal_array = normal_array ** 2
	normal_array = np.sum(normal_array, axis = 1)
	normal_array = np.sqrt(normal_array)

	dot = np.sum(normal1 * normal2)
	cos_angle = dot/(normal_array[0] * normal_array[1])

	if np.abs(cos_angle) <= 1:
		return m.acos(cos_angle)
	else:
		return m.acos(int(cos_angle))

"""def getSpheresAndAngles(coral_name, thickness_dict, line_dir = vtk.DIR_LINE):
	basic_info = copy.deepcopy(thickness_dict[coral_name])
	del thickness_dict
	max_thick = getMaximumThickness(coral_name, basic_info)
	poly_line = readVTK(coral_name, line_dir)
	
	min_thick = getMinimumThickness(poly_line, basic_info)
	term_thick = getTerminalThickness(poly_line, basic_info)

	end_angles = getEndBranchAngle(poly_line, basic_info, min_thick['db_loc'])
	return {**max_thick, **min_thick, **term_thick, **end_angles}


def MaxThickEnds(poly_line, end_branches): # U
	n_branches = len(end_branches)
	ids = np.zeros(n_branches, dtype = int)
	thickness = np.zeros(n_branches)
	for i,br in enumerate(end_branches):
		branch = poly_line.GetCell(br)
		n_points = branch.GetNumberOfPoints()
		pid = branch.GetPointId(n_points- 1)
		thick = poly_line.GetPointData().GetScalars().GetValue(pid)
		ids[i] = pid
		thickness[i] = thick
	return {'da': thickness, 'da_loc': ids}


def getMinimumThickness(poly_line, basic_info): 

	end_branch_ids = copy.deepcopy(basic_info['end_branches'])
	root_id = copy.deepcopy(basic_info['root_point'])
	end_point_ids = copy.deepcopy(basic_info['end_points'])

	del basic_info

	root_loc = poly_line.GetPoint(root_id)
	del root_id, end_point_ids
	n_branch = poly_line.GetNumberOfCells()
	min_thick = np.zeros(n_branch)
	min_id = np.zeros(n_branch, dtype = int)
	
	for i in range(n_branch):
		if i in end_branch_ids:
			ordered = True 
		else:
			ordered = False
		branch = poly_line.GetCell(i)
		bid, thick = getMinimumThicknessBranch(branch, poly_line, root_loc, ordered)

		min_thick[i] = thick 
		min_id[i] = bid

	return {'db': min_thick, 'db_loc':min_id}


def getMinimumThicknessBranch(branch, poly_line, root_loc, ordered = False): 
	n_points = branch.GetNumberOfPoints()
	junc_id = branch.GetPointId(n_points - 1)

	if not ordered:
		dist = np.zeros(2)
		for i,pid in enumerate([0, n_points -1]):
			jid = branch.GetPointId(pid)
			jloc = poly_line.GetPoint(jid)
			dist[i] = vtk.vtkMath.Distance2BetweenPoints(jloc, root_loc)
		min_dist = np.argmin(dist)
		if min_dist == 0:
			junc_id = branch.GetPointId(0)

	junc_loc = poly_line.GetPoint(junc_id)
	junc_d = poly_line.GetPointData().GetScalars().GetValue(junc_id)


	min_sphere = np.zeros(n_points)
	
	for i in range(branch.GetNumberOfPoints()):
		pid = branch.GetPointId(i)
		if pid == junc_id:
			min_sphere[i] = 10000
			continue
		p_loc = poly_line.GetPoint(pid)
		p_d = poly_line.GetPointData().GetScalars().GetValue(pid)

		dist = m.sqrt(vtk.vtkMath.Distance2BetweenPoints(junc_loc, p_loc))
		rad_sum = (junc_d + p_d)/2
		min_sphere[i] = np.abs(dist - rad_sum)

	br_pid = np.argmin(min_sphere)

	pid = branch.GetPointId(br_pid)

	thick = poly_line.GetPointData().GetScalars().GetValue(pid)
	return pid, thick


def getTerminalThickness(poly_line,basic_info):
	end_points = basic_info['end_points']
	thickness = basic_info['medial_thickness']
	return {'dc': thickness[end_points], 'dc_loc':end_points}


def getEndBranchAngle(poly_line, basic_info, db_locs):

	end_branches = copy.deepcopy(basic_info['end_branches'])

	del basic_info
	angles = np.zeros(len(end_branches))
	norms = np.zeros((len(end_branches),2,3))
	juncs = np.zeros(len(end_branches), dtype = int)
	for i,br in enumerate(end_branches):

		angles[i],norms[i],juncs[i] = getBranchAngle(br, db_locs, poly_line)
	return {'angle':angles, 'angle_normals': norms, 'angle_loc': juncs}



def getBranchAngle(branch_id, sphere_ids, poly_line):

	branch = poly_line.GetCell(branch_id)
	n_points = branch.GetNumberOfPoints()
	junc_id = branch.GetPointId(n_points - 1)

	junc_loc = np.asarray(poly_line.GetPoint(junc_id))
	br_loc = np.asarray(poly_line.GetPoint(sphere_ids[branch_id]))

	br_normal = br_loc - junc_loc

	in_cells = vtk.vtkIdList()
	poly_line.GetPointCells(junc_id,in_cells)
	n_cells = in_cells.GetNumberOfIds()
	angles = np.zeros(n_cells)
	norms = np.zeros((n_cells,3))

	for i in range(n_cells):

		neigh_id = in_cells.GetId(i)
		neigh_loc = np.asarray(poly_line.GetPoint(sphere_ids[neigh_id]))
		neigh_normal = neigh_loc - junc_loc
		x = -1
		while np.all((neigh_normal == 0)):
			
			neigh_branch = poly_line.GetCell(neigh_id)
			sp_id = neigh_branch.GetPointId(int(neigh_branch.GetNumberOfPoints()/2)+x)
			neigh_loc = np.asarray(poly_line.GetPoint(sp_id))
			neigh_normal = neigh_loc - junc_loc
			x +=1

		norms[i] = neigh_normal
		angles[i] = calcAngleNormals(br_normal, neigh_normal)

	norms = norms[np.where(angles > 0.001)]
	angles = angles[np.where(angles > 0.0001)]
	norm_set = np.asarray([br_normal, norms[np.argmin(angles)]])

	return np.min(angles),norm_set,junc_id

def calcAngleNormals(normal1, normal2):

	normal_array = np.zeros((2,3))
	normal_array[0] = normal1
	normal_array[1] = normal2
	normal_array = normal_array ** 2
	normal_array = np.sum(normal_array, axis = 1)
	normal_array = np.sqrt(normal_array)

	dot = np.sum(normal1 * normal2)
	cos_angle = dot/(normal_array[0] * normal_array[1])

	if np.abs(cos_angle) <= 1:
		return m.acos(cos_angle)
	else:
		return m.acos(int(cos_angle))


def getMaximumThickness(coral_name, basic_info): 
	thickness = basic_info['medial_thickness']
	junction_ids = basic_info['junctions']
	return {'da': thickness[junction_ids], 'da_loc':junction_ids}








## Alternative of terminal thickness (see PhD Chris)
def getTerminalThicknessKrys(poly_line,basic_info):
	end_branches = basic_info['end_branches']
	n_end_br = len(end_branches)
	dcs = np.zeros(n_end_br)
	dc_loc = np.zeros(n_end_br, dtype = int)
	for i,br in enumerate(end_branches):
		branch = poly_line.GetCell(br)
		dc, dc_id = getTerminalThicknessBranch(poly_line, branch)
		dcs[i] = dc
		dc_loc[i] = dc_id
	return {'dc': dcs, 'dc_loc':dc_loc}



def getTerminalThicknessBranch(poly_line, branch):
	np = branch.GetNumberOfPoints()
	cur_thick = 0
	for i in range(np):
		pid = branch.GetPointId(i)
		thick = poly_line.GetPointData().GetScalars().GetValue(pid)
		if thick > cur_thick:
			cur_thick = thick
		else:

			return cur_thick,pid

	return cur_thick, pid
"""










