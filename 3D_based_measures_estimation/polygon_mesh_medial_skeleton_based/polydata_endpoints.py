import vtk
import numpy as np 
import vtk.util.numpy_support as nps 
import copy
from .basic_skeleton_measures import getBasicSkelMeasures

def setScalarsToIds(poly_data):
	n_points = poly_data.GetNumberOfPoints()
	x = np.arange(n_points)
	poly_data.GetPointData().SetScalars(nps.numpy_to_vtk(x))


def reverseEndBranchesTemp(poly_line, end_branches, end_points):
	for i,br in enumerate(end_branches):

		branch = poly_line.GetCell(br)
	
		if branch.GetPointId(0) != end_points[i]:
			poly_line.ReverseCell(br)

def getTerminalEndpointTemp(poly_line, branch, ep_id):
	np = branch.GetNumberOfPoints()
	cur_thick = 0
	for i in range(np):
		pid = branch.GetPointId(i)
		thick = poly_line.GetPointData().GetScalars().GetValue(pid)
		if thick > cur_thick:
			cur_thick = thick
		else:
			return pid
	return pid

def getBranchNormalDirection(poly_line, branch):
	n_points = branch.GetNumberOfPoints()
	first = branch.GetPointId(0)
	last = branch.GetPointId(int(n_points/2))
	table = np.zeros((2,3))
	
	table[0] = poly_line.GetPoint(last)
	table[1] = poly_line.GetPoint(first)
	norm = table[1] - table[0]
	return norm

def clip_polydata(poly_data, cent, norm,thick):
	
	poly_copy = vtk.vtkPolyData()
	poly_copy.DeepCopy(poly_data)


	cyl = vtk.vtkCylinder()
	cyl.SetCenter(cent)

	cyl.SetAxis(norm)

	cyl.SetRadius(thick/2)
	
	plane = vtk.vtkPlane()
	plane.SetOrigin(cent)
	plane.SetNormal(norm)

	clipper = vtk.vtkExtractPolyDataGeometry()
	clipper.SetInputData(poly_copy)
	clipper.SetImplicitFunction(cyl)
	clipper.ExtractInsideOn()
	
	clipper.Update()
	poly_copy = clipper.GetOutput()
	

	
	clipper = vtk.vtkExtractPolyDataGeometry()
	clipper.SetInputData(poly_copy)

	clipper.SetImplicitFunction(plane)
	clipper.ExtractInsideOff()
	clipper.Update()
	poly_copy = clipper.GetOutput()

	connect = vtk.vtkPolyDataConnectivityFilter()
	connect.SetInputData(poly_copy)
	connect.SetClosestPoint(cent)
	connect.SetExtractionModeToClosestPointRegion()
	connect.Update()
	poly_copy = connect.GetOutput()
	try:
		ids = nps.vtk_to_numpy(poly_copy.GetPointData().GetScalars())
	except:
		ids = []

	return ids

def selectGoodPoints(endpoints, n_points):
	mean = np.mean(n_points)
	std = np.std(n_points)
	up_lim = mean + 3 *std
	final_points = []
	for i in range(len(n_points)):
		if n_points[i] < up_lim:
			final_points.extend(endpoints[i])
	final_points = np.asarray(list(set(final_points)), dtype =int)
	return final_points




def getEndPointsPolyData(poly_data, poly_line):
	
	basic_dict = getBasicSkelMeasures(poly_line)
	ep_ids = basic_dict['end_points']
	eb_ids = basic_dict['end_branches']
	max_thick = basic_dict['max_thick']
	# prepare polydata
	poly_ids = []

	n_points = np.zeros(len(ep_ids))
	for i in range(len(ep_ids)):

		bid = eb_ids[i]
		branch = poly_line.GetCell(bid)
		end_point_loc = poly_line.GetPoint(ep_ids[i])
		norm = getBranchNormalDirection(poly_line, branch)
		thick = max_thick[bid]
		ids = clip_polydata(poly_data, end_point_loc, norm, thick)

		n_points[i] = len(ids)
		poly_ids.append(copy.deepcopy(ids))

	final_points = selectGoodPoints(poly_ids, n_points)
	return {'point_ids':final_points}