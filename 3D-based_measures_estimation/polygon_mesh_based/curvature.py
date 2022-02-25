import vtk
import numpy as np
import math as m
import vtk.util.numpy_support as nps

def getCurvature(poly_data):
	"""
	Calculates the Gaussian (K), Mean (H), minimum curvature (k1) and maximum
	curvature for each vertex of a polygon mesh. This is a VTK based 
	implementation of the discrete scheme proposed in [1]. 

	Each vertex is associated with an area (Amixed, fig. 4 of [1]). K curvature 
	is obtained using the angle deficit method (eq. 9 of [1]). H is obtained by
	first calculating the Mean Curvature Normal Operator (MCNO, eq. 8 of [1]), 
	where H=||MCNO||/2. The sign of H is obtained by comparing its direction
	with the normal vector. k1 and k2 are obtained by as in section 5.1 of [1].

	In this implentation curvature is obtained by looping over the faces rather
	than over the edges as this seemed more suitable for the way polygon meshes
	are structured within VTK.

	[1] Meyer, Mark and Desbrun, Mathieu and Schroder, Peter and Barr, Alan H.
	Discrete differential-geometry operators for triangulated 2-manifolds. 
	Visualization and mathematics III, 2003 35-57.

	Args:
		poly_data (vtk.vtkPolyData): polygon mesh
	
	Returns:
		Dict: dictionary of numpy arrays (float) with keys
			- Gauss: K of each vertex (ordered by vertex ID)
			- Mean: H of each vertex (ordered by vertex ID)
			- Minimum: k2 of each vertex (ordered by vertex ID)
			- Maximum: k1 of each vertex (ordered by vertex ID)
			- Areas: Amixed of each vertex (ordered by vertex ID)

	"""
	n_points = poly_data.GetNumberOfPoints()
	n_face = poly_data.GetNumberOfCells()

	areas = np.zeros(n_points)

	# angle deficit = 2pi - sum(adjecent angles)
	angles = np.ones(n_points) * 2 * m.pi
	defic = np.zeros((n_points,3))
	angles_tri = np.zeros(3)

	areas_tri = np.zeros(3)
	points_tri = np.zeros((3,3))
	length_sq = np.zeros(3)

	# loop over every face
	for face_id in range(n_face):
		
		if not face_id % 100000:
			print(face_id)
		
		cell = poly_data.GetCell(face_id)
		points = cell.GetPoints()
		
		# get the coordinates of the vertices
		for i in range(3):
			points_tri[i] = points.GetPoint(i)

		# get the vectors of the face
		vecs = points_tri - points_tri[[1,2,0]]
		
		# get the three angles and squared side lengths
		for i in range(3):
			length_sq[i] = vtk.vtkMath.Distance2BetweenPoints(points_tri[i], 
				points_tri[i-2])
			angles_tri[i] = vtk.vtkMath.AngleBetweenVectors(-vecs[i], vecs[i-1])

		# check if angles are obtuse
		obtuse = (angles_tri > m.pi/2)

		# calculate cotanges of every angle
		cots = np.cos(angles_tri)/np.sin(angles_tri)
		
		# if there is an obtuse angle
		if np.sum(obtuse):
			area = cell.ComputeArea()

			# divide area of face to each vertex in fixed proportions (*8)
			for i in range(3):
				if obtuse[i]:
					areas_tri[i] = area*4
				else:
					areas_tri[i] = area*2
		else:
			# if no obtuse area use the voronoi area to divide
			for i in range(3):
				areas_tri[i] = cots[i-1] * length_sq[i] + cots[i-2] * length_sq[i-1]
			
		for i in range(3):
			pid = cell.GetPointId(i)
			defic[pid] += cots[i-1] * vecs[i] + cots[i-2] * -vecs[i-1]
			areas[pid] += areas_tri[i]
			angles[pid] -= angles_tri[i]


	# K = angle deficit devided by associated vertex area
	gauss = 8*angles/areas
	MCNO = (4*defic.T/areas.T).T
	mean_c = np.linalg.norm(MCNO, axis = 1)/2
	areas /= 8 
	del angles, defic


	# get outward normal vector of each vertex
	poly_normals = vtk.vtkTriangleMeshPointNormals()
	poly_normals.SetInputData(poly_data)
	poly_normals.Update()
	poly_data = poly_normals.GetOutput()
	del poly_normals
	normals = poly_data.GetPointData().GetArray("Normals")

	for n in range(n_points):
		# if H = 0, do not compare to norm
		if not mean_c[n]:
			continue

		norm = normals.GetTuple(n)
		
		# normalize the MCNO
		meandir = 2*MCNO[n]/mean_c[n]

		# change the sign of the MCNO if the general direction is the same as 
		# the normal vector 
		if np.sum(norm * meandir) < 0:
			mean_c[n] *= -1
	
	del normals, poly_data	
	
	kprep = mean_c * mean_c - gauss
	
	# calculate k1 and k2
	k_min = np.zeros(n_points)
	k_max = np.zeros(n_points)
	for i in range(n_points):
		kp = kprep[i]
		h = mean_c[i]
		
		# if no real solution exists using the standard formula take H
		if kp < 0:
			k_min[i] = h
			k_max[i] = h
		# otherwise take H +- sqrt(K^2-H)
		else:		
			k_min[i] = h - m.sqrt(kp)
			k_max[i] = h + m.sqrt(kp)

	return {'Gauss':gauss, 'Mean': mean_c, 'Minimum': k_min, 'Maximum': k_max, \
		'areas':areas}

