"""
surface_volume.py

This module is used to calculate measures directly related to the surface area
and volume of a polygon mesh (the sphericity and the S/V-ratio).
"""


import vtk
import math as m

def getSurfaceVolumeMeasures(poly_data):
	"""
	Get measures of a polygon mesh that are directly related to its surface 
	area and volume.

	Args:
		poly_data (vtk.vtkPolyData): polygon mesh 

	Returns:
		Dict: dictionary of float values with keys 
			- SA: surface area
			- V: volume
			- SV_ratio: surface-area-to-volume ratio
			- sphercity: sphericity (compactness)
	"""
	
	# calculate surface area and volume
	sa, v = getSurfaceVolume(poly_data)

	# S/V-ratio = SA/V
	sv_ratio = sa/v

	# calculate sphericity
	sphery = getSphericity(sa,v)

	return {'SA':sa, 'V':v, 'SV_ratio': sv_ratio, 'sphericity': sphery}

def getSurfaceVolume(poly_data):
	"""
	Get the surface area and volume of a polygon mesh using vtkMassProperties.

	Args:
		poly_data (vtkPolyData): polygon mesh 

	Returns:
		float: surface area
		float: volume
	"""

	# create filter for the mesh
	mass_filter = vtk.vtkMassProperties()
	mass_filter.SetInputData(poly_data)
	mass_filter.Update()

	# obtain surface area and volume of mesh
	surface = mass_filter.GetSurfaceArea()
	volume = mass_filter.GetVolume()
	return surface, volume


def getSphericity(sa, v):
	"""
	Calculates sphericity based on surface area and volume.

	Args:
		sa (float): surface area
		v (float): volume

	Returns:
		float: sphericity
	"""

	# sphericity = pi^(1/3) (6V)^(2/3) / SA
	sphere_sa_1 = m.pi**(1/3)
	sphere_sa_2 = (6 * v) ** (2/3)
	sphery = sphere_sa_1 * sphere_sa_2/sa
	return sphery