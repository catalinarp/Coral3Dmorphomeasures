import vtk
from helpers import local_directories as ldir

def readVTK(file_name, input_dir = ldir.POLYDATA):
	"""
	Reader reads a VTK file (.vtk) into vtkPolyData object.

	Args:
		file_name (str): filename of VTK file (without .vtk)
		input_dir (str): directory of VTK file
	
	Returns:
		vtk.vtkPolyData: a vtkPolyData object of the .vtk file 
	"""
	
	reader = vtk.vtkPolyDataReader()
	poly_data = readIntoVTK(reader, file_name, input_dir)
	return poly_data

def readIntoVTK(reader, file_name, input_dir, suffix = 'vtk'):
	reader.SetFileName(f"{input_dir}/{file_name}.{suffix}")
	reader.Update()
	return reader.GetOutput()

