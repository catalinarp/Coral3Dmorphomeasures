import vtk
from helpers import local_directories as ldir
import pickle
import csv

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



def readMHD(file_name, input_dir = ldir.DIR_IMG):
	"""
	Reader reads metaimage into vtkImageData object

	Parameters
	----------
	file_name : str
		file_name of meta image header (without .mhd)

	Returns
	-------
	vtkImageData
		an vtkImageData object of the meta image
	"""
	reader = vtk.vtkMetaImageReader()
	image_data = readIntoVTK(reader,file_name, input_dir, 'mhd')
	return reader.GetOutput()

## writers

def writeFromVTK(writer,file_name, data,output_dir,suffix):
	writer.SetFileName(f"{output_dir}/{file_name}.{suffix}")
	writer.SetInputData(data)
	writer.Write()


def writeVTK(file_name, poly_data, output_dir = ldir.DIR_POLYDATA,suffix=""):
	"""
	Writer writes VTK-file (.vtk) of a vtkPolyData object
	
	Parameters
	----------
	file_name : str
		desired file_name and directory for file (without .vtk)
	poly_data : vtkPolyData
		object to save
	
	Returns
	-------
	None

	"""
	if suffix:
		file_name = f"{file_name}_{suffix}"
	writerVTK = vtk.vtkPolyDataWriter()
	writeFromVTK(writerVTK, file_name, poly_data, output_dir, 'vtk')

def writeMHD(file_name, image_data, output_dir = ldir.DIR_IMG):
	"""
	Writer writes meta image header (.MHD) and corresponding image (.zraw) of a
	vtkImageData object

	Parameters
	----------
	file_name : str
		file_name of .mhd and .zraw file (without .mhd)
	image_data : vtkImageData
		object to save

	Returns
	-------
	None
	"""

	writerMHD = vtk.vtkMetaImageWriter()
	writeFromVTK(writerMHD, file_name, image_data, output_dir, 'mhd')

def writePNG(file_name, image, output_dir = ldir.DIR_PHOTO):
	"""
	Writes picture (.png)

	Parameters
	----------
	file_name : str
		file_name of picture (without .png)
	image_data : vtkImage
		image to save

	Returns
	-------
	None
	"""
	writerPNG = vtk.vtkPNGWriter()
	print(file_name)
	writeFromVTK(writerPNG, file_name, image, output_dir, 'png')


#### CSV readers and writers
## readers
def readCSVAsDict(file_name, in_dir = ldir.DIR_RESULTS):
	"""
	Loads csv-file as dictionary, where column 'id' is the key and 
	the value is the other headers as key-value pairs
	
	Parameters
	----------
	file_name : str
		name of csv file (witout .csv)
	indir : str
		directory of the csv file
	
	Returns
	-------
	dict
		dictonary of all the data in csv file per id

	"""

	with open(f"{in_dir}/{file_name}.csv") as csvfile:
		reader = csv.DictReader(csvfile)
		out_dict = {}
		for row in reader:
			name = row['id']
			out_dict[name] = row

	return out_dict



def writeDictAsCSV(measure_dict, file_name, outdir = ldir.DIR_RESULTS):
	with open(f'{outdir}/{file_name}.csv', mode = 'w') as writefile:
		for key in measure_dict:
			headers = measure_dict[key].keys()
			break
		fieldnames = ['id']
		fieldnames.extend(headers)
		writer = csv.DictWriter(writefile, delimiter=',', fieldnames = fieldnames )
		writer.writeheader()
		for key in measure_dict:
			measures = measure_dict[key]
			measures['id'] = key
			writer.writerow(measures)


def writePickle(measure_dict, file_name, output_dir = ldir.DIR_RESULTS):
	with open(f'{output_dir}/{file_name}.pickle', 'wb') as handle:
    	 pickle.dump(measure_dict, handle)

def readPickle(file_name, input_dir = ldir.DIR_RESULTS):
	with open(f"{input_dir}/{file_name}.pickle", 'rb') as handle:
		data = pickle.load(handle)
	return data