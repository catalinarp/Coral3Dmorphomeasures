import vtk 
from pre_skel import *

def makeVoxelSkeleton(in_dir, filename):
	
	skel_process = subprocess.Popen([".skel_itk",in_dir,filename])
	skel_process.communicate()
	voxel_skel = readMHD(f"{filename}_skel", in_dir)
	
	return voxel_skel

def skeletonize(poly_data, spacing = .05, skel_dir = vtk.DIR_TEMP):
	"""
	
	"""
	#smooth poly_data
	smoothPolyData(poly_data)

	# voxelize smooth polydata
	vox_img = voxelizePolyData(poly_data, spacing = spacing)
	fn = 'vox'
	writeMHD(fn, vox_img, skel_dir)

	print('voxelized')
	time.sleep(10)
	del vox_img
	# extract voxel skeleton
	voxel_skel = makeVoxelSkeleton(skel_dir, fn)

	# make a VTK line object of voxel skeleton
	sg = SkelGraph()
	sg.SetInputData(voxel_skel)
	
	sg.Update()
	return sg.GetOutput()

