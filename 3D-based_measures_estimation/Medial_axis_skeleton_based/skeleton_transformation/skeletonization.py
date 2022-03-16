
def makeVoxelSkeleton(in_dir, filename):
	
	skel_process = subprocess.Popen(["./skel_itk",in_dir,filename])
	skel_process.communicate()
	voxel_skel = readMHD(f"{filename}_skel", in_dir)
	
	return voxel_skel

def skeletonize(poly_data, spacing = .05, skel_dir = vtk.DIR_TEMP):
	"""
	
	"""
	
	# voxelize smooth polydata
	vox_img = voxelizePolyData(poly_data, spacing = spacing)
	fn = 'vox'
	writeMHD(fn, vox_img, skel_dir)
	# filename = saveVoxelize(poly_data, 'tmp', out_dir = skel_dir, spacing = spacing)

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

def smoothPolyData(poly_data, max_iter = 100, passband = 0.005):
	smooth_filter = vtk.vtkWindowedSincPolyDataFilter()
	smooth_filter.SetInputData(poly_data)
	smooth_filter.SetNumberOfIterations(max_iter)
	smooth_filter.SetPassBand(passband)
	smooth_filter.SetBoundarySmoothing(True)
	smooth_filter.Update()
	return smooth_filter.GetOutput()


def voxelizePolyData(poly_data, spacing = .05, tolerance=0):
	"""
	Transforms a polygon mesh into a voxel image. This code 

	Args:
		poly_data (vtk.vtkPolyData): closed polygon mesh
		spacing (float
	"""
	if isinstance(spacing,float) or isinstance(spacing, int):
		spacing = [spacing] * 3
	
	# Get some basic info of poly data
	poly_bounds = np.asarray(poly_data.GetBounds())
	min_bounds = poly_bounds[::2]
	max_bounds = poly_bounds[1::2]

	# Prepare image to cut in
	white_img = vtk.vtkImageData()
	white_img.SetSpacing(spacing)

	dims = np.asarray([m.ceil((max_bounds[i]-min_bounds[i])/spacing[i]) + 1 for i in range(3)])
	white_img.SetDimensions(dims)
	ext = dims -1
	white_img.SetExtent(0, ext[0], 0, ext[1], 0, ext[2])
	white_img.SetOrigin(min_bounds)
	white_img.AllocateScalars(vtk.VTK_UNSIGNED_CHAR, 1)

	in_val = 255
	out_val = 0
	n_pix = white_img.GetNumberOfPoints()
	in_val_array = np.ones(n_pix) * in_val
	white_img.GetPointData().SetScalars(nps.numpy_to_vtk(in_val_array))
	
	# Make an image stencil of the polydata
	stencil_filter = vtk.vtkPolyDataToImageStencil()
	stencil_filter.SetTolerance(tolerance)
	stencil_filter.SetInputData(poly_data)
	stencil_filter.SetOutputOrigin(min_bounds)
	stencil_filter.SetOutputSpacing(spacing)
	stencil_filter.SetOutputWholeExtent(white_img.GetExtent())
	stencil_filter.Update()

	# cut the image
	img_stencil = vtk.vtkImageStencil()
	img_stencil.SetInputData(white_img)
	img_stencil.SetStencilConnection(stencil_filter.GetOutputPort())
	img_stencil.ReverseStencilOff()
	img_stencil.SetBackgroundValue(out_val)
	img_stencil.Update()

	return img_stencil.GetOutput()