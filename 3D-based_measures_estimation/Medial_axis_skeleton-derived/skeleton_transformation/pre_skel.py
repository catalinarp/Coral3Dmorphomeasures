import vtk

def smoothPolyData(poly_data, max_iter = 100, passband = 0.005):
	smooth_filter = vtk.vtkWindowedSincPolyDataFilter()
	smooth_filter.SetInputData(poly_data)
	smooth_filter.SetNumberOfIterations(max_iter)
	smooth_filter.SetPassBand(passband)
	smooth_filter.SetBoundarySmoothing(True)
	smooth_filter.Update()
	return smooth_filter.GetOutput()
