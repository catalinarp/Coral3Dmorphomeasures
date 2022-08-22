def getEndCurvature(curvature, end_points, measure = 'gauss'): 
	return curvature[measure][end_points['point_ids']]