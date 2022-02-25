import vtk
import  vtk.util.numpy_support as nps 
import numpy as np

class SkelVoxToGraph(vtk.vtkProgrammableFilter):
	"""
	This filter transforms a voxel skeleton to a skeleton graph. 
	"""
	def __init__(self) :
		"""
		inialize filter, set excecute method to __MakeGraph and required input and output data types
		"""

		# set input and output data types
		self.GetInputPortInformation(0).Set(vtk.vtkAlgorithm.INPUT_REQUIRED_DATA_TYPE(), 'vtkImageData')
		self.GetOutputPortInformation(0).Set(vtk.vtkDataObject.DATA_TYPE_NAME(), "vtkPolyData")
		# set execution method to __MakeGraph
		self.SetExecuteMethod(self.__Makegraph)
		self.offs6 = []
		self.offs18 = []
		self.offs18_6 = []
		self.offs26 = []
		self.offs26_18 = []
	
	def __Makegraph(self) :
		# get input and output, make output data empty
		input = self.GetInput()
		output = self.GetOutput()
		output.Reset()
		self.offs6 = []
		self.offs18 = []
		self.offs18_6 = []
		self.offs26 = []
		self.offs26_18 = []
		self.gen_offsets(input)

		tp = vtk.vtkThresholdPoints()
		tp.SetInputData(input)
		tp.ThresholdByUpper(-.5)
		tp.Update()
		# set the points of the output to be coordinates of the pixels that have a value of minimally 1
		# (that is probably color)
		output.SetPoints(tp.GetOutput().GetPoints())

		tp.RemoveAllInputs()
		del tp

		# set all the scalars of the point data to be equal to the original id of the points
		ids = vtk.vtkIdTypeArray()
		thicks = vtk.vtkFloatArray()
		ids.SetName('OriginalIds')
		numPts = output.GetNumberOfPoints()
		scal = input.GetPointData().GetScalars()
		for i in range(numPts) :
			id = input.FindPoint(output.GetPoint(i))
			ids.InsertNextValue(id)
			thickness = scal.GetComponent(id,0)
			thicks.InsertNextValue(thickness)
		thicks.SetName('MedialThickness')

		output.GetPointData().SetScalars(ids)
		del ids

		# build the skeleton graph
		self.buildgraph(input, output,thicks)


	def nbors(self, id, offs, scal, processed = []) :
		good_nbors = []
		for offset in offs :
			nbor_id = id+offset

			if (scal.GetComponent(nbor_id, 0) >= 0) :
				good_nbors += [nbor_id]

		return good_nbors


	def existspath(self, source, dest, connmap) :

		for cur_neigh in connmap[source]:
			if dest in connmap[cur_neigh]:
				return True

		return False

	def buildgraph(self, skelimage, skelpd,thicks) :
		
		# get the increments and the dimensions of the input data
		inc = skelimage.GetIncrements()
		dim = skelimage.GetDimensions()
		ca = vtk.vtkCellArray()
		
		# obtian all the scalars of the input data
		scal = skelimage.GetPointData().GetScalars()
		# array of pdpointid -> imgpointid info

		# get all the scalars of the output data
		pointdata = skelpd.GetPointData()
		skelptids = pointdata.GetScalars()
		points = nps.vtk_to_numpy(skelpd.GetPoints().GetData())

		vtk_points = vtk.vtkPoints()
		vtk_points.SetData(nps.numpy_to_vtk(points))
		skelpd.SetPoints(vtk_points)

		pd2img = {}
		img2pd = {}
		connto = {}

		# make a dictionaries to translate between the ids of the input and those of the output
		for id in range(skelptids.GetNumberOfTuples()) :
			pd2img[id] = skelptids.GetValue(id)
			img2pd[skelptids.GetValue(id)] = id
			connto[id] = []


		numtup = skelptids.GetNumberOfTuples()
		for offs in self.offs6, self.offs18_6, self.offs26_18 :
			
			for id in range(numtup) :
				# obtain all the neighbours of the cell
				nborlist = [img2pd[x] for x in self.nbors(pd2img[id], offs, scal)]
				for nb in nborlist:
					if not self.existspath(id, nb, connto) :
						cl = connto[id]
						if nb not in cl:
							cl.append(nb)
							cl.sort()
						cl = connto[nb]
						if id not in cl:
							cl.append(id)
							cl.sort()


		
		nonreg = [id for id,conns in connto.items() if len(conns) != 2]
		finalpts = []
		branches = []
		startpoints = []
		for startpt in nonreg :
			for branchstart in connto[startpt] :
				if branchstart not in finalpts and branchstart not in startpoints:

					branch = [startpt, branchstart]
					curpt = branchstart
					while len(connto[curpt]) == 2:
						curpt = connto[curpt]
						if curpt[0] in branch:
							curpt = curpt[1]
						else :
							curpt = curpt[0]
						if curpt in branch and curpt != branch[0] :
							break
						branch.append(curpt)
					finalpts.append(branch[-2])
					branches.append(branch)
			startpoints.append(startpt)

		outca = vtk.vtkCellArray()
		for branch in branches :
			outca.InsertNextCell(len(branch))
			for point in branch :
				outca.InsertCellPoint(point)

		pointdata.SetScalars(thicks)
		skelpd.SetLines(outca)
		del ca
		del outca

	def gen_offsets(self, img) :
		"""
		Generate offsets for locating neighbors in image.
		"""

		# get increments and dimensions
		inc = img.GetIncrements()
		dim = img.GetDimensions()


		for x in range(-1, 2) :
			for y in range(-1, 2) :
				for z in range(-1, 2) :
					pt0 = [0,0,0]
					pt1 = [x,y,z]
					if (pt0 == pt1) :
						continue
					
					# calculate the squared distance between [0,0,0] and values
					dist = vtk.vtkMath.Distance2BetweenPoints(pt0, pt1)
					offset = x*inc[0] + y*inc[1] + z*inc[2]
					if (dist == 1) :
						self.offs6 += [offset]
					if (dist <= 2) :
						self.offs18 += [offset]
					if (dist == 2) :
						self.offs18_6 += [offset]
					if (dist == 3) :
						self.offs26_18 += [offset]
					if (dist <= 3):
						self.offs26 += [offset]