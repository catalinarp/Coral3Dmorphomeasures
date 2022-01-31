# <b>3D-based measures estimation</b>

This folder contains the files and scripts used to estimate 3D-based measures from either triangulated polygon meshes, medial axis skeleton graphs or a combination of both.

### Polygon mesh-based measures
The resulting triangulated polygon meshes were analysed and visualized using the Visualization Toolkit v9.1.0 (VTK; Schroeder et al., 2006) in Python v3.8 (Van Rossum and Drake, 2009) and the Insight Toolkit v5.1.2 (ITK; Ibáñez et al., 2003) in c++20 (ISO/IEC, 2020). 

### Medial axis skeleton-derived measures
To capture the topological branching structure of corals and facilitate the estimation of measures related to this type of morphology, the medial axis skeleton was extracted from the previously rendered 3D models using a voxel thinning algorithm (Lee et al., 1994; Homann, 2007).

### Polygon mesh and medial axis skeleton graph-based measures
The medial skeleton axis and the smoothed polygon meshes were obtained for each specimen and used to estimate additional measures form the coral specimens.
<br>
<br>


![Figure_S2](https://user-images.githubusercontent.com/11543453/151872422-cfbde3a8-de5d-4eae-b550-f73368691e3e.png)
<br>
<br>
<b>Workflow chart | Summary of the approximations used for 3D-based measures estimation</b><br>
Different types of measures were obtained extracting information from different shape representations. 
