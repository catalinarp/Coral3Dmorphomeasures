from helpers import load_data, local_directories 
from Medial_axis_skeleton_based import skeleton_distances

coral_name = '15Oki03'
polygon_skel = load_data.readVTK(coral_name, local_directories.LINE)

skeleton_distances.getSkeletonDistances(polygon_skel)

