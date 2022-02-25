from helpers import load_data

coral_name = '15Oki03'
polygon_mesh = load_data.readVTK(coral_name)
print(polygon_mesh.GetNumberOfPoints())