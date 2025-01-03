import osmnx as ox
import shapely
import numpy as np
import subprocess
import os
import math
from stl import mesh
from mpl_toolkits import mplot3d
from matplotlib import pyplot
#C++ Tool from https://github.com/RikilG/Geometry-Algorithms/tree/master/Triangulation with this version https://github.com/RikilG/Geometry-Algorithms/commit/7bdf25e425b93dc6955331a48980a4b4d8051a6d

def fetch_location_data(bbox, location_type):
    # Fetch building footprints within the bounding box
    if location_type == "buildings":
        gdf = ox.features_from_bbox( bbox , tags = {'building': True})
    if location_type == "paths":
        gdf = ox.features_from_bbox( bbox , tags = {'highway': True})
    if location_type == "water":
        gdf = ox.features_from_bbox( bbox , tags = {'natural': ['water','reef'],'landuse': ['basin','salt_pond'],'leisure': ['swimming_pool']})
    if location_type == "green":
        gdf = ox.features_from_bbox( bbox , tags = {'landuse': ['forest','meadow','grass','allotments','flowerbed','orchard','plant_nursery','vineyard','cemetery','recreation_ground','village_green'],'leisure': ['garden','park','pitch'],'natural': ['grassland','scrub','wood']})
    return gdf

def get_building_height(row, default_height=10):
    # Check for various height attributes
    height_attrs = ['height', 'building:height', 'building:levels']
    for attr in height_attrs:
        if attr in row:
            height = row[attr]
            if isinstance(height, (int, float)) and not np.isnan(height):
                if attr == 'building:levels':
                    return height * 3  # Assuming 3 meters per level
                return height
            elif isinstance(height, str):
                try:
                    height_value = float(height.replace('m', '').strip())
                    return height_value
                except ValueError:
                    continue
    return default_height

def create_solid_base(base_size, base_thickness=2):
    # Define vertices for the base (solid block)
    base_vertices = [
        (0, 0, 0),  # Bottom face
        (base_size, 0, 0),
        (base_size, base_size, 0),
        (0, base_size, 0),
        (0, 0, base_thickness),  # Top face (where buildings will sit)
        (base_size, 0, base_thickness),
        (base_size, base_size, base_thickness),
        (0, base_size, base_thickness)
    ]

    # Define faces for the base
    base_faces = [
        [0, 1, 5], [0, 5, 4],  # Sides
        [1, 2, 6], [1, 6, 5],
        [2, 3, 7], [2, 7, 6],
        [3, 0, 4], [3, 4, 7],
        [4, 5, 6], [4, 6, 7],  # Top face
        [0, 1, 2], [0, 2, 3]   # Bottom face
    ]

    return base_vertices, base_faces

def check_if_overlapping_or_empty(geometry1, geometry2):
    if isinstance(geometry1, shapely.geometry.Polygon) and isinstance(geometry2, shapely.geometry.Polygon):
        difference = shapely.difference(geometry1, geometry2)
    else:
        return True
    if isinstance(difference, shapely.geometry.Polygon):
        if difference.area != 0:
            return True
        else:
            return False
    else:
        print(f"Error, Shapes invalid")

def create_triangle(vertices,side1,side2,side3):
    triangle_coords = ((vertices[side1][0],vertices[side1][1]),(vertices[side2][0],vertices[side2][1]),(vertices[side3][0],vertices[side3][1]))
    return shapely.Polygon(triangle_coords)

def create_planar_face(face_indicies, vertices, geometry_scaled):
    faces = []
    triangles_xy = []

    polygon_input = str(len(face_indicies))
    polygon_coords = []
    for i in range(len(face_indicies)):
        polygon_input += os.linesep + str(int(vertices[face_indicies[i]][0]*1000000000)) + " " + str(int(vertices[face_indicies[i]][1]*1000000000))
    cwd = os.path.dirname(os.path.realpath(__file__))
    polygon_input_file = open("cache/polygon_input.txt", "w")
    polygon_input_file.write(polygon_input)
    polygon_input_file.close()
    output = subprocess.check_output(["./a.out", "cache/polygon_input.txt"], cwd=cwd, universal_newlines=True )
    print(output)
    output = output.replace("(","").replace(")","").replace(",","").split()
    for i in range(2):
        del output[0]
    for i in range(14):
        del output[len(output)-1]
    coord = []
    triangle = []
    for i in range(len(output)):
        output[i] = float(output[i])/1000000000
        coord.append(output[i])
        if len(coord) == 2:
            triangle.append(coord)
            coord = []
        if len(triangle) == 3:
            triangles_xy.append(triangle)
            triangle = []

    tolerance = 1e-05
    for triangle in triangles_xy:
        sides = []
        print()
        for point in triangle:
            x,y = point
            for index in face_indicies:
                if math.isclose(x, vertices[index][0], rel_tol=tolerance) and math.isclose(y, vertices[index][1], rel_tol=tolerance):
                    sides.append(index)
        if len(sides) < 3:
            raise Exception("Not enough Sides for the Triangle! Try adjusting tolerance value.")
        if len(sides) > 3:
            raise Exception("More than 3 Sides in the Triangle! Try adjusting tolerance value.")
        faces.append(sides)

        triangle = create_triangle(vertices,sides[0],sides[1],sides[2])
        x,y = triangle.exterior.xy
        pyplot.plot(x,y)
    x,y = geometry_scaled.exterior.xy
    pyplot.plot(x,y)
    pyplot.show()
    return faces

def create_geometry(vertices,indicies):
    geometry_coords = []
    for index in indicies:
        x = vertices[index][0]
        y = vertices[index][1]
        geometry_coords.append((x,y))
    return shapely.Polygon(geometry_coords)

    
def prepare_mesh(gdf, bbox, target_size=180, max_height_mm=40, default_height=10, base_thickness=2, base_generation=True, object_generation=True):
    # Unpack the bounding box
    #north_lat, north_lng, south_lat, south_lng = bbox
    south_lng, south_lat, north_lng, north_lat = bbox


    # Calculate the scale factors for x and y dimensions
    lat_range = north_lat - south_lat
    lng_range = north_lng - south_lng

    # Define the base size as 20% larger than the target area
    base_size = target_size * 1.2    

    # Calculate scaling factors based on the larger base
    scale_x = target_size / lng_range
    scale_y = target_size / lat_range

    # Calculate offsets to center the buildings on the enlarged base
    center_offset_x = (base_size - (scale_x * lng_range)) / 2
    center_offset_y = (base_size - (scale_y * lat_range)) / 2

    vertices = []
    faces = []

    if base_generation == True:
        # Generate and append solid base
        base_vertices, base_faces = create_solid_base(base_size, base_thickness)
        vertices.extend(base_vertices)
        faces.extend(base_faces)

    count = 0 #debugging only
    if object_generation == True:
        # Calculate the maximum building height
        max_building_height = gdf.apply(lambda row: get_building_height(row, default_height), axis=1).max()
        height_scale = max_height_mm / max_building_height

        for idx, row in gdf.iterrows():
            print(f".", end="") #Some Progress Bar
            geometry = row['geometry']
            count += 1 #debugging only
            if count != 3 and count != 10: 
                continue
            if isinstance(geometry, shapely.geometry.Polygon) or isinstance(geometry, shapely.geometry.LineString):
                if isinstance(geometry, shapely.geometry.LineString):
                    geometry = shapely.buffer(geometry, 0.0001)
                if isinstance(geometry, shapely.geometry.Polygon):
                    exterior_coords = list(geometry.exterior.coords)
                # Create vertices for the geometry
                base_index = len(vertices)
                for coord in exterior_coords:
                    x = ((coord[0] - south_lng) * scale_x) + center_offset_x
                    #if x > base_size:
                    #    x = base_size
                    #if x < 0:
                    #    x = 0
                    y = ((coord[1] - south_lat) * scale_y) + center_offset_y
                    #if y > base_size:
                    #    y = base_size
                    #if y < 0:
                    #    y = 0
                    if y < 0 or x < 0:
                        print(f"X: {x} Y: {y}")
                    height = get_building_height(row, default_height) * height_scale
                    #print(f"Building at index {idx} with coordinates {exterior_coords} has height {height}")

                    v_bottom = (x, y, base_thickness)
                    v_top = (x, y, height + base_thickness)
                    vertices.extend([v_bottom, v_top])
                
                # Create side faces
                for i in range(len(exterior_coords) - 1):
                    bottom1 = base_index + 2 * i
                    bottom2 = base_index + 2 * (i + 1)
                    top1 = base_index + 2 * i + 1
                    top2 = base_index + 2 * (i + 1) + 1

                    faces.append([bottom1, bottom2, top1])
                    faces.append([top1, bottom2, top2])

                # Create top and bottom face
                top_face_indices = [base_index + 2 * i + 1 for i in range(len(exterior_coords) - 1)]
                geometry_scaled = create_geometry(vertices,top_face_indices)
                faces = faces + create_planar_face(top_face_indices,vertices,geometry_scaled)

                # Create bottom face
                bottom_face_indices = [base_index + 2 * i for i in range(len(exterior_coords) - 1)]
                geometry_scaled = create_geometry(vertices,bottom_face_indices)
                #faces = faces + create_planar_face(bottom_face_indices,vertices,geometry_scaled)

    vertices = np.array(vertices)
    faces = np.array(faces)

    return vertices, faces

def save_to_stl(vertices, faces, filename):
    mesh_data = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
    for i, face in enumerate(faces):
        for j in range(3):
            mesh_data.vectors[i][j] = vertices[face[j], :]

    # Create a new 3D plot
    #figure = pyplot.figure()
    #axes = figure.add_subplot(projection='3d')
    #axes.add_collection3d(mplot3d.art3d.Poly3DCollection(mesh_data.vectors))
    #scale = mesh_data.points.flatten()
    #axes.auto_scale_xyz(scale, scale, scale)
    #pyplot.show()

    mesh_data.save(filename)

def main():
    target_size = 180
    base_thickness = 2
    max_height_mm = 3
    default_building_height=40
    #bbox = (4.87123, 52.35893, 4.93389, 52.38351)  #Amsterdam
    bbox = (10.85891, 49.27478, 10.86771, 49.27973) #Suddersdorf
    #bbox = (10.863663, 49.277673, 10.864958, 49.278905) #Suddersdorf Weg Test
    #bbox = (11.06375, 49.44759, 11.09048, 49.45976) #Nürnberg Zentrum
    #bbox = min Longitude , min Latitude , max Longitude , max Latitude 

#    #Generation of Base Plate
#    vertices, faces = prepare_mesh(False, bbox, target_size=target_size, max_height_mm=max_height_mm, default_height=default_building_height, base_thickness=base_thickness, base_generation=True, object_generation=False)
#    save_to_stl(vertices, faces, 'export/standalone_base.stl')

    #Generation of Buildings
    gdf = fetch_location_data(bbox, "buildings")
    vertices, faces = prepare_mesh(gdf, bbox, target_size=target_size, max_height_mm=max_height_mm, default_height=default_building_height, base_thickness=base_thickness, base_generation=False, object_generation=True)
    save_to_stl(vertices, faces, 'export/buildings_without_base.stl')

#    #Generation of Paths
#    gdf = fetch_location_data(bbox, "paths")
#    vertices, faces = prepare_mesh(gdf, bbox, target_size=target_size, max_height_mm=max_height_mm*0.2, default_height=default_building_height, base_thickness=base_thickness, base_generation=False, object_generation=True)
#    save_to_stl(vertices, faces, 'export/paths_without_base.stl')
#
#    #Generation of Water
#    gdf = fetch_location_data(bbox, "water")
#    vertices, faces = prepare_mesh(gdf, bbox, target_size=target_size, max_height_mm=max_height_mm*0.2, default_height=default_building_height, base_thickness=base_thickness, base_generation=False, object_generation=True)
#    save_to_stl(vertices, faces, 'export/water_without_base.stl')
#
#    #Generation of "Green Areas" like Forest and Meadow
#    gdf = fetch_location_data(bbox, "green")
#    vertices, faces = prepare_mesh(gdf, bbox, target_size=target_size, max_height_mm=max_height_mm*0.2, default_height=default_building_height, base_thickness=base_thickness, base_generation=False, object_generation=True)
#    save_to_stl(vertices, faces, 'export/greens_without_base.stl')

if __name__ == "__main__":
    main()
