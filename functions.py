import osmnx as ox
import shapely
import numpy as np
import subprocess
import concurrent.futures
import multiprocessing
import os
import math
import shapely.prepared
from stl import mesh
from mpl_toolkits import mplot3d
from matplotlib import pyplot

def generate_bbox_from_coords(lat_south,lon_west,lat_north,lon_east):
    #if given coordinates do not represent a perfect square we will extend the coordinates accordingly
    print(f"given coordinates:         {lat_south,lon_west,lat_north,lon_east}")
    #if coordinates are on the same side of the globe the distance is just the difference
    if (lat_south >= 0 and lat_north >= 0) or (lat_south <= 0 and lat_north <= 0):
        lat_distance = abs(lat_south - lat_north)
    #if coordinates are not on the same side of the globe we need to add the different values
    elif ((lat_south >= 0 and lat_north <= 0) or (lat_south <= 0 and lat_north >= 0)):
        lat_distance = abs(lat_south) + abs(lat_north)

    #if coordinates are on the same side of the globe the distance is just the difference
    if (lon_west >= 0 and lon_east >= 0) or (lon_west <= 0 and lon_east <= 0):
        lon_distance = abs(lon_west - lon_east)
    #if coordinates are not on the same side of the globe we need to add the different values
    elif ((lon_west >= 0 and lon_east <= 0) or (lon_west <= 0 and lon_east >= 0)):
        lon_distance = abs(lon_west) + abs(lon_east)
        #if the coordinates are more than 200 degree apart. We either got over the +/- 180 line or have invalid data anyways.
        #We take care about the +/- 180 line here
        if lon_distance > 200:
            lon_distance = abs(abs(lon_west)-180) + lon_east
    
    #Now we take care of the different distances
    if lat_distance > lon_distance:
        diff_distance = lat_distance - lon_distance
        if lon_west >= 0:
            lon_west += diff_distance/2
        elif lon_west <=0:
            lon_west -= diff_distance/2
        if lon_east >= 0:
            lon_east += diff_distance/2
        elif lon_east <=0:
            lon_east -= diff_distance/2
    elif lon_distance > lat_distance:
        diff_distance = lon_distance - lat_distance
        if lat_south >= 0:
            lat_south += diff_distance/2
        elif lat_south <=0:
            lat_south -= diff_distance/2
        if lat_north >= 0:
            lat_north += diff_distance/2
        elif lat_north <=0:
            lat_north -= diff_distance/2

    print(f"actually used coordinates: {lat_south,lon_west,lat_north,lon_east}")
    bbox = (lon_west, lat_south, lon_east, lat_north)
    return bbox


def fetch_location_data(bbox, location_type):
    # Fetch building footprints within the bounding box
    if location_type == "buildings":
        gdf = ox.features_from_bbox( bbox , tags = {'building': True, 'historic': ['citywalls']})
    if location_type == "paths":
        gdf = ox.features_from_bbox( bbox , tags = {'highway': True, 'man_made': ['pier'], 'railway': True, 'place': ['islet']})
    if location_type == "water":
        gdf = ox.features_from_bbox( bbox , tags = {'natural': ['water','reef'],'landuse': ['basin','salt_pond'], 'leisure': ['swimming_pool']})
    if location_type == "green":
        gdf = ox.features_from_bbox( bbox , tags = {'landuse': ['forest','meadow','grass','allotments','flowerbed','orchard','plant_nursery','vineyard','cemetery','recreation_ground','village_green'],'leisure': ['garden','park','pitch'],'natural': ['grassland','scrub','wood']})
    return gdf

def get_building_height(row, default_height=10):
    default_citywall_height = default_height*1.7
    # Check for various height attributes
    height_attrs = ['height', 'building:height', 'building:levels', 'historic:city_walls']
    for attr in height_attrs:
        if attr in row:
            height = row[attr]
            if isinstance(height, (int, float)) and not np.isnan(height):
                if attr == 'building:levels':
                    height = height * 3  # Assuming 3 meters per level
                    return height
                return height
            elif isinstance(height, str):
                try:
                    height_value = float(height.replace('m', '').strip())
                    if attr == 'building:levels':
                        height = height_value * 3  # Assuming 3 meters per level
                        return height
                    if attr == 'historic:citywalls':
                        return default_citywall_height
                    return height_value
                except ValueError:
                    continue
    return default_height

def create_solid_base(base_size, base_thickness=2, offset=0):
    # Define vertices for the base (solid block)
    base_vertices = [
        (offset + 0, offset + 0, 0),  # Bottom face
        (offset + base_size, offset + 0, 0),
        (offset + base_size, offset + base_size, 0),
        (offset + 0, offset + base_size, 0),
        (offset + 0, offset + 0, base_thickness),  # Top face (where buildings will sit)
        (offset + base_size, offset + 0, base_thickness),
        (offset + base_size, offset + base_size, base_thickness),
        (offset + 0, offset + base_size, base_thickness)
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

def create_triangle(vertices,side1,side2,side3):
    triangle_coords = ((vertices[side1][0],vertices[side1][1]),(vertices[side2][0],vertices[side2][1]),(vertices[side3][0],vertices[side3][1]))
    return shapely.Polygon(triangle_coords)

def create_planar_face(face_indicies, vertices, id=0):
    faces = []
    triangles_xy = []

    polygon_input = str(len(face_indicies))
    polygon_coords = []
    multiplicator = 10000
    for i in range(len(face_indicies)):
        polygon_input += os.linesep + str(int(vertices[face_indicies[i]][0]*multiplicator)) + " " + str(int(vertices[face_indicies[i]][1]*multiplicator))
    cwd = os.path.dirname(os.path.realpath(__file__))
    polygon_input_file = open("cache/polygon_input_" + str(id) + ".txt", "w")
    polygon_input_file.write(polygon_input)
    polygon_input_file.close()
    output = subprocess.check_output(["./a.out", "cache/polygon_input_" + str(id) + ".txt"], cwd=cwd, universal_newlines=True )
    output = output.replace("(","").replace(")","").replace(",","").split()
    for i in range(2):
        del output[0]
    for i in range(14):
        del output[len(output)-1]
    coord = []
    triangle = []
    for i in range(len(output)):
        if output[i].startswith("Unable"):
            print("")
            print("-------------------------------- WARNING --------------------------------")
            print("--------------- Unable to triangulate the following input ---------------")
            print(polygon_input)
            print("-------------------------------------------------------------------------")
            return faces
        output[i] = float(output[i])/multiplicator
        coord.append(output[i])
        if len(coord) == 2:
            triangle.append(coord)
            coord = []
        if len(triangle) == 3:
            triangles_xy.append(triangle)
            triangle = []

    tolerance = 0
    max_cycles = 10000
    step_size = 0.0000001
    for triangle in triangles_xy:
        sides = []
        for point in triangle:
            x,y = point
            cycles = 0
            hits = []
            while len(hits) != 1:
                hits = []
                for index in face_indicies:
                    if math.isclose(x, vertices[index][0], rel_tol=tolerance) and math.isclose(y, vertices[index][1], rel_tol=tolerance):
                        hits.append(index)
                if cycles > 200:
                    print(f"We are at cycle {cycles} of trying to find a valid triangle. Tolerance is {tolerance}.")
                    print(f"we are looking for x: {x} - y: {y}")
                    print(f"so far we found:")
                    for hit in hits:
                        print(f"x: {vertices[hit][0]} - y: {vertices[hit][1]}")
                if cycles > 0.9*max_cycles and len(hits) > 1:
                    #This will be fixed as soon as the external triangulation is done and we can use the proper ids for processing
                    print(f"after {cycles} cycles we still have more than one matching point. We will choose ony by luck.")
                    hits = ([hits[0]])
                if len(hits) < 1:
                    tolerance += step_size
                elif len(hits) > 1:
                    tolerance -= step_size
                if cycles == max_cycles:
                    raise Exception("Could not find a matching point! Try adjusting tolerance value.")
                cycles += 1
            sides.append(hits[0])
        if len(sides) < 3:
            raise Exception("Not enough Sides for the Triangle! Try adjusting tolerance value.")
        elif len(sides) > 3:
            raise Exception("More than 3 Sides in the Triangle! Try adjusting tolerance value.")
        faces.append(sides)
    return faces

def create_geometry(vertices,indicies):
    geometry_coords = []
    for index in indicies:
        x = vertices[index][0]
        y = vertices[index][1]
        geometry_coords.append((x,y))
    return shapely.Polygon(geometry_coords)

def create_vertices_list(exterior_coords, base_thickness, height):
    vertices_list = []
    for coord in exterior_coords:
        v_bottom = (coord[0], coord[1], base_thickness)
        v_top = (coord[0], coord[1], height + base_thickness)
        vertices_list.extend([v_bottom, v_top])
    return vertices_list

def create_side_faces(base_index, exterior_coords_len):
    side_faces = []
    for i in range(exterior_coords_len - 1):
        bottom1 = base_index + 2 * i
        bottom2 = base_index + 2 * (i + 1)
        top1 = base_index + 2 * i + 1
        top2 = base_index + 2 * (i + 1) + 1

        side_faces.append([bottom1, bottom2, top1])
        side_faces.append([top1, bottom2, top2])
    return side_faces

def scale_polygon(exterior_coords, bbox, target_size, base_size):
    #make the coords writable
    exterior_coords = list(exterior_coords)
    # Unpack the bounding box
    #north_lat, north_lng, south_lat, south_lng = bbox
    south_lng, south_lat, north_lng, north_lat = bbox
    # Calculate the scale factors for x and y dimensions
    lat_range = north_lat - south_lat
    lng_range = north_lng - south_lng
    # Calculate scaling factors based on the larger base
    scale_x = target_size / lng_range
    scale_y = target_size / lat_range
    # Calculate offsets to center the buildings on the enlarged base
    center_offset_x = (base_size - (scale_x * lng_range)) / 2
    center_offset_y = (base_size - (scale_y * lat_range)) / 2

    for i in range(len(exterior_coords)):
        exterior_coords[i] = list(exterior_coords[i])
        exterior_coords[i][0] = round(((exterior_coords[i][0] - south_lng) * scale_x) + center_offset_x, 3)
        exterior_coords[i][1] = round(((exterior_coords[i][1] - south_lat) * scale_y) + center_offset_y, 3)
    return shapely.Polygon(exterior_coords)

def cut_polygon(geometry):
    first_index = 0
    geometry_collection = []
    #we are doing this until all interiors are cut out and the original geometry is completely moved to the collection we will return
    interiors = len(list(geometry.interiors))
    while interiors > 0 and geometry.area != 0:
        #if we are dealing with a MultiPolygon due to previous cuts we are using points of the first polygon.
        #Since the MultiPolygon is reduced with every cut at some point it is no longer a MultiPolygon.
        if isinstance(geometry,shapely.MultiPolygon) or isinstance(geometry,shapely.GeometryCollection):
            for geom in geometry.geoms:
                if isinstance(geom,shapely.Polygon):
                    current_geometry_part = geom
                else:
                    if geom.area > 0:
                        print(f"We are not taking care of Geometry with Area {geom.area}. This is most likely a Bug.")
        else:
            current_geometry_part = geometry
        interiors = 0
        #if the first index to use to draw a line is higher than the the amount of indicies we are at the next polygon. So we need to reset it.
        if len(current_geometry_part.exterior.coords)-1 < first_index:
            first_index = 0
        #a line between first and middle vertex which should result in a more or less diagonal cut
        second_index = int(len(current_geometry_part.exterior.coords)/2)
        first_vertex = current_geometry_part.exterior.coords[first_index]
        second_vertex = current_geometry_part.exterior.coords[second_index]
        line = shapely.LineString([first_vertex, second_vertex])
        cutted_geoms = shapely.ops.split(current_geometry_part, line)
        #In some cases we don't cut anything, then we need another position.
        #Therefore we make the polygon more precise and move the starting index by 1
        for geom in cutted_geoms.geoms:
            if len(list(geom.interiors)) == 0 and isinstance(geom,shapely.Polygon):
                geometry_collection.append(geom)
                geometry = shapely.difference(geometry,geom)
            interiors += len(list(geom.interiors))
        if geometry.area > 0 and current_geometry_part.area > 0:
            max_segment_length = 1/(first_index+1)
            print(f"max_segment_length: {max_segment_length} - first_index: {first_index} - area: {current_geometry_part.area}")
            geometry = shapely.segmentize(geometry,max_segment_length)
            first_index += 1
            i = 0
            while (len(current_geometry_part.exterior.coords)-1) <= first_index:
                print(f"max_segment_length: {max_segment_length} - first_index: {first_index} - area: {current_geometry_part.area}")
                max_segment_length = current_geometry_part.length/(10+i)
                geometry = shapely.segmentize(geometry,max_segment_length)
                if isinstance(geometry,shapely.MultiPolygon) or isinstance(geometry,shapely.GeometryCollection):
                    current_geometry_part = geometry.geoms[0]
                else:
                    current_geometry_part = geometry
                if current_geometry_part.area == 0:
                    break
                i += 1
    return geometry_collection

def generate_object_list(gdf,default_height,max_height_mm):
    #Create a List of 3D Geometries (objects) out of the OSM Data
    #Geometries need to be preprocessed (Correct Type, No Holes, vertices in clockwise order, scaling, ...)
    
    object_list = []

    # Calculate the maximum object height
    max_building_height = gdf.apply(lambda row: get_building_height(row, default_height), axis=1).max()
    height_scale = max_height_mm / max_building_height

    for idx, row in gdf.iterrows():
        #We will have a geometry in the the first place [0] and the height in the second place [1]
        object = []
        #Check if we can process the object. Points and other stuff are not implemented (yet).
        if not (isinstance(row['geometry'], shapely.geometry.Polygon) or isinstance(row['geometry'], shapely.geometry.LineString)):
            print(f"Object of Type {idx[0]} with id {idx[1]} is not implemented (yet).")
            continue

        #Get the geometry out of the raw data
        object.append(row['geometry'])
        # If Object is a string convert to polygon
        if isinstance(object[0], shapely.geometry.LineString):
            #manual definition of path width depending on type
            if hasattr(row,"highway"):
                if row.highway == "motorway":
                    object[0] = shapely.buffer(object[0], 0.0001)
                elif row.highway == "trunk":
                    object[0] = shapely.buffer(object[0], 0.000075)
                elif row.highway == "primary":
                    object[0] = shapely.buffer(object[0], 0.00006)
                elif row.highway == "secondary":
                    object[0] = shapely.buffer(object[0], 0.000055)
                elif row.highway == "tertiary":
                    object[0] = shapely.buffer(object[0], 0.000045)
                elif row.highway == "residential":
                    object[0] = shapely.buffer(object[0], 0.000025)
                elif row.highway == "motorway_link":
                    object[0] = shapely.buffer(object[0], 0.00006)
                elif row.highway == "trunk_link":
                    object[0] = shapely.buffer(object[0], 0.00005)
                elif row.highway == "primary_link":
                    object[0] = shapely.buffer(object[0], 0.00005)
                elif row.highway == "secondary_link":
                    object[0] = shapely.buffer(object[0], 0.00004)
                elif row.highway == "tertiary_link":
                    object[0] = shapely.buffer(object[0], 0.00004)
                elif row.highway == "living_street":
                    object[0] = shapely.buffer(object[0], 0.000025)
                elif row.highway == "service":
                    object[0] = shapely.buffer(object[0], 0.000025)
                elif row.highway == "pedestrian":
                    object[0] = shapely.buffer(object[0], 0.000025)
                elif row.highway == "track":
                    object[0] = shapely.buffer(object[0], 0.000025)
                elif row.highway == "footway":
                    object[0] = shapely.buffer(object[0], 0.000025)
                elif row.highway == "bridleway":
                    object[0] = shapely.buffer(object[0], 0.000025)
                elif row.highway == "path":
                    object[0] = shapely.buffer(object[0], 0.000025)
                elif row.highway == "sidewalk":
                    object[0] = shapely.buffer(object[0], 0.000025)
                elif row.highway == "crossing":
                    object[0] = shapely.buffer(object[0], 0.000025)
                elif row.highway == "traffic_island":
                    object[0] = shapely.buffer(object[0], 0.000025)
                elif row.highway == "cycleway":
                    object[0] = shapely.buffer(object[0], 0.000025)
                else:
                    print(f"unclassified path width for type {row.highway}")
                    object[0] = shapely.buffer(object[0], 0.000025)
            elif hasattr(row,"man_made"):
                if row.man_made == "pier":
                    object[0] = shapely.buffer(object[0], 0.000004)
            elif hasattr(row,"railway") and isinstance(row.railway,str):
                if not row.tunnel == "yes":
                    object[0] = shapely.buffer(object[0], 0.00002)
            else:
                object[0] = shapely.buffer(object[0], 0.000025)
        
        #Get Object height
        object.append(get_building_height(row, default_height) * height_scale)

        object_list.append(object)
    return(object_list)
        
def preprocess_objects(object_list,bbox,target_size,scaling_factor):
    parameters = []
    base_size = target_size*scaling_factor
    for object in object_list:
        #Create a List with all parameters for multiprocessing
        parameters.append([object,bbox,target_size,base_size])
    #Call Preprocessing Function in Multiprocessing
    print("starting preprocessing with", multiprocessing.cpu_count(), "CPU Cores")
    with multiprocessing.Pool(multiprocessing.cpu_count()) as p:
        meta_object_list = p.map(cut_order_scale,parameters)
    ###################################
    #This is only for debugging without multiprocessing
    #meta_object_list = []
    #for param in parameters:
    #    meta_object_list.append(cut_order_scale(param))
    ######################################
    preprocessed_objects = []
    for meta_object in meta_object_list:
            for object in meta_object:
                preprocessed_objects.append(object)
    print(f"preprocessing done")
    return preprocessed_objects

def cut_order_scale(args):
    object,bbox,target_size,base_size = args
    processed_objects = []

    #If Polygon has holes, remove them by splitting it into multiple Polygons
    if(len(list(object[0].interiors))):
        geometry_list = cut_polygon(object[0])
    else:
        #nothing to cut
        geometry_list=([object[0]])

    for geometry in geometry_list:
        #check if points of polygon are clockwise ordered
        if shapely.algorithms.cga.signed_area(geometry.exterior) > 0:
            geometry = shapely.Polygon(reversed(geometry.exterior.coords))

        #scale object
        geometry = scale_polygon(geometry.exterior.coords,bbox,target_size,base_size)
        #simplify object
        geometry.simplify(0.1)
        #in case the geometry is invalid we try to fix it
        if not shapely.is_valid(geometry):
            valid_geom = shapely.make_valid(geometry)
            if isinstance(valid_geom, shapely.geometry.MultiPolygon):
                for geom in valid_geom.geoms:
                    if isinstance(geom,shapely.geometry.Polygon):
                        processed_objects.append([geom,object[1]])
            else:
                if isinstance(valid_geom,shapely.geometry.Polygon):
                    processed_objects.append([valid_geom,object[1]])
        else:
            processed_objects.append([geometry,object[1]])
    return processed_objects

def create_add_faces(base_index, exterior_coords,vertices,faces, id=0):
    # Create side faces
    faces.extend(create_side_faces(base_index, len(exterior_coords)))

    # Create top and bottom face
    top_face_indices = [base_index + 2 * i + 1 for i in range(len(exterior_coords) - 1)]
    faces.extend(create_planar_face(top_face_indices,vertices, id))

    # Create bottom face
    bottom_face_indices = [base_index + 2 * i for i in range(len(exterior_coords) - 1)]
    faces.extend(create_planar_face(bottom_face_indices,vertices, id))

    return faces

def prepare_3d_mesh(preprocessed_objects, target_size, scaling_factor, base_thickness=2, base_generation=True, object_generation=True):
    #Generation of 3D Mesh out of 2D Shapes with height
    vertices = []
    faces = []
    base_size = target_size * scaling_factor

    if base_generation == True:
        # Generate solid base
        base_vertices, base_faces = create_solid_base(base_size, base_thickness)
        # only half array will be used since we only need top OR bottom for the 2D object.
        base_geometry = create_geometry(base_vertices,range(int(len(base_vertices)/2)))
        vertices.extend(base_vertices)
        faces.extend(base_faces)
    else:
        # Generate solid base
        base_vertices, base_faces = create_solid_base(target_size, base_thickness, ((scaling_factor - 1) / 2) * target_size)
        # only half array will be used since we only need top OR bottom for the 2D object.
        base_geometry = create_geometry(base_vertices,range(int(len(base_vertices)/2)))

    
    if object_generation == True:
        id = 0
        for object in preprocessed_objects:
            print(f"processing object {id} of {len(preprocessed_objects)}")
            exterior_coords = list(object[0].exterior.coords)
            height = object[1]
            if id > 5000:
                pyplot.plot(*object[0].exterior.xy)
                pyplot.show()
            #Remove any overhangs over the base plate. This is required since some objects start within the bbox but end outside of it.
            #If this is the case we have object which are a lot to big and are not printable.
            #After Cleanup the Vertices of the final object are added to the whole list
            if shapely.intersects(base_geometry, object[0]):
                intersection = shapely.intersection(base_geometry, object[0])
                if isinstance(intersection, shapely.geometry.MultiPolygon):
                    for polygon in list(intersection.geoms):
                        exterior_coords = list(polygon.exterior.coords)
                        #Get Base Index (Before adding new vertices of this object)
                        base_index = len(vertices)
                        vertices.extend(create_vertices_list(exterior_coords, base_thickness, height))
                        faces = create_add_faces(base_index, exterior_coords, vertices, faces)
                else:
                    exterior_coords = list(intersection.exterior.coords)
                    #Get Base Index (Before adding new vertices of this object)
                    base_index = len(vertices)
                    vertices.extend(create_vertices_list(exterior_coords, base_thickness, height))
                    faces = create_add_faces(base_index, exterior_coords, vertices, faces)
            else:
                continue
            id += 1

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

def cut_two_categories(base_category,cutting_category,detect_islands=False):
    new_object_list = []
    for base_object in base_category:
        for cutting_object in cutting_category:
            if shapely.intersects(base_object[0],cutting_object[0]):
                #Island Detection Example:
                #When we remove the water area from the green area and no area is left, the green area IS fully surroundend and therefore is an island.
                if detect_islands and not shapely.difference(cutting_object[0],base_object[0]).area == 0:
                    continue
                base_object[0] = shapely.difference(base_object[0],cutting_object[0])
        if(base_object[0].area == 0):
            continue
        #after all cuts we append the object to the list
        if isinstance(base_object[0], shapely.geometry.MultiPolygon):
            for polygon in base_object[0].geoms:
                new_object_list.append(([polygon,base_object[1]]))
        if isinstance(base_object[0], shapely.geometry.Polygon):
            new_object_list.append(base_object)
    return new_object_list



def cut_all_categories(object_list_buildings,object_list_paths,object_list_water,object_list_greens):
    #If water is in a green area (like a river or a fointan) or a green is fully enclosed by water (an island) we need to cut the other parts. Otherwise fountains get lost in the end-result.
    #Since we need to implement this cutting algorythm anyways we also cut buildings and paths out of the other objects to prevent overlapping and have a nicer print result.
    #Order of objects from top to bottom:
    # buildings (They are high, so no cutting required)
    # paths
    # water (except the green is fully enclosed by water)
    # green

    #remove stuff from water areas
    object_list_water = cut_two_categories(object_list_water,object_list_paths)
    #This is a special case due to the island detection
    object_list_water = cut_two_categories(object_list_water,object_list_greens,detect_islands=True)

    #remove stuff from the green areas
    object_list_greens = cut_two_categories(object_list_greens,object_list_paths)
    object_list_greens = cut_two_categories(object_list_greens,object_list_water)

    return(object_list_buildings,object_list_paths,object_list_water,object_list_greens)