import osmnx as ox
import shapely
import numpy as np
import subprocess
import concurrent.futures
import multiprocessing
import os
import math
import pandas
import geopandas
from stl import mesh
from mpl_toolkits import mplot3d
from matplotlib import pyplot


def square_bbox_from_center_point(lat, lon, distance):
    gs = geopandas.GeoSeries(shapely.Point(lon, lat))
    gdf = geopandas.GeoDataFrame(geometry=gs,crs='EPSG:4326')
    gdf = gdf.to_crs('EPSG:3857')
    res = gdf.buffer(
        distance=distance,
        cap_style=3,
    )
    returnpoly = res.to_crs('EPSG:4326').iloc[0]
    return returnpoly.bounds

def square_bbox_from_vertices(min_lat, min_lon, max_lat, max_lon):
    gs = geopandas.GeoSeries(shapely.Polygon(((min_lon, min_lat),(max_lon, min_lat), (max_lon, max_lat), (min_lon, max_lat))))
    gdf = geopandas.GeoDataFrame(geometry=gs,crs='EPSG:4326')
    bounds = gdf.to_crs('EPSG:3857').bounds
    
    side_length_x = bounds.maxx[0] - bounds.minx[0]
    side_length_y = bounds.maxy[0] - bounds.miny[0]
    side_length_avg = (side_length_x + side_length_y) / 2
    center_x = bounds.minx[0] + ( side_length_x / 2 )
    center_y = bounds.miny[0] + ( side_length_y / 2 )

    gs = geopandas.GeoSeries(shapely.Point(center_x, center_y))
    gdf = geopandas.GeoDataFrame(geometry=gs,crs='EPSG:3857')
    res = gdf.buffer(
        distance=side_length_avg,
        cap_style=3,
    )

    returnpoly = res.to_crs('EPSG:4326').iloc[0]
    return returnpoly.bounds

def fetch_location_data(bbox, location_type):
    # Fetch building footprints within the bounding box
    try:
        if location_type == "buildings":
            gdf = ox.features_from_bbox( bbox , tags = {'building': True, 'historic': ['citywalls']})
        if location_type == "paths":
            gdf = ox.features_from_bbox( bbox , tags = {'highway': True, 'man_made': ['pier'], 'railway': True})
        if location_type == "water":
            gdf = ox.features_from_bbox( bbox , tags = {'natural': ['water','reef'],'landuse': ['basin','salt_pond'], 'leisure': ['swimming_pool']})
        if location_type == "green":
            gdf = ox.features_from_bbox( bbox , tags = {'landuse': ['forest','meadow','grass','allotments','flowerbed','orchard','plant_nursery','vineyard','cemetery','recreation_ground','village_green'],'leisure': ['garden','park','pitch'],'natural': ['grassland','scrub','wood']})
    except: # Return Empty Polygon if location type is not found
        gdf = ([[shapely.Polygon(),0]])
    return gdf

def get_object_height(row, default_height=10):
    #Areas like water and grass need to be on same level as base plate
    #Areas like forest or scrub need to be elevated slightly
    #Paths are also even with base plate
    #TODO: What about bridges?!
    special_area_attributes = ['natural', 'landuse', 'leisure', 'highway', 'man_made', 'railway']
    for attr in special_area_attributes:
        if attr in row:
            test = row[attr]
            if attr == 'natural' and (row[attr] == 'water' or row[attr] == 'reef' or row[attr] == 'grassland' or row[attr] == 'scrub' or row[attr] == 'wood'):
                #Stuff that needs to be the same level as the base plate
                return -1
            elif attr == 'landuse' and (row[attr] == 'basin' or row[attr] == 'salt_pond' or row[attr] == 'grass' or row[attr] == 'allotments' or row[attr] == 'flowerbed' or row[attr] == 'village_green'):
                return -1
            elif attr == 'natural' and (row[attr] == 'scrub'):
                #scrubs are between 1 and 2 meters high
                return 1
            elif attr == 'natural' and (row[attr] == 'wood'):
                #woods probably about 25 meters high
                return 15
            elif attr == 'landuse' and (row[attr] == 'forest'):
                return 15
            elif attr == 'landuse' and (row[attr] == 'orchard'):
                return 4
            elif attr == 'landuse' and (row[attr] == 'plant_nursery' or row[attr] == 'vineyard' or row[attr] == 'recreation_ground'):
                return 2
            elif attr == 'landuse' and (row[attr] == 'meadow' or row[attr] == 'cemetery'):
                return 0.5
            elif attr == 'leisure' and (row[attr] == 'garden' or row[attr] == 'park'):
                return -1
            elif attr == 'leisure' and (row[attr] == 'pitch'):
                return -1
            elif attr == 'highway' and not pandas.isna(row[attr]):
                #Stuff that needs to be the same level as the base plate
                return -1
            elif attr == 'man_made' and (row[attr] == 'pier'):
                #piers are maybe around 2 meter high?
                return 2
            elif attr == 'railway' and not pandas.isna(row[attr]):
                return 0.5

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

def create_vertices_list(exterior_coords, base_thickness, height, height_offset=1):
    vertices_list = []
    for coord in exterior_coords:
        v_bottom = (coord[0], coord[1], base_thickness + height_offset)
        v_top = (coord[0], coord[1], height + base_thickness + height_offset)
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
    geoms_with_interiors = [geometry]
    geoms_without_interiors = []
    #we are doing this until all interiors are cut out and the original geometry is completely moved to the collection we will return
    cycles = 0
    while len(geoms_with_interiors) > 0:
        geom_list_copy = geoms_with_interiors
        geoms_with_interiors = []
        for geom in geom_list_copy:
            min_x = geom.exterior.xy[0][0]
            max_x = geom.exterior.xy[0][0]
            for x in geom.exterior.xy[0]:
                if x < min_x:
                    min_x = x
                if x > max_x:
                    max_x = x
            min_y = geom.exterior.xy[1][0]
            max_y = geom.exterior.xy[1][0]
            for y in geom.exterior.xy[1]:
                if y < min_y:
                    min_y = y
                if y > max_y:
                    max_y = y
            if cycles % 2:
                first_vertex = ([min_x,min_y])
                second_vertex = ([max_x,max_y])
            else:
                first_vertex = ([min_x,max_y])
                second_vertex = ([max_x,min_y])
            line = shapely.LineString([first_vertex, second_vertex])
            for cutted_geom in shapely.ops.split(geom, line).geoms:
                if cutted_geom.area == 0:
                    continue
                if len(list(cutted_geom.interiors)):
                    geoms_with_interiors.append(cutted_geom)
                else:
                    geoms_without_interiors.append(cutted_geom)
        cycles += 1
    return geoms_without_interiors


def generate_object_list(gdf,default_height,height_scale):
    #Create a List of 3D Geometries (objects) out of the OSM Data
    #Geometries need to be preprocessed (Correct Type, No Holes, vertices in clockwise order, scaling, ...)
    
    object_list = []

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
            if hasattr(row,"highway") and not pandas.isna(row.highway):
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
        height = get_object_height(row, default_height)
        if height > 0:
            object_height = height * 1000 * height_scale
            if object_height > 0.5:
                object.append(object_height) #height in meter * 1000 = height in millimeter; height in millimeter gets then scaled down
                object.append(1) #height offset = 1 since it is on top of the base
            else:
                object.append(object_height + 1) #height in meter * 1000 = height in millimeter; height in millimeter gets then scaled down
                object.append(0) #height offset = 0 since it needs to be embedded in base, the embedded millimeter is added to the object height
        else:
            object.append(1) #It is 1 mm deep embedded in base
            object.append(0) #height offset = 0 since it needs to be embedded in base, the embedded millimeter is added to the object height

        if object[0].area > 0:
            object_list.append(object)
    return(object_list)
        
def preprocess_objects(object_list,bbox,target_size,base_scaling_factor):
    parameters = []
    base_size = target_size*base_scaling_factor
    id = 0
    for object in object_list:
        #Create a List with all parameters for multiprocessing
        parameters.append([object,bbox,target_size,base_size,id])
        id += 1
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
    object,bbox,target_size,base_size,id = args
    print(f"preprocessing id: {id}")
    processed_objects = []
    #If Polygon has holes, remove them by splitting it into multiple Polygons
    if len(list(object[0].interiors)):
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
        #in case the geometry is invalid we try to fix it
        if not shapely.is_valid(geometry):
            valid_geom = shapely.make_valid(geometry)
            if isinstance(valid_geom, shapely.geometry.MultiPolygon):
                for geom in valid_geom.geoms:
                    if isinstance(geom,shapely.geometry.Polygon):
                        processed_objects.append([geom,object[1],object[2]])
            else:
                if isinstance(valid_geom,shapely.geometry.Polygon):
                    processed_objects.append([valid_geom,object[1],object[2]])
        else:
            processed_objects.append([geometry,object[1],object[2]])
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

def prepare_3d_mesh(preprocessed_objects, target_size, base_scaling_factor, base_thickness, base_generation=True, object_generation=False):
    #Generation of 3D Mesh out of 2D Shapes with height
    vertices = []
    faces = []
    base_size = target_size * base_scaling_factor

    if base_generation == True:
        # Generate solid base
        base_vertices, base_faces = create_solid_base(base_size, base_thickness)
        # only half array will be used since we only need top OR bottom for the 2D object.
        base_geometry = create_geometry(base_vertices,range(int(len(base_vertices)/2)))
        vertices.extend(base_vertices)
        faces.extend(base_faces)
    else:
        # Generate solid base
        base_vertices, base_faces = create_solid_base(target_size, base_thickness, ((base_scaling_factor - 1) / 2) * target_size)
        # only half array will be used since we only need top OR bottom for the 2D object.
        base_geometry = create_geometry(base_vertices,range(int(len(base_vertices)/2)))

    
    if object_generation == True:
        id = 0
        for object in preprocessed_objects:
            print(f"processing object {id} of {len(preprocessed_objects)}")
            exterior_coords = list(object[0].exterior.coords)
            height = object[1]
            height_offset = object[2]
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
                        vertices.extend(create_vertices_list(exterior_coords, base_thickness, height, height_offset))
                        faces = create_add_faces(base_index, exterior_coords, vertices, faces)
                else:
                    exterior_coords = list(intersection.exterior.coords)
                    #Get Base Index (Before adding new vertices of this object)
                    base_index = len(vertices)
                    vertices.extend(create_vertices_list(exterior_coords, base_thickness, height, height_offset))
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
            base_low =  base_object[2]
            base_high = base_object[1] + base_object[2]
            cutting_low = cutting_object[2]
            cutting_high = cutting_object[1] + cutting_object[2]
            if not (
                (base_low >= cutting_low and base_low <= cutting_high) or
                (base_high >= cutting_low and base_high <= cutting_high) or
                (cutting_low >= base_low and cutting_low <= base_high) or
                (cutting_high >= base_low and cutting_high <= base_high)
            ):
                continue
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
                new_object_list.append(([polygon,base_object[1],base_object[2]]))
        if isinstance(base_object[0], shapely.geometry.Polygon):
            new_object_list.append(base_object)
    return new_object_list



def cut_all_categories(object_list_buildings,object_list_paths,object_list_water,object_list_greens,object_list_base):
    #If water is in a green area (like a river or a fointan) or a green is fully enclosed by water (an island) we need to cut the other parts. Otherwise fountains get lost in the end-result.
    #Since we need to implement this cutting algorythm anyways we also cut buildings and paths out of the other objects to prevent overlapping and have a nicer print result.
    #Order of objects from top to bottom:
    # buildings (They are high, so no cutting required)
    # paths
    # water (except the green is fully enclosed by water)
    # green
    # base

    #remove stuff from water areas
    object_list_water = cut_two_categories(object_list_water,object_list_paths)
    #This is a special case due to the island detection
    object_list_water = cut_two_categories(object_list_water,object_list_greens,detect_islands=True)

    #remove stuff from the green areas
    object_list_greens = cut_two_categories(object_list_greens,object_list_paths)
    object_list_greens = cut_two_categories(object_list_greens,object_list_water)

    #embed stuff into base plate
    #object_list_base = cut_two_categories(object_list_base,object_list_buildings)
    object_list_base = cut_two_categories(object_list_base,object_list_paths)
    object_list_base = cut_two_categories(object_list_base,object_list_water)
    object_list_base = cut_two_categories(object_list_base,object_list_greens)

    return(object_list_buildings,object_list_paths,object_list_water,object_list_greens,object_list_base)