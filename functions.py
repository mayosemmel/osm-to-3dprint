"""This module provides processing functions for osm-to-3dprint."""

import multiprocessing
import numpy as np
import osmnx as ox
import shapely
import pandas
import geopandas
from stl import mesh
from triangulation import ear_clipping_triangulate

def truncate_float(float_number, decimal_places):
    """Function truncating float number to given decimal places."""
    multiplier = 10 ** decimal_places
    return int(float_number * multiplier) / multiplier

def square_bbox_from_center_point(lat, lon, distance):
    """Genereate a square on a map with a given center coordiante."""
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
    """Make sure the vertices are resulting in a square. If not the quare will be slightly adjusted."""
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
    """Fetch building footprints within the bounding box."""
    try:
        if location_type == "buildings":
            gdf = ox.features_from_bbox( bbox, tags = {
                'building': True, 'historic': ['citywalls']})
        elif location_type == "paths":
            gdf = ox.features_from_bbox( bbox, tags = {
                'highway': True, 'man_made': ['pier'], 'railway': True})
        elif location_type == "water":
            gdf = ox.features_from_bbox( bbox, tags = {
                'natural': ['water','reef'],'landuse': ['basin','salt_pond'], 'leisure': ['swimming_pool']})
        elif location_type == "green":
            gdf = ox.features_from_bbox( bbox,tags ={
                'landuse': ['forest','meadow','grass','allotments','flowerbed','orchard','plant_nursery','vineyard','cemetery','recreation_ground',
                'village_green'],'leisure': ['garden','park','pitch'],'natural': ['grassland','scrub','wood']})
        else:
            #invalid location type -> returning empty polygon
            gdf = ([[shapely.Polygon(),0]])
    except ox._errors.InsufficientResponseError: # pylint: disable=protected-access
        # Return Empty Polygon if location type is not found
        gdf = ([[shapely.Polygon(),0]])
    return gdf

def get_special_area_height(area_type, sub_area_type):
    """This functions defines the height of special areas."""
    height = 0

    if area_type == 'highway':
        #Stuff that needs to be the same level as the base plate
        height =  -1
    elif area_type == 'railway':
        height =  0.5

    if sub_area_type in ('water', 'reef', 'grassland', 'scrub', 'wood', 'garden', 'park', 'swimming_pool', 'pitch', 'basin', 'salt_pond', 'grass',
                         'allotments', 'flowerbed', 'village_green'):
        #Stuff that needs to be the same level as the base plate
        height = -1
    elif sub_area_type in ('meadow', 'cemetery'):
        height =  0.5
    elif sub_area_type == 'scrub':
        #scrubs are between 1 and 2 meters high
        height = 1
    elif sub_area_type in ('pier', 'plant_nursery', 'vineyard', 'recreation_ground'):
        height =  2
    elif sub_area_type == 'orchard':
        height =  4
    elif sub_area_type in ('wood', 'forest'):
        #woods probably about 15 meters high
        height = 15

    return height

def get_object_height(row, default_height=10):
    """
    Query or Define object heights
    Areas like water and grass need to be on same level as base plate
    Areas like forest or scrub need to be elevated slightly
    Paths are also even with base plate
    #TODO: What about bridges?!
    """

    special_area_attributes = {'natural', 'landuse', 'leisure', 'highway', 'man_made', 'railway'}
    attribute = special_area_attributes.intersection(row.index)
    if attribute:
        attribute = attribute.pop()
        height = get_special_area_height(attribute, row[attribute])
        if height != 0:
            return height

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
            if isinstance(height, str):
                try:
                    height_value = float(height.replace('m', '').strip())
                    if attr == 'building:levels':
                        height = height_value * 3  # Assuming 3 meters per level
                    elif attr == 'historic:citywalls':
                        height = default_citywall_height
                    else:
                        height = height_value
                    return height
                except ValueError:
                    continue
    return default_height

def create_solid_base(base_size, base_thickness=2, offset=0):
    """Generate the Base-Square with given size and thickness"""
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

def create_planar_face(vertices,z=0):
    """Convert vertices to shapely polygon and export as triangles."""
    vertices_converted = []
    for vertex in vertices:
        vertices_converted.append(vertex[z])

    polygon = shapely.Polygon(vertices_converted)
    return ear_clipping_triangulate(polygon)

def create_geometry(vertices,indicies):
    """Convert custom vertices/indicies format to a shapely polygon."""
    geometry_coords = []
    for index in indicies:
        x = vertices[index][0]
        y = vertices[index][1]
        geometry_coords.append((x,y))
    return shapely.Polygon(geometry_coords)

def create_vertices_list(exterior_coords, base_thickness, height, height_offset=1):
    """Define upper and lower layer and convert to Vertices List"""
    #Vertices are organized as follows
    #bottom_coords=[x,y,z]
    #top_coords=[x,y,z]
    #coord_pair=[bottom_coords,top_coords]
    #vertices_list=[coord_pair,coord_pair,...]
    vertices_list = []
    precision = 4
    for coord in exterior_coords:
        v_bottom = (round(coord[0], precision), round(coord[1], precision), round(base_thickness + height_offset, precision))
        v_top = (round(coord[0], precision), round(coord[1], precision), round(height + base_thickness + height_offset, precision))
        vertices_list.append([v_bottom, v_top])
    return vertices_list

def create_side_faces(vertices):
    """Generate side faces as zig-zag triangles"""
    side_faces = []
    for i in range(len(vertices) - 1):
        bottom1 = vertices[i][0]
        bottom2 = vertices[i+1][0]
        top1 = vertices[i][1]
        top2 = vertices[i+1][1]

        side_faces.append([bottom1, bottom2, top1])
        side_faces.append([top1, bottom2, top2])
    return side_faces

def scale_polygon(exterior_coords, bbox, target_size, base_size):
    """Scale Polygon to model size."""
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
    center_offset_x = (base_size - target_size) / 2
    center_offset_y = (base_size - target_size) / 2

    for i, coord in enumerate(exterior_coords):
        exterior_coords[i] = list(coord)
        exterior_coords[i][0] = round(((exterior_coords[i][0] - south_lng) * scale_x) + center_offset_x, 6)
        exterior_coords[i][1] = round(((exterior_coords[i][1] - south_lat) * scale_y) + center_offset_y, 6)
    return shapely.Polygon(exterior_coords)

def cut_polygon(geometry):
    """Polygon get split as many times as required to cut all holes open. Returns list of polygons without holes."""
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
                if len(list(cutted_geom.interiors)) > 0:
                    geoms_with_interiors.append(cutted_geom)
                else:
                    geoms_without_interiors.append(cutted_geom)
        cycles += 1
    return geoms_without_interiors

def get_path_width(path_subtype):
    """This functions defines the width of path. Values are manually set."""
    if path_subtype == "pier":
        width = 0.000004
    elif path_subtype in ("living_street", "service", "pedestrian", "track", "footway", "bridleway", "path", "sidewalk", "crossing", "traffic_island",
                          "cycleway", "residential"):
        width = 0.000025
    elif path_subtype in ("tertiary_link", "secondary_link"):
        width = 0.00004
    elif path_subtype == "tertiary":
        width = 0.000045
    elif path_subtype in ("primary_link", "trunk_link"):
        width = 0.00005
    elif path_subtype == "secondary":
        width = 0.000055
    elif path_subtype in ("primary", "motorway_link"):
        width = 0.00006
    elif path_subtype == "trunk":
        width = 0.000075
    elif path_subtype == "motorway":
        width = 0.0001
    else:
        print(f"Unclassified path width for type {path_subtype}! Using Default Value!")
        width = 0.000025
    return width

def generate_object_list(gdf,default_height,height_scale):
    """
    Create a List of 3D Geometries (objects) out of the OSM Data
    Geometries need to be preprocessed (Correct Type, No Holes, vertices in clockwise order, scaling, ...)
    """

    object_list = []

    for idx, row in gdf.iterrows():
        #We will have a geometry in the the first place [0] and the height in the second place [1]
        geo_object = []
        #Check if we can process the object. Points and other stuff are not implemented (yet).
        if not isinstance(row['geometry'], (shapely.geometry.Polygon, shapely.geometry.LineString)):
            print(f"Object of Type {idx[0]} with id {idx[1]} is not implemented (yet).")
            continue

        #Get the geometry out of the raw data
        geo_object.append(row['geometry'])
        # If Object is a string convert to polygon
        if isinstance(geo_object[0], shapely.geometry.LineString):
            #manual definition of path width depending on type
            if hasattr(row,"highway") and not pandas.isna(row.highway):
                geo_object[0] = shapely.buffer(geo_object[0], get_path_width(row.highway))
            elif hasattr(row,"man_made") and not pandas.isna(row.man_made):
                geo_object[0] = shapely.buffer(geo_object[0], get_path_width(row.man_made))
            elif hasattr(row,"railway") and not pandas.isna(row.railway) and not (row.tunnel == "yes" or row.railway == 'razed'):
                geo_object[0] = shapely.buffer(geo_object[0], 0.00002)
            else:
                geo_object[0] = shapely.buffer(geo_object[0], 0.000025)

        #Get Object height
        height = get_object_height(row, default_height)
        if height > 0:
            object_height = height * 1000 * height_scale
            if object_height > 0.5:
                geo_object.append(object_height) #height in meter * 1000 = height in millimeter; height in millimeter gets then scaled down
                geo_object.append(1) #height offset = 1 since it is on top of the base
            else:
                geo_object.append(object_height + 1) #height in meter * 1000 = height in millimeter; height in millimeter gets then scaled down
                geo_object.append(0) #height offset = 0 since it needs to be embedded in base, the embedded millimeter is added to the object height
        else:
            geo_object.append(1) #It is 1 mm deep embedded in base
            geo_object.append(0) #height offset = 0 since it needs to be embedded in base

        if geo_object[0].area > 0:
            object_list.append(geo_object)
    return object_list

def preprocess_objects_meta(object_list,bbox,target_size,base_scaling_factor, scale=True):
    """Meta function for preprocessing. Primarily needed for multithreading."""

    base_size = target_size*base_scaling_factor

    #Call Preprocessing Function in Multiprocessing
    print("starting preprocessing with", multiprocessing.cpu_count(), "CPU Cores")
    #Be careful if reenabling this! It seems to trigger some kind of race condition where SOMETIMES a logic error in the end-result happens.
    #with multiprocessing.Pool(multiprocessing.cpu_count()) as p:
    #    meta_object_list = p.map(preprocess_objects,object_list)
    #    p.close()
    #    p.join()
    ###################################
    #This is only for debugging without multiprocessing
    meta_object_list = []
    for single_object in object_list:
        meta_object_list.append(preprocess_objects(single_object))
    ######################################
    preprocessed_objects = []
    for meta_object in meta_object_list:
        for geo_object in meta_object:
            if scale:
                geo_object[0] = scale_polygon(geo_object[0].exterior.coords,bbox,target_size,base_size)
            preprocessed_objects.append(geo_object)
    print("preprocessing done")
    return preprocessed_objects

def order_points_clockwise(geometry_list):
    """check if points of polygon are clockwise ordered"""
    clockwise_geometries = []
    for geometry in geometry_list:
        if shapely.algorithms.cga.signed_area(geometry.exterior) > 0:
            geometry = shapely.Polygon(reversed(geometry.exterior.coords))
        clockwise_geometries.append(geometry)
    return clockwise_geometries

def preprocess_objects(geo_object):
    """Preprocess geometries to have the correct type, no holes, vertices in clockwise order, scaling, ...)"""
    print('.', end='')
    processed_objects = []
    #If Polygon has holes, remove them by splitting it into multiple Polygons
    if len(list(geo_object[0].interiors)) > 0:
        geometry_list = cut_polygon(geo_object[0])
    else:
        #nothing to cut
        geometry_list = [geo_object[0]]

    clockwise_list = order_points_clockwise(geometry_list)

    for geometry in clockwise_list:
        #in case the geometry is invalid we try to fix it
        if not shapely.is_valid(geometry):
            valid_geom = shapely.make_valid(geometry)
            if isinstance(valid_geom, shapely.geometry.MultiPolygon):
                for geom in valid_geom.geoms:
                    if isinstance(geom,shapely.geometry.Polygon):
                        processed_objects.append([geom,geo_object[1],geo_object[2]])
            else:
                if isinstance(valid_geom,shapely.geometry.Polygon):
                    processed_objects.append([valid_geom,geo_object[1],geo_object[2]])
        else:
            processed_objects.append([geometry,geo_object[1],geo_object[2]])
    return processed_objects

def create_faces(vertices):
    """create faces for all sides of the 3d-part"""
    faces = []

    # Create side faces
    faces.extend(create_side_faces(vertices))

    # Create top face
    faces.extend(create_planar_face(vertices,z=0))

    # Create bottom face
    faces.extend(create_planar_face(vertices,z=1))

    return faces

def make_geometries_valid(invalid_objects):
    """in case the geometry is invalid we try to fix it"""
    valid_objects = []
    for geo_object in invalid_objects:
        if not shapely.is_valid(geo_object[0]):
            valid_geom = shapely.make_valid(geo_object[0])
            if isinstance(valid_geom, shapely.geometry.MultiPolygon):
                for geom in valid_geom.geoms:
                    if isinstance(geom,shapely.geometry.Polygon):
                        valid_objects.append([geom,geo_object[1],geo_object[2]])
            elif isinstance(valid_geom,shapely.geometry.Polygon):
                valid_objects.append([valid_geom,geo_object[1],geo_object[2]])
        else:
            valid_objects.append(geo_object)
    return valid_objects

def create_base_geometry(target_size, base_scaling_factor):
    """Generate the base geometry as shapely polygon"""
    base_size = target_size * base_scaling_factor

    base_geometry = shapely.Polygon((
        (0, 0),
        (0, base_size),
        (base_size, base_size),
        (base_size, 0)))
    return base_geometry

def create_inner_base_geometry(target_size, base_scaling_factor):
    """Generate the inner base geometry as shapely polygon"""
    offset_inner_base = (target_size * (base_scaling_factor - 1)) / 2
    inner_base_geometry = shapely.Polygon((
        (offset_inner_base, offset_inner_base),
        (offset_inner_base, target_size + offset_inner_base),
        (target_size + offset_inner_base, target_size + offset_inner_base),
        (target_size + offset_inner_base, offset_inner_base)))
    return inner_base_geometry

def generate_base(target_size, base_scaling_factor, base_thickness):
    """Generate solid base"""
    faces = []

    base_geometry = create_base_geometry(target_size, base_scaling_factor)
    inner_base_geometry = create_inner_base_geometry(target_size, base_scaling_factor)

    vertices = create_vertices_list(base_geometry.exterior.coords, 0, base_thickness, height_offset=0)
    faces.extend(create_faces(vertices))

    # Frame around Base.
    base_frame = shapely.difference(base_geometry, inner_base_geometry)
    # This line is to cut the base diagonal to have a polygon without holes.
    line = shapely.LineString([[0,0], [target_size * base_scaling_factor,target_size * base_scaling_factor]])
    base_frame_objects = []
    for base_frame_part in shapely.ops.split(base_frame, line).geoms:
        base_frame_objects.append([base_frame_part,1,0])
    for geo_object in base_frame_objects:
        exterior_coords = list(geo_object[0].exterior.coords)
        vertices = create_vertices_list(exterior_coords, base_thickness, geo_object[1], geo_object[2])
        faces.extend(create_faces(vertices))

    return faces

def remove_overhangs_over_base(target_size, base_scaling_factor, base_thickness, geo_object):
    """Remove any overhangs over the base plate. This is required since some objects start within the bbox but end outside of it.
    If this is the case we have object which are a lot to big and are not printable.
    After Cleanup the Faces of the final object are added to the whole list"""
    faces = []
    inner_base_geometry = create_inner_base_geometry(target_size, base_scaling_factor)
    height = geo_object[1]
    height_offset = geo_object[2]
    if shapely.intersects(inner_base_geometry, geo_object[0]):
        intersection = shapely.intersection(inner_base_geometry, geo_object[0])
        if isinstance(intersection, shapely.geometry.MultiPolygon):
            for polygon in list(intersection.geoms):
                exterior_coords = list(polygon.exterior.coords)
                vertices = create_vertices_list(exterior_coords, base_thickness, height, height_offset)
                faces.extend(create_faces(vertices))
        elif isinstance(intersection, shapely.geometry.Polygon):
            exterior_coords = list(intersection.exterior.coords)
            vertices = create_vertices_list(exterior_coords, base_thickness, height, height_offset)
            faces.extend(create_faces(vertices))
    return faces

def prepare_3d_mesh(preprocessed_objects, target_size, base_scaling_factor, base_thickness, base_generation=True, object_generation=False):
    """Generation of 3D Mesh out of 2D Shapes with height"""
    faces = []

    if base_generation:
        # Generate solid base
        faces.extend(generate_base(target_size, base_scaling_factor, base_thickness))

    if object_generation:
        #in case the geometry is invalid we try to fix it
        preprocessed_valid_objects = make_geometries_valid(preprocessed_objects)

        for geo_object in preprocessed_valid_objects:
            faces.extend(remove_overhangs_over_base(target_size, base_scaling_factor, base_thickness, geo_object))

    faces = np.array(faces)

    return faces

def save_to_stl(faces, filename):
    """Save generated mesh to a .stl file"""
    mesh_data = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
    for i, face in enumerate(faces):
        for j in range(3):
            mesh_data.vectors[i][j] = face[j]

    # Create a new 3D plot
    #figure = pyplot.figure()
    #axes = figure.add_subplot(projection='3d')
    #axes.add_collection3d(mplot3d.art3d.Poly3DCollection(mesh_data.vectors))
    #scale = mesh_data.points.flatten()
    #axes.auto_scale_xyz(scale, scale, scale)
    #pyplot.show()

    mesh_data.save(filename)

def cut_two_categories(base_category,cutting_category,detect_islands=False):
    """Cut the object of two categories with each other. This prevents overlapping objects."""
    new_object_list = []
    for base_object in base_category:
        for cutting_object in cutting_category:
            base_low =  base_object[2]
            base_high = base_object[1] + base_object[2]
            cutting_low = cutting_object[2]
            cutting_high = cutting_object[1] + cutting_object[2]
            if not (
                (base_low > cutting_low and base_low < cutting_high) or
                (base_high > cutting_low and base_high < cutting_high) or
                (cutting_low > base_low and cutting_low < base_high) or
                (cutting_high > base_low and cutting_high < base_high) or
                (cutting_low == base_low and cutting_high == base_high)
            ):
                continue
            if shapely.intersects(base_object[0],cutting_object[0]):
                # Island Detection Example:
                # When we remove the water area from the green area and no area is left,
                # the green area IS fully surroundend and therefore is an island.
                if detect_islands and not shapely.difference(cutting_object[0],base_object[0]).area == 0:
                    continue
                base_object[0] = shapely.difference(base_object[0],cutting_object[0])
        print('.', end='')
        if base_object[0].area == 0:
            continue
        #after all cuts we append the object to the list
        if isinstance(base_object[0], shapely.geometry.MultiPolygon):
            for polygon in base_object[0].geoms:
                new_object_list.append(([polygon,base_object[1],base_object[2]]))
        if isinstance(base_object[0], shapely.geometry.Polygon):
            new_object_list.append(base_object)
    print('.')
    return new_object_list



def cut_all_categories(object_list_buildings,object_list_paths,object_list_water,object_list_greens,object_list_base):
    """Cut all given categories in the correct order."""
    # If water is in a green area (like a river or a fointan) or a green is fully enclosed by water (an island) we need to cut the other parts.
    # Otherwise fountains get lost in the end-result.
    # Since we need to implement this cutting algorythm anyways we also cut buildings and paths out of the other objects to prevent overlapping and
    # have a nicer print result.

    # Order of objects from top to bottom:
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

    #sometimes buildings and paths are overlapping
    object_list_paths = cut_two_categories(object_list_paths,object_list_buildings)


    #embed stuff into base plate
    #object_list_base = cut_two_categories(object_list_base,object_list_buildings)
    object_list_base = cut_two_categories(object_list_base,object_list_paths)
    object_list_base = cut_two_categories(object_list_base,object_list_water)
    object_list_base = cut_two_categories(object_list_base,object_list_greens)

    return(object_list_buildings,object_list_paths,object_list_water,object_list_greens,object_list_base)
