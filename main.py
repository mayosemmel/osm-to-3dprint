"""
Starting Point of OpenStreetMap to 3D Print file Generation. All Parameters are defined in this module.
"""

import math
import shapely
import geopy.distance
import functions as f


### TODO ###
#
# multithreading/multitasking for triangulation
# refactor organization of vertices
# triangulation in 3 dimensional objects
# Ongoing: Code Cleanup according to PyLint



def main():
    """
    "main" function, will be called as first function in the program.
    """
    #size of print in mm excluding "frame" overhang
    target_size = 180
    #Thickness of Base Plate in mm, needs to be bigger than 1 mm
    base_thickness = 2
    #everything bigger than one will create a "frame" around the actual landscape
    base_scaling_factor = 1.1
    #How much higher than correct scaling should the buildings be?
    #With bigger Parts of Land this factor should be higher
    height_scaling_factor = 1.2
    default_building_height = 9 # in meter

    #Get BBox from vertices
    #bbox = f.square_bbox_from_vertices( min Latitude, min Longitude, max Latitude, max Longitude )
    #bbox = f.square_bbox_from_vertices(52.35893, 4.87123, 52.38351, 4.93389)  #Amsterdam
    #bbox = f.square_bbox_from_vertices(49.27478, 10.85891, 49.27973, 10.86771) #Suddersdorf
    bbox = f.square_bbox_from_vertices(51.757883, -1.266515, 51.759302, -1.263503) #Oxford University (Polygon with Holes)
    #bbox = f.square_bbox_from_vertices(49.44759, 11.06375, 49.45976, 11.09048) #Nürnberg Zentrum
    #bbox = f.square_bbox_from_vertices(49.40804, 11.07375, 49.42298, 11.11181) #Nürnberg Rangierbahnhof
    #bbox = f.square_bbox_from_vertices(49.43102, 11.09297, 49.41910, 11.11254) #Debugging
    #bbox = f.square_bbox_from_vertices(49.38656, 11.03946, 49.41215, 11.07800) # Nürnberg Hafen
    #bbox = f.square_bbox_from_vertices(49.56984, 10.58769, 49.58768, 10.63133) # Neustadt Aisch

    #Get BBox from center point and square size
    #bbox = f.square_bbox_from_center_point(Latitude, Longitude, Square Size in Meter)
    #bbox = f.square_bbox_from_center_point(48.32950556656733, 10.90461275575229, 1000) #Augsburg
    #bbox = f.square_bbox_from_center_point(49.453675, 11.077115, 1000) #Nürnberg Zentrum
    #bbox = f.square_bbox_from_center_point(48.76336, 11.42484, 2000) #Ingolstadt

    #Define what should be generated
    base_plate = True
    buildings = True
    paths = True
    water = True
    green = True

    #Get the Scaling Factor of the whole model, use this block if the scaling should be correct (same factor for lat/long and height).
    long1, lat1, long2, lat2 = bbox
    diagonale = geopy.distance.geodesic((lat1, long1), (lat2, long2)).km
    side_length = diagonale / math.sqrt(2)
    model_scaling_factor = (target_size / 1000) / (side_length * 1000) # mm / 1000 = meter ## km * 1000 = meter
    height_scale = model_scaling_factor * height_scaling_factor

    #Preparation of 2D Object with height as metadata for later processing
    #Ojects are organized as follows
    #object[0] -> shapely Polygon
    #object[1] -> height in mm
    #object[2] -> height offset to the base in mm

    #Generation of Base Plate
    #We have a Part of the Base Plate which is solid, this will be generated automatically
    #The first mm on the upper side needs to be cut by paths and stuff like that so it is on the same level as the base plate

    if base_plate:
        if base_thickness <= 1:
            raise ValueError("Invalid Base Thickness. Needs to be bigger than 1 mm.")
        base = shapely.Polygon((
            (0, 0),
            (target_size * base_scaling_factor, 0),
            (target_size * base_scaling_factor, target_size * base_scaling_factor),
            (0, target_size * base_scaling_factor)))
        object_list_base = ([[base,1,0]])
        object_list_base = f.preprocess_objects_meta(object_list_base,bbox,target_size,base_scaling_factor,scale=False)
        base_thickness = base_thickness - 1
    else:
        #for further processing we need some "empty" Polygon as dummy data
        object_list_base = ([[shapely.Polygon(),0,0]])

    #Generation of Buildings
    if buildings:
        gdf = f.fetch_location_data(bbox, "buildings")
        object_list_buildings = f.generate_object_list(gdf,default_building_height,height_scale)
        object_list_buildings = f.preprocess_objects_meta(object_list_buildings,bbox,target_size,base_scaling_factor,scale=True)
    else:
        #for further processing we need some "empty" Polygon as dummy data
        object_list_buildings = ([[shapely.Polygon(),0,0]])

    #Generation of Paths
    if paths:
        gdf = f.fetch_location_data(bbox, "paths")
        object_list_paths = f.generate_object_list(gdf,default_building_height,height_scale)
        object_list_paths = f.preprocess_objects_meta(object_list_paths,bbox,target_size,base_scaling_factor,scale=True)
    else:
        #for further processing we need some "empty" Polygon as dummy data
        object_list_paths = ([[shapely.Polygon(),0,0]])

    #Generation of Water
    if water:
        gdf = f.fetch_location_data(bbox, "water")
        object_list_water = f.generate_object_list(gdf,default_building_height,height_scale)
        object_list_water = f.preprocess_objects_meta(object_list_water,bbox,target_size,base_scaling_factor,scale=True)
    else:
        #for further processing we need some "empty" Polygon as dummy data
        object_list_water = ([[shapely.Polygon(),0,0]])

    #Generation of "Green Areas" like Forest and Meadow
    if green:
        gdf = f.fetch_location_data(bbox, "green")
        object_list_greens = f.generate_object_list(gdf,default_building_height,height_scale)
        object_list_greens = f.preprocess_objects_meta(object_list_greens,bbox,target_size,base_scaling_factor,scale=True)
    else:
        #for further processing we need some "empty" Polygon as dummy data
        object_list_greens = ([[shapely.Polygon(),0,0]])

    #If water is in a green area (like a river or a fointan) or a green is fully enclosed by water (an island) we need to cut the other parts.
    #Otherwise fountains get lost in the end-result.
    #Since we need to implement this cutting algorythm anyways we also cut buildings and paths out of the other objects to prevent overlapping and
    #have a nicer print result.
    print("starting to cut the layer categories with each other")
    all_categories = f.cut_all_categories(object_list_buildings,object_list_paths,object_list_water,object_list_greens,object_list_base)
    object_list_buildings,object_list_paths,object_list_water,object_list_greens,object_list_base = all_categories


    #We need Preprocessing again because we did a cut to severeal objects. Therefore there could be interiors and other problems again.
    #Generation of Base Plate
    if base_plate:
        preprocessed_base = f.preprocess_objects_meta(object_list_base,bbox,target_size,base_scaling_factor, scale=False)
        faces = f.prepare_3d_mesh(
            preprocessed_base, target_size, base_scaling_factor, base_thickness, base_generation=True, object_generation=True)
        f.save_to_stl(faces, 'export/base.stl')
        print("generation of base plate completed")

    #Generation of Buildings
    if buildings and len(object_list_buildings) > 0:
        preprocessed_buildings = f.preprocess_objects_meta(object_list_buildings,bbox,target_size,base_scaling_factor,scale=False)
        faces = f.prepare_3d_mesh(
            preprocessed_buildings, target_size, base_scaling_factor, base_thickness, base_generation=False, object_generation=True)
        f.save_to_stl(faces, 'export/buildings.stl')
        print("generation of buildings completed")

    #Generation of Paths
    if paths and len(object_list_paths) > 0:
        preprocessed_paths = f.preprocess_objects_meta(object_list_paths,bbox,target_size,base_scaling_factor,scale=False)
        faces = f.prepare_3d_mesh(
            preprocessed_paths, target_size, base_scaling_factor, base_thickness, base_generation=False, object_generation=True)
        f.save_to_stl(faces, 'export/paths.stl')
        print("generation of paths completed")

    #Generation of Water
    if water and len(object_list_water) > 0:
        preprocessed_water = f.preprocess_objects_meta(object_list_water,bbox,target_size,base_scaling_factor,scale=False)
        faces = f.prepare_3d_mesh(
            preprocessed_water, target_size, base_scaling_factor, base_thickness, base_generation=False, object_generation=True)
        f.save_to_stl(faces, 'export/water.stl')
        print("generation of water completed")

    #Generation of "Green Areas" like Forest and Meadow
    if green and len(object_list_greens) > 0:
        preprocessed_greens = f.preprocess_objects_meta(object_list_greens,bbox,target_size,base_scaling_factor,scale=False)
        faces = f.prepare_3d_mesh(
            preprocessed_greens, target_size, base_scaling_factor, base_thickness, base_generation=False, object_generation=True)
        f.save_to_stl(faces, 'export/greens.stl')
        print("generation of greens completed")


if __name__ == "__main__":
    main()
