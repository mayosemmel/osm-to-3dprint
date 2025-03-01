import osmnx as ox
import shapely
import numpy as np
import subprocess
import concurrent.futures
import multiprocessing
import os
import math
import shapely.prepared
from functions import *
from stl import mesh
from mpl_toolkits import mplot3d
from matplotlib import pyplot
#C++ Tool from https://github.com/mayosemmel/Geometry-Algorithms/tree/master/Triangulation


### TODO ###
#
# fix cut_polygon. We need to completely remove all interiors without leaving the function.
# maybe implement triangulation by myself?
# What about a bbox which is not a square? -> aspect ratio



def main():
    #size of print in mm excluding "frame" overhang
    target_size = 180
    #Thickness of Base Plate in mm
    base_thickness = 2
    #everything bigger than one will create a "frame" around the actual landscape
    scaling_factor = 1.1
    #This highly depends on the astetics of the city/village you are trying to print. In Big Cities a value ~25 is mostly good. For villages values ~10 are good.
    max_building_height_mm = 10
    default_building_height = 9
    path_height = 0.6
    water_height = 0.4
    green_height = 0.6

    #bbox = (4.87123, 52.35893, 4.93389, 52.38351)  #Amsterdam
    #bbox = (10.85891, 49.27478, 10.86771, 49.27973) #Suddersdorf
    #bbox = (-1.266515, 51.757883, -1.263503, 51.759302) #Oxford University (Polygon with Holes)
    #bbox = (11.06375, 49.44759, 11.09048, 49.45976) #Nürnberg Zentrum
    #bbox = (11.07375, 49.40804, 11.11181, 49.42298) #Nürnberg Rangierbahnhof
    #bbox = generate_bbox_from_coords(49.35939, 10.87368, 49.37005, 10.89294) #Buchschwabach
    #bbox = (11.03946, 49.38656, 11.07800, 49.41215) # Nürnberg Hafen
    #bbox = (10.58769, 49.56984, 10.63133, 49.58768) # Neustadt Aisch
    #bbox = min Longitude , min Latitude , max Longitude , max Latitude 
    
    #Get BBox from center point and square size
    #bbox = square_poly(Latitude, Longitude, Square Size in Meter)
    #bbox = square_poly(48.32950556656733, 10.90461275575229, 1000) #Augsburg
    bbox = square_poly(49.453675, 11.077115, 1000) #Nürnberg Zentrum
    
    #Define what should be generated
    base_plate = True
    buildings = True
    paths = True
    water = True
    green = True
    
    #Preparation of 2D Object with height as metadata for later processing

    #Generation of Buildings
    if buildings:
        gdf = fetch_location_data(bbox, "buildings")
        object_list_buildings = generate_object_list(gdf,default_building_height,max_building_height_mm)
    else:
        #for further processing we need some "empty" Polygon as dummy data
        object_list_buildings = ([[shapely.Polygon(),0]])
    
    #Generation of Paths
    if paths:
        gdf = fetch_location_data(bbox, "paths")
        object_list_paths = generate_object_list(gdf,default_building_height,path_height)
    else:
        #for further processing we need some "empty" Polygon as dummy data
        object_list_paths = ([[shapely.Polygon(),0]])
    
    #Generation of Water
    if water:
        gdf = fetch_location_data(bbox, "water")
        object_list_water = generate_object_list(gdf,default_building_height,water_height)
    else:
        #for further processing we need some "empty" Polygon as dummy data
        object_list_water = ([[shapely.Polygon(),0]])
    
    #Generation of "Green Areas" like Forest and Meadow
    if green:
        gdf = fetch_location_data(bbox, "green")
        object_list_greens = generate_object_list(gdf,default_building_height,green_height)
    else:
        #for further processing we need some "empty" Polygon as dummy data
        object_list_greens = ([[shapely.Polygon(),0]])

    #If water is in a green area (like a river or a fointan) or a green is fully enclosed by water (an island) we need to cut the other parts. Otherwise fountains get lost in the end-result.
    #Since we need to implement this cutting algorythm anyways we also cut buildings and paths out of the other objects to prevent overlapping and have a nicer print result.
    print(f"starting to cut the layer categories with each other")
    object_list_buildings,object_list_paths,object_list_water,object_list_greens = cut_all_categories(object_list_buildings,object_list_paths,object_list_water,object_list_greens)

    #Generation of Base Plate
    if base_plate:
        vertices, faces = prepare_3d_mesh(False,target_size, scaling_factor, base_thickness, base_generation=True, object_generation=False)
        save_to_stl(vertices, faces, 'export/standalone_base.stl')
        print(f"generation of base plate completed")

    #Generation of Buildings
    if buildings and len(object_list_buildings) > 0:
        preprocessed_buildings = preprocess_objects(object_list_buildings,bbox,target_size,scaling_factor)
        vertices, faces = prepare_3d_mesh(preprocessed_buildings, target_size, scaling_factor, base_generation=False, object_generation=True)
        save_to_stl(vertices, faces, 'export/buildings_without_base.stl')
        print(f"generation of buildings completed")

    #Generation of Paths
    if paths and len(object_list_paths) > 0:
        preprocessed_paths = preprocess_objects(object_list_paths,bbox,target_size,scaling_factor)
        vertices, faces = prepare_3d_mesh(preprocessed_paths, target_size, scaling_factor, base_generation=False, object_generation=True)
        save_to_stl(vertices, faces, 'export/paths_without_base.stl')
        print(f"generation of paths completed")

    #Generation of Water
    if water and len(object_list_water) > 0:
        preprocessed_water = preprocess_objects(object_list_water,bbox,target_size,scaling_factor)
        vertices, faces = prepare_3d_mesh(preprocessed_water, target_size, scaling_factor, base_generation=False, object_generation=True)
        save_to_stl(vertices, faces, 'export/water_without_base.stl')
        print(f"generation of water completed")

    #Generation of "Green Areas" like Forest and Meadow
    if green and len(object_list_greens) > 0:
        preprocessed_greens = preprocess_objects(object_list_greens,bbox,target_size,scaling_factor)
        vertices, faces = prepare_3d_mesh(preprocessed_greens, target_size, scaling_factor, base_generation=False, object_generation=True)
        save_to_stl(vertices, faces, 'export/greens_without_base.stl')
        print(f"generation of greens completed")


if __name__ == "__main__":
    main()
