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
# Width of Paths
# Priority of Layers -> Water cuts out in Green, Paths cut out Water and Green
# Performance!! Multithreading is working, but still it is quite slow
# maybe implement triangulation by myself?
# What about a bbox which is not a square? -> aspect ratio
# calculate base size only on one place


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
    default_citywall_height = 16
    path_height = 1
    water_height = 0.5
    green_height = 0.75
    base_size = target_size * scaling_factor

    #bbox = (4.87123, 52.35893, 4.93389, 52.38351)  #Amsterdam
    #bbox = (10.85891, 49.27478, 10.86771, 49.27973) #Suddersdorf
    bbox = (-1.266515, 51.757883, -1.263503, 51.759302) #Oxford University (Polygon with Holes)
    #bbox = (11.06375, 49.44759, 11.09048, 49.45976) #Nürnberg Zentrum
    #bbox = (11.07375, 49.40804, 11.11181, 49.42298) #Rangierbahnhof
    #bbox = (10.58769, 49.56984, 10.63133, 49.58768) # Neustadt Aisch
    #bbox = min Longitude , min Latitude , max Longitude , max Latitude 

    #Define what should be generated
    base_plate = True
    buildings = True
    paths = True
    water = True
    green = True
    
    #Generation of Base Plate
    if base_plate:
        vertices, faces = prepare_mesh(False,target_size=target_size, base_size=target_size*scaling_factor, base_thickness=base_thickness, base_generation=True, object_generation=False, scaling_factor=scaling_factor)
        save_to_stl(vertices, faces, 'export/standalone_base.stl')
        print(f"generation of base plate completed")

    #Generation of Buildings
    if buildings:
        gdf = fetch_location_data(bbox, "buildings")

        # Calculate the maximum object height
        max_building_height = gdf.apply(lambda row: get_building_height(row, default_building_height), axis=1).max()
        height_scale = max_building_height_mm / max_building_height
        
        preprocessed_objects = preprocess_objects_meta(gdf,bbox,target_size,base_size=target_size*scaling_factor,default_height=default_building_height,default_citywall_height=default_citywall_height,height_scale=height_scale)
        vertices, faces = prepare_mesh(preprocessed_objects, target_size=target_size, base_size=base_size, base_generation=False, object_generation=True, scaling_factor=scaling_factor)
        save_to_stl(vertices, faces, 'export/buildings_without_base.stl')
        print(f"generation of buildings completed")

    #Generation of Paths
    if paths:
        height_scale = path_height / default_building_height
        gdf = fetch_location_data(bbox, "paths")
        preprocessed_objects = preprocess_objects_meta(gdf,bbox,target_size,base_size=target_size*scaling_factor,default_height=default_building_height,default_citywall_height=default_citywall_height,height_scale=height_scale)
        vertices, faces = prepare_mesh(preprocessed_objects, target_size=target_size, base_size=base_size, base_generation=False, object_generation=True, scaling_factor=scaling_factor)
        save_to_stl(vertices, faces, 'export/paths_without_base.stl')
        print(f"generation of paths completed")

    #Generation of Water
    if water:
        height_scale = water_height / default_building_height
        gdf = fetch_location_data(bbox, "water")
        preprocessed_objects = preprocess_objects_meta(gdf,bbox,target_size,base_size=target_size*scaling_factor,default_height=default_building_height,default_citywall_height=default_citywall_height,height_scale=height_scale)
        vertices, faces = prepare_mesh(preprocessed_objects, target_size=target_size, base_size=base_size, base_generation=False, object_generation=True, scaling_factor=scaling_factor)
        save_to_stl(vertices, faces, 'export/water_without_base.stl')
        print(f"generation of water completed")

    #Generation of "Green Areas" like Forest and Meadow
    if green:
        height_scale = green_height / default_building_height
        gdf = fetch_location_data(bbox, "green")
        preprocessed_objects = preprocess_objects_meta(gdf,bbox,target_size,base_size=target_size*scaling_factor,default_height=default_building_height,default_citywall_height=default_citywall_height,height_scale=height_scale)
        vertices, faces = prepare_mesh(preprocessed_objects, target_size=target_size, base_size=base_size, base_generation=False, object_generation=True, scaling_factor=scaling_factor)
        save_to_stl(vertices, faces, 'export/greens_without_base.stl')
        print(f"generation of greens completed")


if __name__ == "__main__":
    main()
