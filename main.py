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
# triangulation errors
# Priority of Layers -> Water cuts out in Green, Paths cut out Water and Green
# Performance!! Multithreading is working, but still it is quite slow
# maybe implement triangulation by myself?


def main():
    target_size = 180
    base_thickness = 2
    #This highly depends on the astetics of the city/village you are trying to print. In Big Cities a value ~25 is mostly good. For villages values ~10 are good.
    max_building_height_mm = 10
    default_building_height= 9
    default_citywall_height= 16
    #bbox = (4.87123, 52.35893, 4.93389, 52.38351)  #Amsterdam
    #bbox = (10.85891, 49.27478, 10.86771, 49.27973) #Suddersdorf
    #bbox = (-1.266515, 51.757883, -1.263503, 51.759302) #Oxford University (Polygon with Holes)
    bbox = (11.06375, 49.44759, 11.09048, 49.45976) #Nürnberg Zentrum
    #bbox = min Longitude , min Latitude , max Longitude , max Latitude 

    #Define what should be generated
    base_plate = False
    buildings = True
    paths = False
    water = False
    green = False

    #Generation of Base Plate
    if base_plate:
        vertices, faces = prepare_mesh(False, bbox, target_size=target_size, max_height_mm=max_building_height_mm, default_height=default_building_height, base_thickness=base_thickness, base_generation=True, object_generation=False)
        save_to_stl(vertices, faces, 'export/standalone_base.stl')
        print(f"generation of base plate completed")

    #Generation of Buildings
    if buildings:
        gdf = fetch_location_data(bbox, "buildings")
        vertices, faces = prepare_mesh(gdf, bbox, target_size=target_size, max_height_mm=max_building_height_mm, default_height=default_building_height, base_thickness=base_thickness, base_generation=False, object_generation=True, default_citywall_height=default_citywall_height)
        save_to_stl(vertices, faces, 'export/buildings_without_base.stl')
        print(f"generation of buildings completed")

    #Generation of Paths
    if paths:
        gdf = fetch_location_data(bbox, "paths")
        vertices, faces = prepare_mesh(gdf, bbox, target_size=target_size, max_height_mm=1, default_height=default_building_height, base_thickness=base_thickness, base_generation=False, object_generation=True)
        save_to_stl(vertices, faces, 'export/paths_without_base.stl')
        print(f"generation of paths completed")

    #Generation of Water
    if water:
        gdf = fetch_location_data(bbox, "water")
        vertices, faces = prepare_mesh(gdf, bbox, target_size=target_size, max_height_mm=0.5, default_height=default_building_height, base_thickness=base_thickness, base_generation=False, object_generation=True)
        save_to_stl(vertices, faces, 'export/water_without_base.stl')
        print(f"generation of water completed")

    #Generation of "Green Areas" like Forest and Meadow
    if green:
        gdf = fetch_location_data(bbox, "green")
        vertices, faces = prepare_mesh(gdf, bbox, target_size=target_size, max_height_mm=0.75, default_height=default_building_height, base_thickness=base_thickness, base_generation=False, object_generation=True)
        save_to_stl(vertices, faces, 'export/greens_without_base.stl')
        print(f"generation of greens completed")


if __name__ == "__main__":
    main()
