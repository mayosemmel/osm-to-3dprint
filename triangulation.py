"""This module provides multiple options for triangulation of complex polygons."""

import copy
import sys
import shapely

def split_list(a_list):
    """split a list in half and return both halfs"""
    half = len(a_list)//2
    return a_list[:half], a_list[half:]

def monotone_subdivision(polygon, level=0, direction=0):
    """Make polygon monotone. This is important for triangulation."""
    #get the bounds of the polygon
    min_x,min_y,max_x,max_y = shapely.bounds(polygon)
    #We will iterate through all points and draw a vertical line at each x or y-position.
    #If this line cuts and more than two polygons are the result we will keep the cut.
    #Otherwise we will not cut at this point and continue with the next point.

    if level > 100:
        print(f"Warning! Recursion level is at {level}!")
    direction = level % 2
    return_polygons = []
    is_monotone = True

    #we are splitting the list in half and combine them again so the cuts are more middled and therefore likely less cuts will be done.
    coords1, coords2 = split_list(polygon.exterior.coords)
    coords2.extend(coords1)

    for coord in coords2:
        if direction == 0:
            current_x = coord[0]
            cutting_line = shapely.LineString([[current_x,min_y],[current_x,max_y]])
        else:
            current_y = coord[1]
            cutting_line = shapely.LineString([[min_x,current_y],[max_x,current_y]])
        cutted_geoms = shapely.ops.split(polygon, cutting_line).geoms
        if len(cutted_geoms) > 2 or (level > 10 and len(cutted_geoms) > 1):
            is_monotone = False
            for geom in cutted_geoms:
                if isinstance(geom, shapely.geometry.Polygon) and geom.area > 0:
                    return_polygons.extend(monotone_subdivision(geom,level=level+1))
            break
    if is_monotone:
        return_polygons.append(polygon)

    return return_polygons


def ear_clipping_triangulate(geometry, level=0):
    """triangulate polygon with ear clipping algorythm"""
    #check if points of polygon are clockwise ordered
    if shapely.algorithms.cga.signed_area(geometry.exterior) > 0:
        geometry = shapely.Polygon(reversed(geometry.exterior.coords))
    vertices = list(geometry.exterior.coords)
    #shapely closes polygons with a point on the same position as the first point. This additional point needs to be removed.
    if vertices[0] == vertices[len(vertices)-1]:
        vertices.pop(len(vertices)-1)
    triangles = []

    if len(vertices) < 3:
        print("Unable to complete triangulation, less than three points given.")
        return triangles
    if geometry.area == 0:
        return triangles

    infinite_loop = False
    initial_vertices = copy.deepcopy(vertices)

    while len(vertices) > 0:
        if infinite_loop:
            if level >= sys.getrecursionlimit()*0.9:
                print(initial_vertices)
            triangles = []
            for geom in monotone_subdivision(geometry, level=level+1):
                triangles.extend(ear_clipping_triangulate(geom, level=level+1))
            break
        infinite_loop = True

        for current, _ in enumerate(vertices):
            if len(vertices) < 3:
                if len(vertices) > 0:
                    vertices.pop(current)
                infinite_loop = False
                break

            #index of previous point
            if current == 0:
                previous = len(vertices) - 1
            else:
                previous = current - 1
            #index of next point
            if current == len(vertices) - 1:
                next = 0
            else:
                next = current + 1
            #create a triangle out of the previous, current and next point
            triangle = shapely.Polygon((vertices[previous],vertices[current], vertices[next]))
            #check if triangle is a line (points coliniar) or a line -> area=0
            #check if triangle points are clockwise ordered
            if triangle.area == 0 or shapely.algorithms.cga.signed_area(triangle.exterior) < 0:
                #check if any other point is in triangle
                point_in_triangle = False
                for vertex in initial_vertices:
                    vertex = shapely.Point(vertex)
                    if triangle.contains(vertex):
                        point_in_triangle = True
                        break

                if not point_in_triangle:
                    #Check if two points are the same point
                    #If we have three points which are clockwise ordered or coliniar we add the triangle to the triangle list
                    if not (vertices[previous] == vertices[current] or vertices[next] == vertices[current] or vertices[next] == vertices[previous]):
                        triangles.append((vertices[previous],vertices[current],vertices[next]))
                    #remove vertice from array
                    vertices.pop(current)
                    infinite_loop = False
                    break
    return triangles
