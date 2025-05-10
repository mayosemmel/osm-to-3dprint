import shapely

def split_list(a_list):
    half = len(a_list)//2
    return a_list[:half], a_list[half:]

def monotoneSubdivision(polygon):
    #get the bounds of the polygon, we just need the y values
    min_x,min_y,max_x,max_y = shapely.bounds(polygon)
    #We will iterate through all points and draw a vertical line at each x-position.
    #If this line cuts and more than two polygons are the result we will keep the cut.
    #Otherwise we will not cut at this point and continue with the next point.

    return_polygons = []
    isMonotone = True

    #we are splitting the list in half and combine them again so the cuts are more middled and therefore likely less cuts will be done.
    coords1, coords2 = split_list(polygon.exterior.coords)
    coords2.extend(coords1)

    for coord in coords2:
        current_x = coord[0]
        cutting_line = shapely.LineString([[current_x,min_y],[current_x,max_y]])
        cutted_geoms = shapely.ops.split(polygon, cutting_line).geoms
        if len(cutted_geoms) > 2:
            isMonotone = False
            for geom in cutted_geoms:
                if isinstance(geom, shapely.geometry.Polygon) and geom.area > 0:
                    return_polygons.extend(monotoneSubdivision(geom))
            break
    if isMonotone:
        return_polygons.append(polygon)
    
    return return_polygons


def earClippingTriangulate(geometry):
    vertices = list(geometry.exterior.coords)
    triangles = []

    if len(vertices) < 3:
        print("Unable to completely triangulate, less than three points given.")
        return triangles
    
    infinite_loop = False
    initial_vertices = vertices

    while( len(vertices) > 0):
        if infinite_loop == True:
            raise Exception("Unable to complete triangulation.")
        infinite_loop = True

        for current in range(len(vertices)):
            if len(vertices) < 3:
                if len(vertices) > 0:
                    vertices.pop(current)
                infinite_loop = False
                break

            #index of previous point
            test = len(vertices)
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


