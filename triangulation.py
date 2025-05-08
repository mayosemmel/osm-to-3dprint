import shapely

def earClippingTriangulate(vertices):
    #make vertices writable
    vertices = list(vertices)
    triangles = []

    if len(vertices) < 3:
        print("Unable to completely triangulate, likely to be a 8 shape or self intersecting polygon")
        return triangles
    
    infinite_loop = False
    initial_vertices = vertices

    while( len(vertices) > 0):
        if infinite_loop == True:
            raise Exception("Unable to compelete triangulation")
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


