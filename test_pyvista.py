"""This module is for testing out the possibilities of pyvista."""

import pyvista as pv
import shapely
from functions import *

#U-Shape
#vertices_out = [(0, 0, 5), (10, 0, 0), (10, 10, 0), (0, 10, 5), (0, 6, 5), (6, 6, 2), (6, 4, 2), (0, 4, 5), (0, 0, 5)]
#vertices_in = [[(9, 9, 0.5), (8, 9, 1), (8, 8, 1), (9, 8, 0.5)]]
vertices_out = [(0, 0, 0), (10, 0, 0), (10, 10, 0), (0, 10, 0), (0, 6, 0), (6, 6, 0), (6, 4, 0), (0, 4, 0), (0, 0, 0)]
geometry = shapely.Polygon(vertices_out)

mesh = pv.PolyData()
for polygon in cut_polygon(geometry):

    coords = list(polygon.exterior.coords)
    length = len(coords)
    faces = [length]
    faces.extend(range(length))

    mesh = mesh.merge(pv.PolyData(coords, faces))


vertices_out = [(0, 0, 5), (10, 0, 5), (10, 10, 5), (0, 10, 5), (0, 6, 5), (6, 6, 5), (6, 4, 5), (0, 4, 5), (0, 0, 5)]
geometry = shapely.Polygon(vertices_out)

for polygon in cut_polygon(geometry):

    coords = list(polygon.exterior.coords)
    length = len(coords)
    faces = [length]
    faces.extend(range(length))

    mesh = mesh.merge(pv.PolyData(coords, faces))


mesh = mesh.triangulate()

plotter = pv.Plotter()
actor = plotter.add_mesh(mesh)
actor = plotter.show_bounds(
    grid='front',
    location='outer',
    all_edges=True
)
plotter.show()

#mesh.plot(show_edges=True, line_width=5)
