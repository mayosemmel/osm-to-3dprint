"""This module is for testing out the possibilities of pyvista."""

import pyvista as pv
import numpy as np

vertices = [(0, 0, 1), (10, 0, 1), (10, 10, 1), (0, 10, 1), (0, 8, 1), (8, 8, 1), (8, 2, 1), (0, 2, 1), (0, 0, 1)]
#points = np.array([[0, 0, 0], [1, 0, 0], [1, 0.5, 0], [0, 0.5, 0]])
# faces = np.hstack([[3, 0, 2, 3]])
faces = np.hstack([[9, 0, 1, 2, 3, 4, 5, 6, 7, 8]])
#faces = [(4, 0, 1, 2, 3), (3, 2, 4, 5)]
mesh = pv.PolyData(vertices, faces)
mesh.plot(show_edges=True, line_width=5)

mesh = mesh.triangulate()

plotter = pv.Plotter()
actor = plotter.add_mesh(mesh)
actor = plotter.show_bounds(
    grid='front',
    location='outer',
    all_edges=True
)
#plotter.show()

mesh.plot(show_edges=True, line_width=5)
