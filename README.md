# Fork
This is a Fork from [https://github.com/nimbus-flight/osm-to-3dprint](https://github.com/nimbus-flight/osm-to-3dprint) with several adjustments and improvements.
The two major improvements are compatibility with current (December 2024) Python Libraries and Multi Color Export in multiple files.
Part of these improvements is the correct triangulation of shapes. In the original code the solution did not work correctly in a lot of situations
Since the original Version wasn't even working the sponsoring links are removed in this readme.

# osm-to-3dprint
Export OpenStreetMap (OSM) tiles, convert them to STL format, and import them directly into your 3D slicing software.
- building height data may be missing in certain cities, so check for other height attributes in the osm data
```# Check for various height attributes
    height_attrs = ['height', 'building:height', 'building:levels']
```

# Benefits
- Optimized for 3D Printing: The exported STL files are designed to have no non-manifold edges, which ensures that the model is ready for 3D printing without the need for repairs.
- Small File Size: Despite covering extensive areas, the exported STL files are compact in size. For example, the entire downtown area of San Francisco (buildings.stl) is around ~17.5 MB (around 8.86 square kilometers).
- Multi File Export of Water, Paths, Baseplate and Buildings for Multi Color Printing

# Features
- Customizable Parameters: Easily adjust parameters such as target size, maximum building height, and base thickness to suit your specific needs.
- Ready for Slicing: The generated STL files are optimized for 3D printing, ensuring minimal preparation time and reducing potential errors.

## Installation
```python3 -m pip install -r requirements.txt```

This Project relies on a C++ Algorythm for Triangulation of the areas. I could not find any working library for python. 
Please follow the instructions on GitHub for integration of the Tool. We expect it in the same directory as "main.py" with filename a.out.

This project relies on several libraries, including OSM, overpass-api, osmnx, shapely, trimesh, and numpy-stl.

## Sample Usage
The whole project needs to be consideres work in progress. Most of the times the main branch should work but please consider using older commits if not.

Adjust the GPS coordinates to define the bounding box of the area you wish to export. There are two options to do that, either use a center and a square size or define a box by coordinates.

#bbox = square_poly(Latitude, Longitude, Square Size in Meter)
bbox = square_poly(48.32950556656733, 10.90461275575229, 1000) #Augsburg Example
```
#bbox = min Longitude , min Latitude , max Longitude , max Latitude
bbox = (4.87123, 52.35893, 4.93389, 52.38351) # Amsterdam Example
```

You can also modify parameters like target_size, max_height, and base_thickness to customize the export.

Run the script:
```python3 main.py```

## Example Output
Here's are samples visualized in Bambu Studio:
Amsterdam:
![image](https://github.com/user-attachments/assets/94267752-b349-49f3-b246-a426d121780d)

Nürnberg Rangierbahnhof:
![image](https://github.com/user-attachments/assets/0635b480-f29c-4d2a-baa7-708164dce163)

For Bambu Studio select all STL File at once, then you will be asked if you want to import it as one object. Afterwards you can select multiple colors.

![Screenshot from 2024-12-27 23-30-40](https://github.com/user-attachments/assets/e4eb4981-c2ce-4eba-a010-1145df96813a)



## Acknowledgments / Licensing
The original Project is licensed with MIT License, so this one will be too. BUT the original project also had the following phrase regarding ChatGPT. Code out of Language Models is not license free and depending on situation might be copied from other projects. No new AI Code was added but the original parts remain. Please consider this when using the Code.
Most of the original code is replaced in the meantime. The ChatGPT Problem should be almost non-existent.



