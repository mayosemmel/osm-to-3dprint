# Fork
This is a Fork from [https://github.com/nimbus-flight/osm-to-3dprint](https://github.com/nimbus-flight/osm-to-3dprint) with several adjustments and improvements.
The two major improvements are compatibility with current (December 2024) Python Libraries and Multi Color Export in multiple files.
Since the original Version wasn't even working the sponsoring links are removed in this readme.

# osm-to-3dprint
Export OpenStreetMap (OSM) tiles, convert them to STL format, and import them directly into your 3D slicing software.
- building height data may be missing in certain cities, so check for other height attributes in the osm data
```# Check for various height attributes
    height_attrs = ['height', 'building:height', 'building:levels']
```

# Benefits
- Larger Area Exports: Unlike Cadmapper, which limits free exports to 1 square kilometer, osm-to-3dprint allows you to export much larger areas without any cost.
- Optimized for 3D Printing: The exported STL files are designed to have no non-manifold edges, which ensures that the model is ready for 3D printing without the need for repairs.
- Small File Size: Despite covering extensive areas, the exported STL files are compact in size. For example, the entire downtown area of San Francisco (buildings.stl) is around ~17.5 MB (around 8.86 square kilometers).
- Relatively Fast - less than one minute to run main.py and get a city file for any given city.
- Multi File Export of Water, Paths, Baseplate and Buildings for Multi Color Printing

# Features
- Customizable Parameters: Easily adjust parameters such as target size, maximum building height, and base thickness to suit your specific needs.
- Ready for Slicing: The generated STL files are optimized for 3D printing, ensuring minimal preparation time and reducing potential errors.

## Installation
```python3 -m pip install -r requirements.txt```

This project relies on several libraries, including OSM, overpass-api, osmnx, shapely, trimesh, and numpy-stl.

## Sample Usage
Adjust the GPS coordinates to define the bounding box of the area you wish to export:

```
# Example bounding box: min Longitude , min Latitude , max Longitude , max Latitude 
bbox = (4.87123, 52.35893, 4.93389, 52.38351) # Amsterdam Example
```

You can also modify parameters like target_size, max_height, and base_thickness to customize the export.

Run the script:
```python3 main.py```

## Known Restrictions
The Generation of Streets is not working properly yet. Also Nature like stuff is not mapped correctly. Therefore small villages don't look good (yet).

## Example Output
Here's a sample of Amsterdam visualized in Bambu Studio:

![Screenshot from 2024-12-27 23-16-30](https://github.com/user-attachments/assets/a9a5ded4-53f0-4c70-ae7f-d6aba429af53)

For Bambu Studio select all STL File at once, then you will be asked if you want to import it as one object. Afterwards you can select multiple colors.

![Screenshot from 2024-12-27 23-30-40](https://github.com/user-attachments/assets/e4eb4981-c2ce-4eba-a010-1145df96813a)



## Acknowledgments / Licensing
The original Project is licensed with MIT License, so this one will be too. BUT the original project also had the following phrase regarding ChatGPT. Code out of Language Models is not license free and depending on situation might be copied from other projects. No new AI Code was added but the original parts remain. Please consider this when using the Code.

*Special thanks to [ChatGPT](https://www.openai.com/chatgpt) by OpenAI for assisting in the development of the codebase.*



