import json
import re
import gradio as gr
from typing import List, Tuple, Union
import pandas as pd
import os
import shutil
import zipfile
from fastapi import FastAPI, Request
from app.statistics_shapefile import calculate_statistics_in_polygon
from dotenv import load_dotenv
import tempfile
import geopandas as gpd

def add_stats_to_dbf(dbf_path: str, stats: list, indice: str, output_dir: str) -> Tuple[str, str]:
    indice_dir = os.path.join(output_dir, indice)
    os.makedirs(indice_dir, exist_ok=True)

    gdf = gpd.read_file(dbf_path)
    for stat_dict in stats:
        for obj_id, stat_values in stat_dict.items():
            row_index = gdf[gdf['objectid'] == obj_id].index
            if not row_index.empty:
                for stat_name, value in stat_values.items():
                    if stat_name not in gdf.columns:
                        gdf[stat_name] = None
                    gdf.at[row_index[0], stat_name] = value

    updated_dbf_path = os.path.join(indice_dir, os.path.basename(dbf_path).replace(".dbf", f'_{indice}.dbf'))
    gdf.to_file(updated_dbf_path, driver="ESRI Shapefile")
    csv_path = os.path.join(indice_dir, os.path.basename(dbf_path).replace(".dbf", f'_{indice}.csv'))
    gdf.to_csv(csv_path, index=False)
    return updated_dbf_path, csv_path


def process_shp_data(images: List[str], shp: str) -> str:
    """
    Processes a shapefile and updates it with statistical data based on the provided images.

    Args:
        images (List[str]): A list of file paths for the images.
        shp (str): The file path for the shapefile (in .zip format).

    Returns:
        str: The path to the ZIP file containing updated shapefiles.
    """
    extract_path = tempfile.mkdtemp()  
    output_zip_path = os.path.join(tempfile.gettempdir(), 'updated_shapefile.zip')
    json_path = os.path.join(tempfile.gettempdir(), 'temp_shapefile.json')  
    indices = set()

    for image in images:
        match = re.match(r"(\w+)_\d{4}_\d{2}\.tif", os.path.basename(image))
        if match:
            indices.add(match.group(1))

    try:
        with zipfile.ZipFile(shp, 'r') as zip_ref:
            zip_ref.extractall(extract_path)

        shp_file = dbf_file = None
        for file in os.listdir(extract_path):
            if file.endswith('.shp'):
                shp_file = os.path.join(extract_path, file)
            elif file.endswith('.dbf'):
                dbf_file = os.path.join(extract_path, file)

        if not shp_file or not dbf_file:
            raise FileNotFoundError("No .shp or .dbf file found in ZIP.")

        gdf = gpd.read_file(shp_file)
        first_column_name = gdf.columns[0]
        gdf.to_file(json_path, driver='GeoJSON')

        with open(json_path) as f:
            geojson_data = json.load(f)
        for indice in indices:
            stats = []
            for feature in geojson_data['features']:
                geometry = feature['geometry']
                polygon_id = feature['properties'][first_column_name]
                stats.append(calculate_statistics_in_polygon(geometry, images, polygon_id, indice))
            indice_dir = os.path.join(extract_path, indice)
            os.makedirs(indice_dir, exist_ok=True)
            updated_dbf_path, csv_path = add_stats_to_dbf(dbf_file, stats, indice, extract_path)

        with zipfile.ZipFile(output_zip_path, 'w') as zipf:
            for indice in indices:
                indice_dir = os.path.join(extract_path, indice)
                for folder_name, _, filenames in os.walk(indice_dir):
                    for filename in filenames:
                        file_path = os.path.join(folder_name, filename)
                        zipf.write(file_path, os.path.relpath(file_path, extract_path))

        return output_zip_path

    finally:
        shutil.rmtree(extract_path)
        os.remove(json_path)


io =gr.Interface(
        fn=process_shp_data,
        theme="soft",
        inputs=[
            gr.File(label="Upload Images", file_count="multiple", type="filepath"),
            gr.File(label="Upload Shapefile ZIP", type="filepath"),   
        ],
        outputs=[
            gr.File(label="Download Updated Shapefiles")        ],
        title="Zonal Statistics",
        description="Calculate the statistics (mean value and standard deviation) of a given index over selected zones. The index is marked by an image and the zones are delineated by an overlapping polygon shapefile. Upload images and an overlapping zipped shapefile to calculate the statistics. The input image nomenclature must follow the format INDEX_YYYY_MM.TIF. The output file is an updated shapefile with the statistics as columns in the attribute table. ",
        flagging_mode="never",       
        analytics_enabled=False
    )