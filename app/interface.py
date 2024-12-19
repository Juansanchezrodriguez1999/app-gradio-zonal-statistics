from datetime import datetime
import json
import re
import gradio as gr
from typing import List, Tuple
import pandas as pd
import os
import shutil
import zipfile
from fastapi import FastAPI, Request
from app.sigpac_to_geometry import sigpac_to_geometry
from sigpac_tools.find import find_from_cadastral_registry
from app.plots import all_statistics
from app.plots import plot_statistics
from app.cut_from_geometry import cut_from_geometry

from app.plots import temporal_means

from dotenv import load_dotenv
import tempfile
import geopandas as gpd
from app.statistics_shapefile import calculate_statistics_in_polygon

from app.download_merge import download_tif_files
from app.get_tiles import get_tiles_polygons
from app.generate_map import generate_map_from_geojson
import geopandas as gpd
from shapely.geometry import shape
from gradio_calendar import Calendar

import folium

# Coordenadas aproximadas del centro de Andalucía
lat = 37.5443
lon = -4.7278
zoom_start = 7  # Nivel de zoom inicial

# Crear el mapa
initial_map = folium.Map(location=[lat, lon], zoom_start=zoom_start, tiles="OpenStreetMap")

'''def create_zip(file_paths: List[str]) -> str:
    """
    Creates a ZIP file containing the provided image files.

    Args:
        file_paths (List[str]): List of file paths to include in the zip file.

    Returns:
        str: Path to the created zip file.
    """
    try:
        output_zip_path = os.path.join(tempfile.gettempdir(), 'cropped_images.zip')
        with zipfile.ZipFile(output_zip_path, 'w') as zipf:
            for file in file_paths:
                zipf.write(file, os.path.basename(file))
        return output_zip_path
    except Exception as e:
        raise Exception(f"Error creating zip file: {str(e)}")

def process_sigpac_data(province: int, municipality: int, polygon: int, parcel: int, precinct: int, format: str, images: List[str]) -> str:
    """
    Processes images by cutting them according to SIGPAC geometry and returns a ZIP file with cropped images.

    Args:
        province (int): Province code.
        municipality (int): Municipality code.
        polygon (int): Polygon code.
        parcel (int): Parcel code.
        precinct (int): Precinct code.
        format (str): Output image format (e.g., 'tif', 'jp2').
        images (List[str]): List of image file paths to process.

    Returns:
        str: Path to the ZIP file containing the cropped images.
    """
    try:
        geometry, metadata = sigpac_to_geometry(province, municipality, polygon, parcel, precinct)
        geojson_data = {
            "type": "FeatureCollection",
            "features": [
                {
                    "type": "Feature",
                    "geometry": geometry,
                    "properties": metadata
                }
            ]
        }
        
        geojson_path = "geometry.json"
        with open(geojson_path, "w") as geojson_file:
            json.dump(geojson_data, geojson_file)

        cropped_images = cut_from_geometry(geometry, format, images)
        output_zip_path = create_zip(cropped_images)
        return output_zip_path,geojson_path
    except FileNotFoundError as e:
        raise FileNotFoundError(f"File not found: {str(e)}")
    except Exception as e:
        raise Exception(f"An error occurred: {str(e)}")'''

def process_catastral_data(catastral_registry: int, images: List[str]) -> Tuple[str, str]:
    """
    Processes images by cutting them according to SIGPAC geometry and returns a ZIP file with cropped images and geometry in GeoJSON format.

    Args:
        catastral_registry (int): Cadastral registry number.
        format (str): Output image format (e.g., 'tif', 'jp2').
        images (List[str]): List of image file paths to process.

    Returns:
        Tuple[str, str]: Paths to the ZIP file containing cropped images and the GeoJSON file with geometry.
    """
    unique_formats = list(set(f.split('.')[-1].lower() for f in images if isinstance(f, str) and '.' in f))
    if len(unique_formats) > 1:
        raise ValueError(f"Unsupported format. You must upload images in one unique format.")

    try:
        geometry, metadata = find_from_cadastral_registry(catastral_registry)
        geojson_data = {
            "type": "FeatureCollection",
            "features": [
                {
                    "type": "Feature",
                    "geometry": geometry,
                    "properties": metadata
                }
            ]
        }
        total_html_files = []
        total_geojson_files = []
        combined_df = pd.DataFrame()  
        temporal_df = pd.DataFrame()  
        geojson_path = "geometry.json"
        with open(geojson_path, "w") as geojson_file:
            json.dump(geojson_data, geojson_file)

        images_dir = cut_from_geometry(geometry, unique_formats[0], images)
        indices = list({file.split('/')[-1].split('_')[0].strip() for file in images_dir})
        for indice in indices:
            stats_index = {}
            stats = []
            for feature in geojson_data['features']:
                geometry = feature['geometry']
                polygon_id = str(catastral_registry)
                feature['objectID'] = polygon_id
                stats.append(calculate_statistics_in_polygon(geometry, images_dir, polygon_id, indice))
            
            stats_index[indice] = stats

            for key, polygons in stats_index.items():
                records = []
                for polygon in polygons:
                    for polygon_id, metrics in polygon.items():
                        record = {'polygon_id': polygon_id, 'indice': key} 
                        record.update(metrics)
                        records.append(record)

                df = pd.DataFrame(records)
                df_result, csv_path, monthly_means = all_statistics(df, indice)
                combined_df = pd.concat([combined_df, df_result], ignore_index=True)

            '''temporal_df = temporal_means(combined_df[combined_df["indice"] == indice])
            temporal_df["indice"] = indice 
            temporal_dict = temporal_df.groupby('polygon_id').apply(
                lambda x: {
                    f"{row['mes']}/{row['años']}".replace(" ", ""): {
                        'mean': row['media'],
                        'std': row['desviacion']
                    }
                    for _, row in x.iterrows()
                }
            ).to_dict()'''

            convinced_dict = combined_df[combined_df["indice"] == indice].groupby('polygon_id').apply(
                lambda x: {
                    f"{row['mes']}-{row['anio']}".replace(" ", ""): {
                        'median': row['mediana'],
                        'mean': row['media'],
                        'std': row['desviacion']
                    }
                    for _, row in x.iterrows()
                }
            ).to_dict()

            for feature in geojson_data['features']:
                '''if feature["objectID"] in temporal_dict:
                    feature["temporalStatistics"] = temporal_dict[feature["objectID"]]'''
                if feature["objectID"] in convinced_dict:
                    feature["zonalStatistics"] = convinced_dict[feature["objectID"]]

            updated_file_name = f'Geojson_actualizado_{indice}.geojson'
            with open(updated_file_name, 'w') as file:
                json.dump(geojson_data, file, indent=4)
            total_geojson_files.append(updated_file_name)

        unique_polygons = combined_df["polygon_id"].unique()
        for polygon in unique_polygons:
            df_polygon = combined_df[combined_df["polygon_id"] == polygon].drop(columns=["polygon_id"])
            total_html_files.extend(plot_statistics(df_polygon, ["zonal"], polygon))     

        zip_output_geojson = os.path.join(tempfile.mkdtemp(), 'Geojson_actualizado.zip')
        with zipfile.ZipFile(zip_output_geojson, 'w') as zipf:
            for geojson_file in total_geojson_files:
                zipf.write(geojson_file, os.path.basename(geojson_file))

        zip_output_plots = os.path.join(tempfile.mkdtemp(), 'Gráficos.zip')
        with zipfile.ZipFile(zip_output_plots, 'w') as zipf:
            for html_file in total_html_files:
                zipf.write(html_file, os.path.basename(html_file))
        return zip_output_plots,zip_output_geojson
    except FileNotFoundError as e:
        raise FileNotFoundError(f"File not found: {str(e)}")
    except Exception as e:
        raise Exception(f"An error occurred: {str(e)}")
    
def process_catastral_data_sentinel(catastral_registry: int, indexes: list, date_start: str, date_end:str) -> str:
    """
    Processes images by cutting them according to SIGPAC geometry and returns a ZIP file with cropped images and geometry in GeoJSON format.

    Args:
        catastral_registry (int): Cadastral registry number.
        format (str): Output image format (e.g., 'tif', 'jp2').
        images (List[str]): List of image file paths to process.

    Returns:
        Tuple[str, str]: Paths to the ZIP file containing cropped images and the GeoJSON file with geometry.
    """
    year_start = date_start.strftime("%Y")
    month_start = date_start.strftime("%B")
    year_end = date_end.strftime("%Y")
    month_end = date_end.strftime("%B")
    years = [year_start,year_end]
    months = [month_start,month_end]
    indexes = [indexes.upper()]    
    total_html_files = []
    total_geojson_files = []
    combined_df = pd.DataFrame()  
    temporal_df = pd.DataFrame()  

    geometry, metadata = find_from_cadastral_registry(catastral_registry)
    geojson_data = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "geometry": geometry,
                "properties": metadata
            }
        ]
    }
    features = geojson_data["features"]
    geometries = [shape(feature["geometry"]) for feature in features]
    properties = [feature["properties"] for feature in features]

    gdf = gpd.GeoDataFrame(properties, geometry=geometries)
    gdf = gdf.set_crs(geojson_data["features"][0]["geometry"].get("CRS", "")
)

    zones_utm = get_tiles_polygons(gdf)
    list_zones_utm = list(zones_utm)
    images_dir = download_tif_files(list_zones_utm, years, indexes, months) 
    unique_formats = list(set(f.split('.')[-1].lower() for f in images_dir if isinstance(f, str) and '.' in f))
    if len(unique_formats) > 1:
        raise ValueError(f"Unsupported format. You must upload images in one unique format.")

    main_map = generate_map_from_geojson(geojson_data)
        
    cropped_images = cut_from_geometry(geometry,unique_formats[0],images_dir)

    for indice in indexes:
        stats_index = {}
        stats = []
        
        for feature in geojson_data['features']:
            geometry = feature['geometry']
            polygon_id = str(catastral_registry)
            feature['objectID'] = polygon_id
            stats.append(calculate_statistics_in_polygon(geometry, images_dir, polygon_id, indice))
        
        stats_index[indice] = stats

        for key, polygons in stats_index.items():
            records = []
            for polygon in polygons:
                for polygon_id, metrics in polygon.items():
                    record = {'polygon_id': polygon_id, 'indice': key} 
                    record.update(metrics)
                    records.append(record)

            df = pd.DataFrame(records)
            df_result, csv_path, monthly_means = all_statistics(df, indice)
            combined_df = pd.concat([combined_df, df_result], ignore_index=True)

        '''temporal_df = temporal_means(combined_df[combined_df["indice"] == indice])
        temporal_df["indice"] = indice 
        temporal_dict = temporal_df.groupby('polygon_id').apply(
            lambda x: {
                f"{row['mes']}/{row['años']}".replace(" ", ""): {
                    'mean': row['media'],
                    'std': row['desviacion']
                }
                for _, row in x.iterrows()
            }
        ).to_dict()'''

        convinced_dict = combined_df[combined_df["indice"] == indice].groupby('polygon_id').apply(
            lambda x: {
                f"{row['mes']}-{row['anio']}".replace(" ", ""): {
                    'median': row['mediana'],
                    'mean': row['media'],
                    'std': row['desviacion']
                }
                for _, row in x.iterrows()
            }
        ).to_dict()

        for feature in geojson_data['features']:
            '''if feature["objectID"] in temporal_dict:
                feature["temporalStatistics"] = temporal_dict[feature["objectID"]]'''
            if feature["objectID"] in convinced_dict:
                feature["zonalStatistics"] = convinced_dict[feature["objectID"]]

        updated_file_name = f'Geojson_actualizado_{indice}.geojson'
        with open(updated_file_name, 'w') as file:
            json.dump(geojson_data, file, indent=4)
        total_geojson_files.append(updated_file_name)

    unique_polygons = combined_df["polygon_id"].unique()
    for polygon in unique_polygons:
        df_polygon = combined_df[combined_df["polygon_id"] == polygon].drop(columns=["polygon_id"])
        total_html_files.extend(plot_statistics(df_polygon, ["zonal"], polygon))     

    zip_output_geojson = os.path.join(tempfile.mkdtemp(), 'Geojson_actualizado.zip')
    with zipfile.ZipFile(zip_output_geojson, 'w') as zipf:
        for geojson_file in total_geojson_files:
            zipf.write(geojson_file, os.path.basename(geojson_file))

    zip_output_plots = os.path.join(tempfile.mkdtemp(), 'Gráficos.zip')
    with zipfile.ZipFile(zip_output_plots, 'w') as zipf:
        for html_file in total_html_files:
            zipf.write(html_file, os.path.basename(html_file))
    
    return zip_output_plots, zip_output_geojson, main_map._repr_html_()

def process_geojson_data(geojson: str, images: List[str]) -> str:
    """
    Processes images by cutting them according to the provided GeoJSON geometry and returns a ZIP file with cropped images.

    Args:
        geojson (str): Path to the GeoJSON file.
        format (str): Output image format (e.g., 'tif', 'jp2').
        images (List[str]): List of image file paths to process.

    Returns:
        str: Path to the ZIP file containing the cropped images.
    """
    indices = set()
    unique_formats = list(set(f.split('.')[-1].lower() for f in images if isinstance(f, str) and '.' in f))
    if len(unique_formats) > 1:
        raise ValueError(f"Unsupported format. You must upload images in one unique format.")

    total_html_files = []
    total_geojson_files = []
    combined_df = pd.DataFrame()  
    temporal_df = pd.DataFrame()  

    try:
        with open(geojson, 'r') as file:
            geojson_data = json.load(file)

        for image in images:
                match = re.match(r"(\w+)_\d{4}_\d{2}\.tif", os.path.basename(image))
                if match:
                    indices.add(match.group(1))

        for indice in indices:
            stats_index = {}
            stats = []
            for feature in geojson_data['features']:
                geometry = feature['geometry']
                polygon_id = datetime.now().strftime("%Y%m%d%H%M%S%f")[:-3]
                feature['objectID'] = polygon_id
                images_dir = cut_from_geometry(geometry, unique_formats[0], images)
                stats.append(calculate_statistics_in_polygon(geometry, images_dir, polygon_id, indice))
            
            stats_index[indice] = stats

            for key, polygons in stats_index.items():
                records = []
                for polygon in polygons:
                    for polygon_id, metrics in polygon.items():
                        record = {'polygon_id': polygon_id, 'indice': key} 
                        record.update(metrics)
                        records.append(record)

                df = pd.DataFrame(records)
                df_result, csv_path, monthly_means = all_statistics(df, indice)
                combined_df = pd.concat([combined_df, df_result], ignore_index=True)

            '''temporal_df = temporal_means(combined_df[combined_df["indice"] == indice])
            temporal_df["indice"] = indice 
            temporal_dict = temporal_df.groupby('polygon_id').apply(
                lambda x: {
                    f"{row['mes']}/{row['años']}".replace(" ", ""): {
                        'mean': row['media'],
                        'std': row['desviacion']
                    }
                    for _, row in x.iterrows()
                }
            ).to_dict()'''

            convinced_dict = combined_df[combined_df["indice"] == indice].groupby('polygon_id').apply(
                lambda x: {
                    f"{row['mes']}-{row['anio']}".replace(" ", ""): {
                        'median': row['mediana'],
                        'mean': row['media'],
                        'std': row['desviacion']
                    }
                    for _, row in x.iterrows()
                }
            ).to_dict()

            for feature in geojson_data['features']:
                '''if feature["objectID"] in temporal_dict:
                    feature["temporalStatistics"] = temporal_dict[feature["objectID"]]'''
                if feature["objectID"] in convinced_dict:
                    feature["zonalStatistics"] = convinced_dict[feature["objectID"]]

            updated_file_name = f'Geojson_actualizado_{indice}.geojson'
            with open(updated_file_name, 'w') as file:
                json.dump(geojson_data, file, indent=4)
            total_geojson_files.append(updated_file_name)

        unique_polygons = combined_df["polygon_id"].unique()
        for polygon in unique_polygons:
            df_polygon = combined_df[combined_df["polygon_id"] == polygon].drop(columns=["polygon_id"])
            total_html_files.extend(plot_statistics(df_polygon, ["zonal"], polygon))     

        zip_output_geojson = os.path.join(tempfile.mkdtemp(), 'Geojson_actualizado.zip')
        with zipfile.ZipFile(zip_output_geojson, 'w') as zipf:
            for geojson_file in total_geojson_files:
                zipf.write(geojson_file, os.path.basename(geojson_file))

        zip_output_plots = os.path.join(tempfile.mkdtemp(), 'Gráficos.zip')
        with zipfile.ZipFile(zip_output_plots, 'w') as zipf:
            for html_file in total_html_files:
                zipf.write(html_file, os.path.basename(html_file))
        return zip_output_plots,zip_output_geojson
    except FileNotFoundError as e:
        raise FileNotFoundError(f"File not found: {str(e)}")
    except Exception as e:
        raise Exception(f"An error occurred: {str(e)}")

def process_geojson_data_sentinel(geojson: dict, indexes: list, date_start: str, date_end:str) -> str:
    """
    Processes images based on GeoJSON and returns a ZIP file with cropped images.

    Args:
        geojson (dict): GeoJSON data.
        years (list): List of years for data.
        indexes (list): List of indexes to apply.
        months (list): List of months for data.

    Returns:
        str: Path to the ZIP file with cropped images.
    """
    year_start = date_start.strftime("%Y")
    month_start = date_start.strftime("%B")
    year_end = date_end.strftime("%Y")
    month_end = date_end.strftime("%B")
    years = [year_start,year_end]
    months = [month_start,month_end]
    indexes = [indexes.upper()]    
    total_html_files = []
    
    total_geojson_files = []
    combined_df = pd.DataFrame()  
    temporal_df = pd.DataFrame() 
    gdf = gpd.read_file(geojson)
    json_path = os.path.join(tempfile.gettempdir(), 'temp_shapefile.json')
    gdf.to_file(json_path, driver='GeoJSON')
    with open(json_path) as f:
        geojson_data = json.load(f)
    zones_utm = get_tiles_polygons(gdf)
    list_zones_utm = list(zones_utm)
    images_dir = download_tif_files(list_zones_utm, years, indexes, months) 
    unique_formats = list(set(f.split('.')[-1].lower() for f in images_dir if isinstance(f, str) and '.' in f))
    if len(unique_formats) > 1:
        raise ValueError(f"Unsupported format. You must upload images in one unique format.")

    main_map = generate_map_from_geojson(geojson_data)
    geometry= geojson_data['features'][0]['geometry']

    cropped_images = cut_from_geometry(geometry,unique_formats[0],images_dir)

    for indice in indexes:
        stats_index = {}
        stats = []
        
        for feature in geojson_data['features']:
            geometry = feature['geometry']
            polygon_id = datetime.now().strftime("%Y%m%d%H%M%S%f")[:-3]
            feature['objectID'] = polygon_id
            stats.append(calculate_statistics_in_polygon(geometry, images_dir, polygon_id, indice))
        
        stats_index[indice] = stats

        for key, polygons in stats_index.items():
            records = []
            for polygon in polygons:
                for polygon_id, metrics in polygon.items():
                    record = {'polygon_id': polygon_id, 'indice': key} 
                    record.update(metrics)
                    records.append(record)

            df = pd.DataFrame(records)
            df_result, csv_path, monthly_means = all_statistics(df, indice)
            combined_df = pd.concat([combined_df, df_result], ignore_index=True)

        '''temporal_df = temporal_means(combined_df[combined_df["indice"] == indice])
        temporal_df["indice"] = indice 
        temporal_dict = temporal_df.groupby('polygon_id').apply(
            lambda x: {
                f"{row['mes']}/{row['años']}".replace(" ", ""): {
                    'mean': row['media'],
                    'std': row['desviacion']
                }
                for _, row in x.iterrows()
            }
        ).to_dict()'''

        convinced_dict = combined_df[combined_df["indice"] == indice].groupby('polygon_id').apply(
            lambda x: {
                f"{row['mes']}-{row['anio']}".replace(" ", ""): {
                    'median': row['mediana'],
                    'mean': row['media'],
                    'std': row['desviacion']
                }
                for _, row in x.iterrows()
            }
        ).to_dict()

        for feature in geojson_data['features']:
            '''if feature["objectID"] in temporal_dict:
                feature["temporalStatistics"] = temporal_dict[feature["objectID"]]'''
            if feature["objectID"] in convinced_dict:
                feature["zonalStatistics"] = convinced_dict[feature["objectID"]]

        updated_file_name = f'Geojson_actualizado_{indice}.geojson'
        with open(updated_file_name, 'w') as file:
            json.dump(geojson_data, file, indent=4)
        total_geojson_files.append(updated_file_name)

    unique_polygons = combined_df["polygon_id"].unique()
    for polygon in unique_polygons:
        df_polygon = combined_df[combined_df["polygon_id"] == polygon].drop(columns=["polygon_id"])
        total_html_files.extend(plot_statistics(df_polygon, ["zonal"], polygon))     

    zip_output_geojson = os.path.join(tempfile.mkdtemp(), 'Geojson_actualizado.zip')
    with zipfile.ZipFile(zip_output_geojson, 'w') as zipf:
        for geojson_file in total_geojson_files:
            zipf.write(geojson_file, os.path.basename(geojson_file))

    zip_output_plots = os.path.join(tempfile.mkdtemp(), 'Gráficos.zip')
    with zipfile.ZipFile(zip_output_plots, 'w') as zipf:
        for html_file in total_html_files:
            zipf.write(html_file, os.path.basename(html_file))
    
    return zip_output_geojson, zip_output_plots, main_map._repr_html_()

    
    return output_zip_path,main_map._repr_html_()

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

def process_shp_data(shp: str,images: List[str] ) -> str:
    """
    Processes images by cutting them according to the provided shapefile geometry and returns a ZIP file with cropped images.

    Args:
        shp (str): Path to the shapefile ZIP.
        format (str): Output image format (e.g., 'tif', 'jp2').
        images (List[str]): List of image file paths to process.

    Returns:
        str: Path to the ZIP file containing the cropped images.
    """
    print("imagenbes cargadas")

    print(images)
    extract_path = './temp_extracted_files'
    output_zip_path = os.path.join(tempfile.gettempdir(), 'Shapefile_actualizado.zip')
    updated_shapefiles_path = tempfile.mkdtemp()  
    total_html_files = []
    total_geojson_files = []
    combined_df = pd.DataFrame()  
    temporal_df = pd.DataFrame() 
    json_path = './temp_shapefile.json'
    indices = set()
    for image in images:
        match = re.match(r"(\w+)_\d{4}_\d{2}\.tif", os.path.basename(image))
        if match:
            indices.add(match.group(1))

    if not os.path.exists(extract_path):
        os.makedirs(extract_path)
    
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
        
        if shp_file:
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
                    unique_formats = list(set(f.split('.')[-1].lower() for f in images if isinstance(f, str) and '.' in f))
                    if len(unique_formats) > 1:
                        raise ValueError(f"Unsupported format. You must upload images in one unique format.")

                    images_dir = cut_from_geometry(geometry, unique_formats[0], images)
                    print("despues de cortee")
                    print(images_dir)
                    stats.append(calculate_statistics_in_polygon(geometry, images_dir, polygon_id, indice))
                indice_dir = os.path.join(extract_path, indice)
                os.makedirs(indice_dir, exist_ok=True)
                updated_dbf_path, csv_path = add_stats_to_dbf(dbf_file, stats, indice, extract_path)

            for index_folder in os.listdir(extract_path):
                folder_path = os.path.join(extract_path, index_folder)
                updated_folder_path = os.path.join(updated_shapefiles_path, index_folder)
                os.makedirs(updated_folder_path, exist_ok=True)

                if os.path.isdir(folder_path):
                    dbf_files = [file for file in os.listdir(folder_path) if file.endswith('.dbf')]
                    if dbf_files:
                        dbf_path = os.path.join(folder_path, dbf_files[0])
                        dbf_df = gpd.read_file(dbf_path)
                        df_result, csv_path, monthly_means = all_statistics(dbf_df, index_folder)
                        combined_df = pd.concat([combined_df, df_result], ignore_index=True)
                        '''temporal_df = temporal_means(combined_df)
                        
                        for _, row in temporal_df.iterrows():
                            mean_col = f"{row['mes']}{row['años'].replace('-', '').replace(' ', '')}mean"
                            std_col = f"{row['mes']}{row['años'].replace('-', '').replace(' ', '')}std"
                            
                            if mean_col not in dbf_df.columns:
                                dbf_df[mean_col] = None
                            if std_col not in dbf_df.columns:
                                dbf_df[std_col] = None
                            
                            dbf_df.loc[dbf_df[first_column_name] == row['polygon_id'], mean_col] = row['media']
                            dbf_df.loc[dbf_df[first_column_name] == row['polygon_id'], std_col] = row['desviacion']
                        
                        updated_dbf_path = os.path.join(updated_folder_path, os.path.basename(dbf_path))
                        dbf_df.to_file(updated_dbf_path, driver="ESRI Shapefile")
                        
                        for file in os.listdir(folder_path):
                            if file.endswith(('.shp', '.shx', '.prj', '.cpg')): 
                                shutil.copy(os.path.join(folder_path, file), updated_folder_path)'''
            unique_polygons = combined_df["polygon_id"].unique()
            for polygon in unique_polygons:
                df_polygon = combined_df[combined_df["polygon_id"] == polygon].drop(columns=["polygon_id"])
                total_html_files.extend(plot_statistics(df_polygon, ["zonal"], polygon))
            zip_output_plots = os.path.join(tempfile.mkdtemp(), 'Gráficos.zip')
            with zipfile.ZipFile(zip_output_plots, 'w') as zipf:
                for html_file in total_html_files:
                    zipf.write(html_file, os.path.basename(html_file))
            with zipfile.ZipFile(output_zip_path, 'w') as zipf:
                for indice in indices:
                    indice_dir = os.path.join(extract_path, indice)
                    for folder_name, _, filenames in os.walk(indice_dir):
                        for filename in filenames:
                            file_path = os.path.join(folder_name, filename)
                            zipf.write(file_path, os.path.relpath(file_path, extract_path))
        else:
            raise FileNotFoundError("No .shp file found in ZIP.")
        shutil.rmtree(extract_path)
        os.remove(json_path)
        return output_zip_path,zip_output_plots

    except FileNotFoundError as e:
        raise FileNotFoundError(f"File not found: {str(e)}")
    except Exception as e:
        raise Exception(f"An error occurred: {str(e)}")
    
def process_shp_data_sentinel(shp: str, indexes: list, date_start: str, date_end:str) -> str:
    """
    Processes images by cutting them according to the provided shapefile geometry and returns a ZIP file with cropped images.

    Args:
        shp (str): Path to the shapefile ZIP.
        format (str): Output image format (e.g., 'tif', 'jp2').
        images (List[str]): List of image file paths to process.

    Returns:
        str: Path to the ZIP file containing the cropped images.
    """
    year_start = date_start.strftime("%Y")
    month_start = date_start.strftime("%B")
    year_end = date_end.strftime("%Y")
    month_end = date_end.strftime("%B")
    years = [year_start,year_end]
    months = [month_start,month_end]
    indexes = [indexes.upper()]    
    extract_path = tempfile.mkdtemp()
    output_zip_path = os.path.join(tempfile.gettempdir(), 'Shapefile_actualizado.zip')
    updated_shapefiles_path = tempfile.mkdtemp()  
    total_html_files = []
    total_geojson_files = []
    combined_df = pd.DataFrame()  
    temporal_df = pd.DataFrame() 
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
    gdf = gpd.read_file(shp_file)
    json_path = os.path.join(tempfile.gettempdir(), 'temp_shapefile.json')
    gdf.to_file(json_path, driver='GeoJSON')
    first_column_name = gdf.columns[0]
    with open(json_path) as f:
        geojson_data = json.load(f)
    zones_utm = get_tiles_polygons(gdf)
    list_zones_utm = list(zones_utm)
    images_dir = download_tif_files(list_zones_utm, years, indexes, months) 
    unique_formats = list(set(f.split('.')[-1].lower() for f in images_dir if isinstance(f, str) and '.' in f))
    if len(unique_formats) > 1:
        raise ValueError(f"Unsupported format. You must upload images in one unique format.")


    for indice in indexes:
        stats = []
        for feature in geojson_data['features']:
            geometry = feature['geometry']
            polygon_id = feature['properties'][first_column_name]
            stats.append(calculate_statistics_in_polygon(geometry, images_dir, polygon_id, indice))
        indice_dir = os.path.join(extract_path, indice)
        os.makedirs(indice_dir, exist_ok=True)
        updated_dbf_path, csv_path = add_stats_to_dbf(dbf_file, stats, indice, extract_path)

    for index_folder in os.listdir(extract_path):
        folder_path = os.path.join(extract_path, index_folder)
        updated_folder_path = os.path.join(updated_shapefiles_path, index_folder)
        os.makedirs(updated_folder_path, exist_ok=True)

        if os.path.isdir(folder_path):
            dbf_files = [file for file in os.listdir(folder_path) if file.endswith('.dbf')]
            if dbf_files:
                dbf_path = os.path.join(folder_path, dbf_files[0])
                dbf_df = gpd.read_file(dbf_path)
                df_result, csv_path, monthly_means = all_statistics(dbf_df, index_folder)
                combined_df = pd.concat([combined_df, df_result], ignore_index=True)
                '''temporal_df = temporal_means(combined_df)
                
                for _, row in temporal_df.iterrows():
                    mean_col = f"{row['mes']}{row['años'].replace('-', '').replace(' ', '')}mean"
                    std_col = f"{row['mes']}{row['años'].replace('-', '').replace(' ', '')}std"
                    
                    if mean_col not in dbf_df.columns:
                        dbf_df[mean_col] = None
                    if std_col not in dbf_df.columns:
                        dbf_df[std_col] = None
                    
                    dbf_df.loc[dbf_df[first_column_name] == row['polygon_id'], mean_col] = row['media']
                    dbf_df.loc[dbf_df[first_column_name] == row['polygon_id'], std_col] = row['desviacion']
                
                updated_dbf_path = os.path.join(updated_folder_path, os.path.basename(dbf_path))
                dbf_df.to_file(updated_dbf_path, driver="ESRI Shapefile")
                
                for file in os.listdir(folder_path):
                    if file.endswith(('.shp', '.shx', '.prj', '.cpg')): 
                        shutil.copy(os.path.join(folder_path, file), updated_folder_path)'''
    unique_polygons = combined_df["polygon_id"].unique()
    for polygon in unique_polygons:
        df_polygon = combined_df[combined_df["polygon_id"] == polygon].drop(columns=["polygon_id"])
        total_html_files.extend(plot_statistics(df_polygon, ["zonal"], polygon))
    zip_output_plots = os.path.join(tempfile.mkdtemp(), 'Gráficos.zip')
    with zipfile.ZipFile(zip_output_plots, 'w') as zipf:
        for html_file in total_html_files:
            zipf.write(html_file, os.path.basename(html_file))
    with zipfile.ZipFile(output_zip_path, 'w') as zipf:
        for indice in indexes:
            indice_dir = os.path.join(extract_path, indice)
            for folder_name, _, filenames in os.walk(indice_dir):
                for filename in filenames:
                    file_path = os.path.join(folder_name, filename)
                    zipf.write(file_path, os.path.relpath(file_path, extract_path))
    main_map = generate_map_from_geojson(geojson_data)
    
    return output_zip_path,zip_output_plots,main_map._repr_html_()
    
def process_csv_data(csv: str, images: List[str], latitude_column:str, longitude_column:str) -> str:
    """
    Processes images by cutting them according to the provided shapefile geometry and returns a ZIP file with cropped images.

    Args:
        shp (str): Path to the shapefile ZIP.
        format (str): Output image format (e.g., 'tif', 'jp2').
        images (List[str]): List of image file paths to process.

    Returns:
        str: Path to the ZIP file containing the cropped images.
    """    
    df = pd.read_csv(csv)
    unique_formats = list(set(f.split('.')[-1].lower() for f in images if isinstance(f, str) and '.' in f))
    if len(unique_formats) > 1:
        raise ValueError(f"Unsupported format. You must upload images in one unique format.")

    coordinates = df[['Y', 'X']].values.tolist()
    coordinates.append(coordinates[0])
    geojson_data = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [coordinates],
                },
                "properties": {}
            }
        ]
    }
    total_html_files = []
    total_geojson_files = []
    combined_df = pd.DataFrame()  
    temporal_df = pd.DataFrame()  
    geojson_path = "geometry.json"
    with open(geojson_path, "w") as geojson_file:
        json.dump(geojson_data, geojson_file)
    geometry= geojson_data['features'][0]['geometry']
    images_dir = cut_from_geometry(geometry, unique_formats[0], images)
    indices = list({file.split('/')[-1].split('_')[0].strip() for file in images_dir})
    for indice in indices:
        stats_index = {}
        stats = []
        
        for feature in geojson_data['features']:
            geometry = feature['geometry']
            polygon_id = datetime.now().strftime("%Y%m%d%H%M%S%f")[:-3]
            feature['objectID'] = polygon_id
            stats.append(calculate_statistics_in_polygon(geometry, images_dir, polygon_id, indice))
        
        stats_index[indice] = stats

        for key, polygons in stats_index.items():
            records = []
            for polygon in polygons:
                for polygon_id, metrics in polygon.items():
                    record = {'polygon_id': polygon_id, 'indice': key} 
                    record.update(metrics)
                    records.append(record)

            df = pd.DataFrame(records)
            df_result, csv_path, monthly_means = all_statistics(df, indice)
            combined_df = pd.concat([combined_df, df_result], ignore_index=True)

        '''temporal_df = temporal_means(combined_df[combined_df["indice"] == indice])
        temporal_df["indice"] = indice 
        temporal_dict = temporal_df.groupby('polygon_id').apply(
            lambda x: {
                f"{row['mes']}/{row['años']}".replace(" ", ""): {
                    'mean': row['media'],
                    'std': row['desviacion']
                }
                for _, row in x.iterrows()
            }
        ).to_dict()
'''
        convinced_dict = combined_df[combined_df["indice"] == indice].groupby('polygon_id').apply(
            lambda x: {
                f"{row['mes']}-{row['anio']}".replace(" ", ""): {
                    'median': row['mediana'],
                    'mean': row['media'],
                    'std': row['desviacion']
                }
                for _, row in x.iterrows()
            }
        ).to_dict()

        for feature in geojson_data['features']:
            '''if feature["objectID"] in temporal_dict:
                feature["temporalStatistics"] = temporal_dict[feature["objectID"]]'''
            if feature["objectID"] in convinced_dict:
                feature["zonalStatistics"] = convinced_dict[feature["objectID"]]

        updated_file_name = f'Geojson_actualizado_{indice}.geojson'
        with open(updated_file_name, 'w') as file:
            json.dump(geojson_data, file, indent=4)
        total_geojson_files.append(updated_file_name)

    unique_polygons = combined_df["polygon_id"].unique()
    for polygon in unique_polygons:
        df_polygon = combined_df[combined_df["polygon_id"] == polygon].drop(columns=["polygon_id"])
        total_html_files.extend(plot_statistics(df_polygon, ["zonal"], polygon))     

    zip_output_geojson = os.path.join(tempfile.mkdtemp(), 'Geojson_actualizado.zip')
    with zipfile.ZipFile(zip_output_geojson, 'w') as zipf:
        for geojson_file in total_geojson_files:
            zipf.write(geojson_file, os.path.basename(geojson_file))

    zip_output_plots = os.path.join(tempfile.mkdtemp(), 'Gráficos.zip')
    with zipfile.ZipFile(zip_output_plots, 'w') as zipf:
        for html_file in total_html_files:
            zipf.write(html_file, os.path.basename(html_file))
    return zip_output_plots,zip_output_geojson

def process_csv_data_sentinel(csv: str, indexes: list, date_start: str, date_end:str, latitude_column:str, longitude_column:str) -> str:
    """
    Processes images by cutting them according to the provided shapefile geometry and returns a ZIP file with cropped images.

    Args:
        shp (str): Path to the shapefile ZIP.
        format (str): Output image format (e.g., 'tif', 'jp2').
        images (List[str]): List of image file paths to process.

    Returns:
        str: Path to the ZIP file containing the cropped images.
    """
    year_start = date_start.strftime("%Y")
    month_start = date_start.strftime("%B")
    year_end = date_end.strftime("%Y")
    month_end = date_end.strftime("%B")
    years = [year_start,year_end]
    months = [month_start,month_end]
    indexes = [indexes.upper()]    
    format="tif"
    indexes = [x.upper() for x in indexes]  
    total_html_files = []
    total_geojson_files = []
    combined_df = pd.DataFrame()  
    temporal_df = pd.DataFrame() 

    df = pd.read_csv(csv)

    coordinates = df[['Y', 'X']].values.tolist()
    coordinates.append(coordinates[0])
    geojson_data = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [coordinates],
                },
                "properties": {}
            }
        ]
    }

    features = geojson_data["features"]
    geometries = [shape(feature["geometry"]) for feature in features]
    properties = [feature["properties"] for feature in features]

    gdf = gpd.GeoDataFrame(properties, geometry=geometries)
    gdf.set_crs("EPSG:4326", allow_override=True, inplace=True)

    zones_utm = get_tiles_polygons(gdf)
    list_zones_utm = list(zones_utm)
    images_dir = download_tif_files(list_zones_utm, years, indexes, months) 
    main_map = generate_map_from_geojson(geojson_data)
    geometry= geojson_data['features'][0]['geometry']

    cropped_images = cut_from_geometry(geometry,format,images_dir)

    for indice in indexes:
        stats_index = {}
        stats = []
        
        for feature in geojson_data['features']:
            geometry = feature['geometry']
            polygon_id = datetime.now().strftime("%Y%m%d%H%M%S%f")[:-3]
            feature['objectID'] = polygon_id
            stats.append(calculate_statistics_in_polygon(geometry, images_dir, polygon_id, indice))
        
        stats_index[indice] = stats

        for key, polygons in stats_index.items():
            records = []
            for polygon in polygons:
                for polygon_id, metrics in polygon.items():
                    record = {'polygon_id': polygon_id, 'indice': key} 
                    record.update(metrics)
                    records.append(record)

            df = pd.DataFrame(records)
            df_result, csv_path, monthly_means = all_statistics(df, indice)
            combined_df = pd.concat([combined_df, df_result], ignore_index=True)

        '''temporal_df = temporal_means(combined_df[combined_df["indice"] == indice])
        temporal_df["indice"] = indice 
        temporal_dict = temporal_df.groupby('polygon_id').apply(
            lambda x: {
                f"{row['mes']}/{row['años']}".replace(" ", ""): {
                    'mean': row['media'],
                    'std': row['desviacion']
                }
                for _, row in x.iterrows()
            }
        ).to_dict()'''

        convinced_dict = combined_df[combined_df["indice"] == indice].groupby('polygon_id').apply(
            lambda x: {
                f"{row['mes']}-{row['anio']}".replace(" ", ""): {
                    'median': row['mediana'],
                    'mean': row['media'],
                    'std': row['desviacion']
                }
                for _, row in x.iterrows()
            }
        ).to_dict()

        for feature in geojson_data['features']:
            '''if feature["objectID"] in temporal_dict:
                feature["temporalStatistics"] = temporal_dict[feature["objectID"]]'''
            if feature["objectID"] in convinced_dict:
                feature["zonalStatistics"] = convinced_dict[feature["objectID"]]

        updated_file_name = f'Geojson_actualizado_{indice}.geojson'
        with open(updated_file_name, 'w') as file:
            json.dump(geojson_data, file, indent=4)
        total_geojson_files.append(updated_file_name)

    unique_polygons = combined_df["polygon_id"].unique()
    for polygon in unique_polygons:
        df_polygon = combined_df[combined_df["polygon_id"] == polygon].drop(columns=["polygon_id"])
        total_html_files.extend(plot_statistics(df_polygon, ["zonal"], polygon))     

    zip_output_geojson = os.path.join(tempfile.mkdtemp(), 'Geojson_actualizado.zip')
    with zipfile.ZipFile(zip_output_geojson, 'w') as zipf:
        for geojson_file in total_geojson_files:
            zipf.write(geojson_file, os.path.basename(geojson_file))

    zip_output_plots = os.path.join(tempfile.mkdtemp(), 'Gráficos.zip')
    with zipfile.ZipFile(zip_output_plots, 'w') as zipf:
        for html_file in total_html_files:
            zipf.write(html_file, os.path.basename(html_file))
    
    return zip_output_geojson, zip_output_plots, main_map._repr_html_()

'''OLD FUNCTIONS
def geojson_interface() -> gr.Interface:
    """
    Creates the Gradio interface for processing images using a GeoJSON file.

    Returns:
        gr.Interface: The Gradio interface for processing images based on GeoJSON geometries.
    """
    return gr.Interface(
        fn=process_geojson_data,
        inputs=[
            gr.File(label="Upload Images", file_count="multiple", type="filepath"),
            gr.File(label="Upload Geojson", type="filepath"),   
        ],
        outputs=[gr.File(label="Download Updated geojson"),
                 gr.File(label="Output Plots")        
        ],        
        title="Zonal Statistics (User Images)",
        description="Calculate the statistics (mean value and standard deviation) of a given index over selected zones. The index is marked by an image and the zones are delineated by an overlapping polygon shapefile. Upload images and an overlapping zipped shapefile to calculate the statistics. The input image nomenclature must follow the format INDEX_YYYY_MM.TIF. The output file is an updated shapefile with the statistics as columns in the attribute table. ",
        flagging_mode="never",       
        analytics_enabled=False
    )

def geojson_interface_sentinel() -> gr.Interface:
    """
    Creates the Gradio interface for processing images using a GeoJSON file.

    Returns:
        gr.Interface: The Gradio interface for processing images based on GeoJSON geometries.
    """
    return gr.Interface(
        fn=process_geojson_data_sentinel,
       inputs=[
            gr.File(label="Geojson"),
            gr.CheckboxGroup(label="Indexes", choices=['moisture', 'ndvi', 'ndwi', 'ndsi', 'evi', 'osavi', 'evi2', 'ndre', 'ndyi', 'mndwi', 'bri', 'ri', 'bsi', 'cril']),
            gr.Radio(label="Start Year", choices=[2017,2018, 2019, 2020, 2021, 2022, 2023, 2024]),
            gr.Radio(label="Start Month", choices=["January", "February", "March", "April", "May", "June", "July","August", "September", "October", "November", "December"]),
            gr.Radio(label="End Year", choices=[2017,2018, 2019, 2020, 2021, 2022, 2023, 2024]),
            gr.Radio(label="End Month", choices=["January", "February", "March", "April", "May", "June", "July","August", "September", "October", "November", "December"]),

        ],
        outputs=[gr.File(label="Download Updated geometry"),
                 gr.File(label="Download Zonal Statistics Plots"),
            gr.HTML(label="Map view", elem_id="output-map", visible=True, value=initial_map._repr_html_())
],
        title="Zonal Statistics (Sentinel-2 Images)",
        description="Calculate the statistics (mean, median, and standard deviation) of various Sentinel-2 images for selected indices and input geometries. The input can include geometries with cadastral registry, SIGPAC codes, Shapefiles, or GeoJSON files, as well as the desired years and indices. The output is an updated geometry file that includes the computed statistics as new data in the attributes, along with visual plots for analysis.",
        flagging_mode="never",       
        analytics_enabled=False
    )

def shapefile_interface() -> gr.Interface:
    """
    Creates the Gradio interface for processing images using a Shapefile ZIP.

    Returns:
        gr.Interface: The Gradio interface for processing images based on Shapefile geometries.
    """
    return gr.Interface(
        fn=process_shp_data,
        inputs=[
            gr.File(label="Upload Images", file_count="multiple", type="filepath"),
            gr.File(label="Upload Shapefile ZIP", type="filepath"),   
        ],
        outputs=[gr.File(label="Download Updated Shapefiles"),
                 gr.File(label="Output Plots")        
        ],        
        title="Zonal Statistics (User Images)",
        description="Calculate the statistics (mean value and standard deviation) of a given index over selected zones. The index is marked by an image and the zones are delineated by an overlapping polygon shapefile. Upload images and an overlapping zipped shapefile to calculate the statistics. The input image nomenclature must follow the format INDEX_YYYY_MM.TIF. The output file is an updated shapefile with the statistics as columns in the attribute table. ",
        flagging_mode="never",       
        analytics_enabled=False
    )

def shapefile_interface_sentinel() -> gr.Interface:
    """
    Creates the Gradio interface for processing images using a Shapefile ZIP.

    Returns:
        gr.Interface: The Gradio interface for processing images based on Shapefile geometries.
    """
    return gr.Interface(
        fn=process_shp_data_sentinel,
        inputs=[
            gr.File(label="Shapefile"),
            gr.CheckboxGroup(label="Indexes", choices=['moisture', 'ndvi', 'ndwi', 'ndsi', 'evi', 'osavi', 'evi2', 'ndre', 'ndyi', 'mndwi', 'bri', 'ri', 'bsi', 'cril']),
            gr.Radio(label="Start Year", choices=[2017,2018, 2019, 2020, 2021, 2022, 2023, 2024]),
            gr.Radio(label="Start Month", choices=["January", "February", "March", "April", "May", "June", "July","August", "September", "October", "November", "December"]),
            gr.Radio(label="End Year", choices=[2017,2018, 2019, 2020, 2021, 2022, 2023, 2024]),
            gr.Radio(label="End Month", choices=["January", "February", "March", "April", "May", "June", "July","August", "September", "October", "November", "December"]),

        ],
        outputs=[gr.File(label="Download Updated Shapefile"),
                 gr.File(label="Download Zonal Statistics Plots"),
            gr.HTML(label="Map view", elem_id="output-map", visible=True, value=initial_map._repr_html_())
],
        title="Zonal Statistics (Sentinel-2 Images)",
        description="Calculate the statistics (mean value and standard deviation) of a given index over selected zones. The index is marked by an image and the zones are delineated by an overlapping polygon shapefile. Upload images and an overlapping zipped shapefile to calculate the statistics. The input image nomenclature must follow the format INDEX_YYYY_MM.TIF. The output file is an updated shapefile with the statistics as columns in the attribute table. ",
        flagging_mode="never",       
        analytics_enabled=False
    )

def sigpac_interface() -> gr.Interface:
    """
    Creates the Gradio interface for processing images using SIGPAC geometry data.

    Returns:
        gr.Interface: The Gradio interface for processing images using SIGPAC parcel data.
    """
    return gr.Interface(
        fn=process_sigpac_data,
        inputs=[
            gr.Number(label="Province", value=41),
            gr.Number(label="Municipality", value=55),
            gr.Number(label="Polygon", value=10),
            gr.Number(label="Parcel", value=11),
            gr.Number(label="Precinct", value=3),
            gr.Radio(choices=["tif", "jp2"], label="Image Format", value=""),
            gr.File(label="Upload Images", file_count="multiple", type="filepath")
        ],
        outputs=gr.File(label="Download ZIP with Cropped Images"),
        title="SIGPAC Image Cropping",
        description=(
            "Crop images using geographic parcel data from SIGPAC. "
            "Provide specific SIGPAC codes (Province, Municipality, Polygon, Parcel, and Precinct) "
            "to define the geometry of the area to crop. Upload one or more images to be processed "
            "and select the output format (`tif` or `jp2`). "
            "The interface retrieves the parcel geometry based on the provided SIGPAC codes, applies it to the images, "
            "and generates a ZIP file containing the cropped images, ready for download."
        ),
        flagging_mode="never",       
        analytics_enabled=False
    )

def catastral_interface() -> gr.Interface:
    """
    Creates the Gradio interface for processing images using SIGPAC geometry data.

    Returns:
        gr.Interface: The Gradio interface for processing images using SIGPAC parcel data.
    """
    return gr.Interface(
        fn=process_catastral_data,
        inputs=[
            gr.Text(label="Catastral Registry", value="41055A010000110003JL"),
            gr.File(label="Upload Images", file_count="multiple", type="filepath")
        ],
        outputs=[gr.File(label="Download plots"),
                gr.File(label="Output geometry")
        ],
        title="Zonal Statistics (User Images)",
        description="Calculate the statistics (mean, median, and standard deviation) of various Sentinel-2 images for selected indices and input geometries. The input can include geometries with cadastral registry, SIGPAC codes, Shapefiles, or GeoJSON files, as well as the desired years and indices. The output is an updated geometry file that includes the computed statistics as new data in the attributes, along with visual plots for analysis. The input image nomenclature must follow the format INDEX_YYYY_MM.TIF",
        flagging_mode="never",       
        analytics_enabled=False
        )

def catastral_interface_sentinel() -> gr.Interface:
    """
    Creates the Gradio interface for processing images using SIGPAC geometry data.

    Returns:
        gr.Interface: The Gradio interface for processing images using SIGPAC parcel data.
    """
    return gr.Interface(
        fn=process_catastral_data_sentinel,
        inputs=[
            gr.Text(label="Catastral Registry", value="41055A010000110003JL"),
            gr.CheckboxGroup(label="Indexes", choices=['moisture', 'ndvi', 'ndwi', 'ndsi', 'evi', 'osavi', 'evi2', 'ndre', 'ndyi', 'mndwi', 'bri', 'ri', 'bsi', 'cril']),
            gr.Radio(label="Start Year", choices=[2017,2018, 2019, 2020, 2021, 2022, 2023, 2024]),
            gr.Radio(label="Start Month", choices=["January", "February", "March", "April", "May", "June", "July","August", "September", "October", "November", "December"]),
            gr.Radio(label="End Year", choices=[2017,2018, 2019, 2020, 2021, 2022, 2023, 2024]),
            gr.Radio(label="End Month", choices=["January", "February", "March", "April", "May", "June", "July","August", "September", "October", "November", "December"]),

        ],
        outputs=[gr.File(label="Download Updated geometry"),
                 gr.File(label="Download Zonal Statistics Plots"),
            gr.HTML(label="Map view", elem_id="output-map", visible=True, value=initial_map._repr_html_())
],
        title="Zonal Statistics (Sentinel-2 Images)",
        description="Calculate the statistics (mean, median, and standard deviation) of various Sentinel-2 images for selected indices and input geometries. The input can include geometries with cadastral registry, SIGPAC codes, Shapefiles, or GeoJSON files, as well as the desired years and indices. The output is an updated geometry file that includes the computed statistics as new data in the attributes, along with visual plots for analysis.",
        flagging_mode="never",       
        analytics_enabled=False
    )

def csv_interface() -> gr.Interface:
    """
    Creates the Gradio interface for processing images using SIGPAC geometry data.

    Returns:
        gr.Interface: The Gradio interface for processing images using SIGPAC parcel data.
    """
    return gr.Interface(
        fn=process_csv_data,
        inputs=[
            gr.File(label="CSV"),
            gr.File(label="Upload Images", file_count="multiple", type="filepath")
        ],
        outputs=[gr.File(label="Download plots"),
                gr.File(label="Output geometry")
        ],
        title="Zonal Statistics (User Images)",
        description=(
    "Use this tool to calculate zonal statistics of satellite images for an area defined by a CSV file containing coordinates. "
    "The input image nomenclature must follow the format INDEX_YYYY_MM.TIF"
    "The CSV file must include the following columns: Vertex, X, and Y, where X and Y are the longitude and latitude "
    "coordinates, respectively.\n Example format:"
    "<table style='font-size:12px; border-collapse:collapse; width:20%; margin:auto;'>"
    "<thead>"
    "<tr style='border-bottom:1px solid #ddd;'><th>Vertex</th><th>X</th><th>Y</th></tr>"
    "</thead>"
    "<tbody>"
    "<tr><td>1</td><td>37.66314226</td><td>-5.302012025</td></tr>"
    "<tr><td>2</td><td>37.66634019</td><td>-5.22130994</td></tr>"
    "<tr><td>3</td><td>37.72700059</td><td>-5.246265025</td></tr>"
    "<tr><td>4</td><td>37.72265957</td><td>-5.334305567</td></tr>"
    "</tbody>"
    "</table>"
),
        flagging_mode="never",       
        analytics_enabled=False
    )

def csv_interface_sentinel() -> gr.Interface:
    """
    Creates the Gradio interface for processing images using SIGPAC geometry data.

    Returns:
        gr.Interface: The Gradio interface for processing images using SIGPAC parcel data.
    """
    return gr.Interface(
        fn=process_csv_data_sentinel,
        inputs=[
            gr.File(label="CSV"),
            gr.CheckboxGroup(label="Indexes", choices=['moisture', 'ndvi', 'ndwi', 'ndsi', 'evi', 'osavi', 'evi2', 'ndre', 'ndyi', 'mndwi', 'bri', 'ri', 'bsi', 'cril']),
            gr.Radio(label="Start Year", choices=[2017,2018, 2019, 2020, 2021, 2022, 2023, 2024]),
            gr.Radio(label="Start Month", choices=["January", "February", "March", "April", "May", "June", "July","August", "September", "October", "November", "December"]),
            gr.Radio(label="End Year", choices=[2017,2018, 2019, 2020, 2021, 2022, 2023, 2024]),
            gr.Radio(label="End Month", choices=["January", "February", "March", "April", "May", "June", "July","August", "September", "October", "November", "December"]),
        ],
        outputs=[gr.File(label="Download Updated geometry"),
                 gr.File(label="Download Zonal Statistics Plots"),
            gr.HTML(label="Map view", elem_id="output-map", visible=True, value=initial_map._repr_html_())
],
        title="Zonal Statistics (Sentinel-2 Images)",
        description=(
    "Use this tool to calculate zonal statistics of satellite images for an area defined by a CSV file containing coordinates. "
    "The CSV file must include the following columns: Vertex, X, and Y, where X and Y are the longitude and latitude "
    "coordinates, respectively.\n Example format:"
    "<table style='font-size:12px; border-collapse:collapse; width:20%; margin:auto;'>"
    "<thead>"
    "<tr style='border-bottom:1px solid #ddd;'><th>Vertex</th><th>X</th><th>Y</th></tr>"
    "</thead>"
    "<tbody>"
    "<tr><td>1</td><td>37.66314226</td><td>-5.302012025</td></tr>"
    "<tr><td>2</td><td>37.66634019</td><td>-5.22130994</td></tr>"
    "<tr><td>3</td><td>37.72700059</td><td>-5.246265025</td></tr>"
    "<tr><td>4</td><td>37.72265957</td><td>-5.334305567</td></tr>"
    "</tbody>"
    "</table>"
    
),
        flagging_mode="never",       
        analytics_enabled=False
    )

def build_interface() -> gr.Blocks:
    """
    Builds the Gradio interface with separate tabs for SIGPAC and GeoJSON image processing.

    Returns:
        gr.Blocks: The complete Gradio interface with tabs.
    """
    with gr.Blocks(theme="soft") as interface:
        gr.Markdown("# Zonal Statistics")
        gr.Markdown("Zonal statistics general description")
        with gr.Tab("Catastral Registry User"):
            catastral_interface()
        with gr.Tab("Catastral Registry S-2"):
            catastral_interface_sentinel()
        with gr.Tab("Shapefile User"):
            shapefile_interface()
        with gr.Tab("Shapefile S-2"):
            shapefile_interface_sentinel()
        with gr.Tab("SIGPAC"):
            sigpac_interface()
        with gr.Tab("CSV User"):
            csv_interface()
        with gr.Tab("CSV S-2"):
            csv_interface_sentinel()
        with gr.Tab("GeoJSON User"):
            geojson_interface()
        with gr.Tab("GeoJSON S-2"):
            geojson_interface_sentinel()
        
    return interface

io = build_interface()'''

def build_interface() -> gr.Blocks:
    """
    Builds the Gradio interface with separate tabs for SIGPAC and GeoJSON image processing.

    Returns:
        gr.Blocks: The complete Gradio interface with tabs.
    """
    with gr.Blocks(theme="soft") as interface:  
        gr.Markdown("""# Recorte de imágenes
        Esta herramienta sirve para el cálculo de las estadisticas zonales de imágenes satélite. El cálculo se hace sobre la zona de la imágen que está dentro de la geometría proporcionada. 
        Esta herramienta pretende calcular las estadísticas zonales de las imágenes de un sensor específico (satélite, avión, dron) tomadas en un corto periodo de tiempo para obtener unas gráficas y una geometría actualizada con los valores estadísticos.
        Puedes elegir entre hacer las estadícticas de tus propias imágenes o directamente obtener las estadísticas de imágenes de archivo procedentes de capturas de Sentinel-2 indicando la zona deseada y el índice espectral.
        ### Selecciona imágenes fuente de imagen:
        """)
        with gr.Tab("Imágenes propias"):
            with gr.Row():
                with gr.Column():
                    gr.Markdown("""
                    ### Detalles de los inputs                
                    - **Imágenes**: Imágenes que el usuario desea recortar.
                    - **Tipo de geometría**: Indica cual es el archivo de geometría que se va a utilizar.
                    - **Archivo de geometría**: Archivo que contiene la geometría deseada.
                        - **Si seleccionas un archivo CSV**:
                            El archivo CSV debe contener 3 o más puntos de coordenada (lat/lon) para definir una geometría.  
                            En el CSV de ejemplo encontramos estas columnas:
                            - **Vertex**: Índice de fila (no necesario).
                            - **X**: Latitud (coordenadas geográficas).
                            - **Y**: Longitud (coordenadas geográficas).
                            ```csv
                            Vertex,X,Y
                            1,37.66314226,-5.302012025
                            2,37.66634019,-5.22130994
                            3,37.72700059,-5.246265025
                            4,37.72265957,-5.334305567
                            ```                    
                    - **Columna latitud**: Nombre de la columna de latitud en el csv.
                    - **Columna longitud**: Nombre de la columna de longitud en el csv.
                    - **Registro catastral**: Codigo catastral de la zona.
                    ### Detalles de los outputs                
                    - **Gráficos**: Gráficos con las estadísticas zonales de las imágenes y geometrías seleccionadas.
                    - **Geometría**: Geometría actualizada con los datos de las estadísticas zonales.
                    """)

                with gr.Column():
                    imagenes= gr.File(label="Subir imágenes", file_count="multiple", type="filepath")
                    input_file_type = gr.Radio(
                        label="Selecciona la geometría", 
                        choices=['Shapefile', 'CSV', 'Geojson', 'Registro Catastral']
                    )
                    geometry_file = gr.File(label="Sube archivo de geometría", visible=True)

                    latitude_column = gr.Text(label="Inserte el nombre de la columna de la latitud", value="X", visible=False)

                    longitude_column = gr.Text(label="Inserte el nombre de la columna de la longitud", value="Y", visible=False)

                    catastral_registry = gr.Text(label="Inserte código de registro catastral", value="41055A010000110003JL", visible=False)
                
            with gr.Row(): 
                submit_button = gr.Button(value="Procesar")  
                clear_button = gr.Button(value="Limpiar") 
            with gr.Row(): 
                output_plots = gr.File(label="Descargar gráficas")
                output_geometry = gr.File(label="Descargar geometría")
            def clear_inputs():
                """
                Resetea todos los campos del formulario a sus valores iniciales.
                """
                return (
                    gr.update(value=None),
                    gr.update(value=None),  
                    gr.update(value=None),  
                    gr.update(value=""),    
                    gr.update(value=""),    
                    gr.update(value=""),    
                    gr.update(value=None),  
                )
            clear_button.click(
                fn=clear_inputs,
                inputs=[],
                outputs=[
                    imagenes,              
                    input_file_type,        
                    geometry_file,         
                    latitude_column,    
                    longitude_column,       
                    catastral_registry,    
                ]
            )
            def update_visibility(file_type):
                if file_type == "Registro Catastral":
                    return (
                        gr.update(visible=False),  
                        gr.update(visible=True),   
                        gr.update(visible=False),  
                        gr.update(visible=False)   
                    )
                elif file_type == "CSV":
                    return (
                        gr.update(visible=True),   
                        gr.update(visible=False), 
                        gr.update(visible=True),   
                        gr.update(visible=True)    
                    )
                else:
                    return (
                        gr.update(visible=True),   
                        gr.update(visible=False),  
                        gr.update(visible=False), 
                        gr.update(visible=False)  
                    )

            def process_inputs(imagenes,input_file_type, geometry_file, catastral_registry, latitude_column, longitude_column):
                if input_file_type == "Registro Catastral":
                    return process_catastral_data(catastral_registry, imagenes)
                elif input_file_type == "Shapefile":
                    return process_shp_data(geometry_file, imagenes )
                elif input_file_type == "CSV":
                    return process_csv_data(geometry_file, imagenes, latitude_column, longitude_column)
                else:
                    return process_geojson_data(geometry_file, imagenes)

            input_file_type.change(
                fn=update_visibility,
                inputs=[input_file_type],
                outputs=[geometry_file, catastral_registry, latitude_column, longitude_column]
            )

            submit_button.click(
                fn=process_inputs,
                inputs=[imagenes, input_file_type, geometry_file, catastral_registry, latitude_column, longitude_column],
                outputs=[output_plots,output_geometry]
            )
    
        with gr.Tab("Imágenes de archivo"):
            with gr.Row():
                with gr.Column():
                    gr.Markdown("""
                    ### Detalles de los inputs                
                    - **Índice**: Índice espectral.
                    - **Tipo de geometría**: Indica cual es el archivo de geometría que se va a utilizar.
                        - **Archivo de geometría**: Archivo que contiene la geometría deseada.
                        - **Si seleccionas un archivo CSV**:
                            El archivo CSV debe contener 3 o más puntos de coordenada (lat/lon) para definir una geometría.  
                            En el CSV de ejemplo encontramos estas columnas:
                            - **Vertex**: Índice de fila (no necesario).
                            - **X**: Latitud (coordenadas geográficas).
                            - **Y**: Longitud (coordenadas geográficas).
                            ```csv
                            Vertex,X,Y
                            1,37.66314226,-5.302012025
                            2,37.66634019,-5.22130994
                            3,37.72700059,-5.246265025
                            4,37.72265957,-5.334305567
                            ```    
                    - **Columna latitud**: Nombre de la columna de latitud en el csv.
                    - **Columna longitud**: Nombre de la columna de longitud en el csv.
                    - **Registro catastral**: Codigo catastral de la zona.
                    - **Fecha**: Fecha de la imágen (Los satélites no pasan todos los días por lo que se buscará la fecha más cercana a la indicada).
                    ### Detalles de los outputs                
                    - **Gráficos**: Gráficos con las estadísticas zonales de las imágenes y geometrías seleccionadas.
                    - **Geometría**: Geometría actualizada con los datos de las estadísticas zonales.
                    - **Mapa**: Mapa con la geometría recortada.
                    """)

                with gr.Column():
                    indexes = gr.Radio(
                        label="Selecciona el índice", 
                        choices=['moisture', 'ndvi', 'ndwi', 'ndsi', 'evi', 'osavi', 'evi2', 'ndre', 'ndyi', 'mndwi', 'bri', 'ri', 'bsi', 'cril']
                    )
                    input_file_type = gr.Radio(
                        label="Selecciona la geometría", 
                        choices=['Shapefile', 'CSV', 'Geojson', 'Registro Catastral']
                    )
                    geometry_file = gr.File(label="Sube archivo de geometría", visible=True)
                    latitude_column = gr.Text(label="Inserte el nombre de la columna de la latitud", value="X", visible=False)

                    longitude_column = gr.Text(label="Inserte el nombre de la columna de la longitud", value="Y", visible=False)


                    catastral_registry = gr.Text(label="Inserte código de registro catastral", value="41055A010000110003JL", visible=False)

                    selected_date_start = Calendar(type="date", label="Seleccionar fecha inicial")
                    selected_date_end = Calendar(type="date", label="Seleccionar fecha final")
                
                      
                    
            with gr.Row():
                submit_button = gr.Button(value="Procesar")
                clear_button = gr.Button(value="Limpiar") 
            with gr.Row():
                output_plots = gr.File(label="Descargar gráficas")
                output_geometry = gr.File(label="Descargar geometría")
                map_view = gr.HTML(label="Mapa", elem_id="output-map", value=initial_map._repr_html_(), visible=True)

            def clear_inputs():
                """
                Resetea todos los campos del formulario a sus valores iniciales.
                """
                return (
                    gr.update(value=None),
                    gr.update(value=None),  
                    gr.update(value=None),  
                    gr.update(value=""),    
                    gr.update(value=""),    
                    gr.update(value=""),    
                    gr.update(value=None),  
                    gr.update(value=None),  
                    gr.update(value=None),  
                )

            clear_button.click(
                fn=clear_inputs,
                inputs=[],
                outputs=[
                    indexes,              
                    input_file_type,        
                    geometry_file,         
                    latitude_column,    
                    longitude_column,       
                    catastral_registry,    
                    selected_date_start,          
                    selected_date_end,          
                ]
            )
            def update_visibility(file_type):
                if file_type == "Registro Catastral":
                    return (
                        gr.update(visible=False),  
                        gr.update(visible=True),   
                        gr.update(visible=False),  
                        gr.update(visible=False)   
                    )
                elif file_type == "CSV":
                    return (
                        gr.update(visible=True),   
                        gr.update(visible=False), 
                        gr.update(visible=True),   
                        gr.update(visible=True)    
                    )
                else:
                    return (
                        gr.update(visible=True),   
                        gr.update(visible=False),  
                        gr.update(visible=False), 
                        gr.update(visible=False)  
                    )

            def process_inputs(input_file_type, geometry_file, catastral_registry, indexes, selected_date_start, selected_date_end, latitude_column, longitude_column):
                if input_file_type == "Registro Catastral":
                    return process_catastral_data_sentinel(catastral_registry,indexes, selected_date_start, selected_date_end)
                elif input_file_type == "Shapefile":
                    return process_shp_data_sentinel(geometry_file, indexes, selected_date_start, selected_date_end)
                elif input_file_type == "CSV":
                    return process_csv_data_sentinel(geometry_file, indexes, selected_date_start, selected_date_end, latitude_column, longitude_column)
                else:
                    return process_geojson_data_sentinel(geometry_file, indexes, selected_date_start, selected_date_end)

            input_file_type.change(
                fn=update_visibility,
                inputs=[input_file_type],
                outputs=[geometry_file, catastral_registry, latitude_column, longitude_column]
            )

            submit_button.click(
                fn=process_inputs,
                inputs=[input_file_type, geometry_file, catastral_registry, indexes, selected_date_start, selected_date_end, latitude_column, longitude_column],
                outputs=[output_plots,output_geometry, map_view]
            )
    
    return interface

io = build_interface()