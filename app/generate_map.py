import folium
from shapely.geometry import shape, MultiPolygon
from shapely.ops import unary_union

def generate_map_from_geojson(geojson_data: dict) -> str:
    """
    Generates an interactive map from GeoJSON data and saves it as an HTML file.

    Args:
        geojson_data (dict): GeoJSON data to visualize.

    Returns:
        folium.Map: Folium map object with all geometries.
    """
    geometries = [shape(feature["geometry"]) for feature in geojson_data["features"]]
    unified_geometry = unary_union(geometries)
    bounds = unified_geometry.bounds  
    center = [(bounds[1] + bounds[3]) / 2, (bounds[0] + bounds[2]) / 2]
    m = folium.Map(location=center, zoom_start=15)
    m.fit_bounds([[bounds[1], bounds[0]], [bounds[3], bounds[2]]])
    folium.GeoJson(geojson_data, name="Geometries").add_to(m)
    return m