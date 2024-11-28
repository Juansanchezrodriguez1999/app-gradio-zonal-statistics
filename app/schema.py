schema = {
    "name": "agora-datalab/zonal_statistics",
    "description": """# 
Calculate the statistics (mean value and standard deviation) of a given index over selected zones. The index is marked by an image and the zones are delineated by an overlapping polygon shapefile. Upload images and an overlapping zipped shapefile to calculate the statistics. The input image nomenclature must follow the format INDEX_YYYY_MM.TIF. The output file is an updated shapefile with the statistics as columns in the attribute table.## Getting Started
You can access the platform using the credentials generated for you by the platform administrator.
## Contact
If you have any questions, please contact us at [support@agora-datalab.eu](mailto:support@agora-datalab.eu).""",
    "labels": ["web-application", "data-service", "shapefile", "images", "zonal-statistics"],
    "jsonforms:schema": {
        "type": "object",
        "properties": {
            "username": {"type": "string", "readOnly": True},
            "password": {"type": "string", "readOnly": True},
        },
    },
    "jsonforms:uischema": {
        "type": "VerticalLayout",
        "elements": [
            {"type": "Label", "text": "Credentials"},
            {"type": "Control", "scope": "#/properties/username", "label": "Username"},
            {"type": "Control", "scope": "#/properties/password", "label": "Password"},
        ],
    },
    "jsonforms:data": {"username": "", "password": ""},
    "embed": "",
}