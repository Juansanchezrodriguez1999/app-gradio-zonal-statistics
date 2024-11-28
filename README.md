# Zonal Statistics

## Description
This application processes ShapeFile files by calculating mean and median statistics on images based on the geometry of entities within the ShapeFile. The calculated statistics are added to the `.dbf` file of the ShapeFile and saved along with a CSV file.

Images should be in the format `index_year_month.tif` (for example, `ndvi_2023_10.tif`), where "index" is the name of the index to calculate (such as "ndvi").

## Features
1. Loads ShapeFile and multiple images.
2. Calculates mean and median statistics based on the ShapeFile geometry.
3. Creates specific folders for each index and stores `.dbf` and `.csv` files with the statistics.
4. Downloads processed files in a ZIP archive.

## Code Structure
- **add_stats_to_dbf**: Function that adds calculated statistics to a `.dbf` file and saves them along with a `.csv` file.
- **process_shp_data**: Main function that processes the ShapeFile and generates folders with statistics for each index.
- **validate_login**: Function to validate access credentials.
- **shapefile_interface**: User interface created with Gradio.

## Requirements
- Python 3.x
- Libraries: `gradio`, `fastapi`, `pandas`, `geopandas`, `dotenv`, `zipfile`, `tempfile`

## Execution
1. Install dependencies:
    ```bash
    pip install -r requirements.txt
    ```

2. Run the FastAPI server:
    ```bash
    uvicorn main:app --reload
    ```

3. Open the browser at `http://127.0.0.1:8000` to access the application.

## Usage
1. Upload the ShapeFile (in a ZIP file) and images in the format `index_year_month.tif`.
2. Select the image format (tif or jp2).
3. Download the ZIP file with processed data and statistics.