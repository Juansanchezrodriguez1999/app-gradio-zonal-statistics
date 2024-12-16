import os
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.io as pio
import geopandas as gpd
import plotly.graph_objs as go

def all_statistics(df, indice):
    """
    Processes a DBF file to compute the mean and standard deviation for a given index, returning the results as a DataFrame.
    Additionally, calculates the monthly averages across all years and returns a separate DataFrame with these values.

    Args:
        dbf_path (str): Path to the DBF file to be processed.
        indice (str): Index to analyze (e.g., 'NDVI').

    Returns:
        Tuple[pd.DataFrame, str, pd.DataFrame]: 
            - DataFrame containing the calculated statistics for each polygon.
            - Path to the generated CSV file containing the detailed results.
            - DataFrame containing the monthly averages across all years.
    """
    
    result = {'anio': [], 'mes': [], 'indice': [], 'media': [], 'mediana': [], 'desviacion': [], 'polygon_id': []}
    first_column_name = df.columns[0]
    for column in df.columns:
        if column.endswith('_mean'):
            anio_mes = column.split('_')[0]
            mes = anio_mes[:2]
            anio = anio_mes[2:]
            
            mean_vals = pd.to_numeric(df[column], errors='coerce')
            std_vals = pd.to_numeric(df[column.replace('_mean', '_std')], errors='coerce')
            median_vals = pd.to_numeric(df[column.replace('_mean', '_medi')], errors='coerce')
            
            result['anio'].extend([anio] * len(mean_vals))
            result['mes'].extend([mes] * len(mean_vals))
            result['indice'].extend([indice] * len(mean_vals))
            result['media'].extend(mean_vals)
            result['mediana'].extend(median_vals)
            result['desviacion'].extend(std_vals)
            result['polygon_id'].extend(df[first_column_name])

    df_result = pd.DataFrame(result)
    monthly_means = df_result.groupby('mes')['mediana'].mean().reset_index()
    monthly_means.columns = ['mes', 'mediana_temporal']
    
    filepath = "data.csv"
    df_result.to_csv(filepath, index=False)
    
    return df_result, filepath, monthly_means

def temporal_means(combined_df):
    """
    Calcula las medias de las medianas y las medias de las desviaciones 
    agrupando por mes y polygon_id.
    """
    combined_df["anio"] = combined_df["anio"].astype(str).str.zfill(2)  
    combined_df["mes"] = combined_df["mes"].astype(str).str.zfill(2)  

    combined_df["años"] = combined_df["anio"].astype(str)
    print(combined_df)
    result_df = (
        combined_df
        .groupby(["mes", "polygon_id"])
        .agg(
            media_mediana=("mediana", "mean"),
            media_desviacion=("desviacion", "mean"),
            años=("años", lambda x: f"{x.min()}-{x.max()}") 
        )
        .reset_index()
    )
    print(result_df)
    return result_df

def plot_statistics(statistics, plot_type, poligon):
    """
    Generates plots based on calculated statistics and saves them as HTML files.

    Args:
        statistics (pd.DataFrame): DataFrame containing statistical data (mean and standard deviation) by month and year.
        plot_types (list[str]): List of plot types to generate (options: "month", "total", "monthly").
        polygon (str): Identifier for the polygon being analyzed.

    Returns:
        list[str]: List of paths to the generated HTML files.
    """
    months_dict = {
        "01": "January", "02": "February", "03": "March", "04": "April", "05": "May",
        "06": "June", "07": "July", "08": "August", "09": "September", "10": "October",
        "11": "November", "12": "December"
    }

    output_dir = "plots"
    os.makedirs(output_dir, exist_ok=True)
    
    indexes_unicos = statistics["indice"].unique()
    archivos_html = []
    for index in indexes_unicos:
        df_filtered = statistics[statistics["indice"] == index]
        for grafica in plot_type:
            fig_medi = go.Figure()

            if grafica == "mes":
                for mes in sorted(set(months)):
                    df_mes = df_filtered[df_filtered["mes"] == mes]
                    years_mes = df_mes["anio"].tolist()
                    mes_means = df_mes["mediana"].tolist()
                    mes_std = df_mes["desviacion"].tolist()

                    fig_medi.add_trace(go.Scatter(x=years_mes, y=mes_means, mode='lines+markers', name=f"{months_dict[mes]} (mean)"))
                    fig_medi.add_trace(go.Scatter(x=years_mes, y=np.array(mes_means) + np.array(mes_std),
                                                  mode='lines', name=f"{months_dict[mes]} (mean + std)", line=dict(dash='dash')))
                    fig_medi.add_trace(go.Scatter(x=years_mes, y=np.array(mes_means) - np.array(mes_std),
                                                  mode='lines', name=f"{months_dict[mes]} (mean - std)", line=dict(dash='dash')))

            elif grafica == "total":
                years = df_filtered["anio"].tolist()
                months = df_filtered["mes"].tolist()
                means = df_filtered["mediana"].tolist()
                std = df_filtered["desviacion"].tolist()
                labels = [f"{a}-{months_dict[m]}" for a, m in zip(years, months)]
                
                fig_medi.add_trace(go.Scatter(x=labels, y=means, mode='lines+markers', name="median", marker=dict(color='green')))
                fig_medi.add_trace(go.Scatter(x=labels, y=np.array(means) + np.array(std),
                                              mode='lines', name="median + std", line=dict(dash='dash')))
                fig_medi.add_trace(go.Scatter(x=labels, y=np.array(means) - np.array(std),
                                              mode='lines', name="median - std", line=dict(dash='dash')))
                                              
            elif grafica == "temporal":
                monthly_mean = df_filtered.groupby("mes")["mediana"].mean()
                monthly_std = df_filtered.groupby("mes")["desviacion"].mean()  

                fig_medi.add_trace(go.Scatter(x=[months_dict[m] for m in monthly_mean.index], y=monthly_mean.values,
                                              mode='lines+markers', name="Median by Month", marker=dict(color='blue')))
                fig_medi.add_trace(go.Scatter(x=[months_dict[m] for m in monthly_mean.index], y=(monthly_mean + monthly_std.fillna(0)).values,
                                              mode='lines', name="Median + STD", line=dict(dash='dash')))
                fig_medi.add_trace(go.Scatter(x=[months_dict[m] for m in monthly_mean.index], y=(monthly_mean - monthly_std.fillna(0)).values,
                                              mode='lines', name="Median - STD", line=dict(dash='dash')))

            fig_medi.update_layout(
                title=f"{index} - {grafica} - {poligon} median-std",
                xaxis_title="Date (year-month)" if grafica != "temporal" else "Month",
                yaxis_title="Median value",
                hovermode="x unified"
            )

            html_filename = os.path.join(output_dir, f"plot_{index}_{grafica}_{poligon}.html")
            pio.write_html(fig_medi, html_filename)
            archivos_html.append(html_filename)

    return archivos_html
