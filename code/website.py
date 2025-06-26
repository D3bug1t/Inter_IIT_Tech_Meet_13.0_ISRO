import numpy as np
import pandas as pd
import os 
import matplotlib.pyplot as plt
import dash_leaflet as dl
import dash_leaflet.express as dlx
from dash import Dash, html, Output, Input
from dash_extensions.javascript import arrow_function, assign
import json
import geopandas
from matplotlib.colors import Normalize
import matplotlib.colors as mcolors

url  = 'https://cartocdn-gusc.global.ssl.fastly.net/opmbuilder/api/v1/map/named/opm-moon-basemap-v0-1/all/{z}/{x}/{y}.png'
element =  input("Enter element to create interactive website: ")

gdf =  geopandas.read_file(f'shape_files/{element}Si_Ratio_Shpfile.shp')
gdf.set_crs('epsg:4326', inplace=True)
if not os.path.exists('GeoJson_files'):
    os.mkdir('GeoJson_files')

gdf.to_file(f'GeoJson_files/{element}.json')

with open(f'GeoJson_files/{element}.json') as f:
    red_data  = json.load(f)

vals=  [0]*len(red_data['features'])
for i in range(len(red_data['features'])):
    vals[i] =  red_data['features'][i]['properties'][f'Average_{element}']

Q1 = np.percentile(vals,25)
Q3 = np.percentile(vals, 75)
IQR = Q3 - Q1

    # Define bounds for outliers
lower_bound = Q1 - 1.5 * IQR
upper_bound = Q3 + 1.5 * IQR


norm  = Normalize(vmin = lower_bound, vmax=  upper_bound)
colors = [0]*len(red_data['features'])
for i in range(len(red_data['features'])):
    colors[i] = mcolors.rgb2hex(plt.cm.turbo(norm(red_data['features'][i]['properties'][f'Average_{element}'])))

def get_info(feature=None):
    header = [html.H4("Ratio")]
    if not feature:
        return header + [html.P("Hover over a pixel!")]
    return header + [
                     "Ratio : {:.3f}".format(feature["properties"][f"Average_{element}"])]

colorscale  =dict(zip(vals, colors))
style = dict(weight=1, opacity=0.7, color=None,edge_width = 0, fillOpacity=0.8, fillColor = 'red')

c=  [mcolors.rgb2hex(plt.cm.turbo(norm(i/100))) for i in range(int(100*round(lower_bound,1)),int(100*round(upper_bound,1)),10)]
colorbar = dl.Colorbar( colorscale=c, width=300, height=30, position="bottomleft", min=lower_bound, max=upper_bound)
# from branca.colormap import LinearColormap

# colormap = LinearColormap(["green", "yellow", "red"], vmin = 0.6, vmax = 1.2)
# Geojson rendering logic, must be JavaScript as it is executed in clientside.
style_handle = assign("""function(feature, context){
    const {classes, colorscale, style, colorProp} = context.hideout
    style.fillColor = colorscale[feature.properties[colorProp]];  // set the fill color according to the class    
    return style;
}""")
# Create geojson.
geojson = dl.GeoJSON(data=red_data,  # url to geojson file
                     style=style_handle,  # how to style each polygon
                     zoomToBounds=True,  # when true, zooms to bounds when data changes (e.g. on load)
                     zoomToBoundsOnClick=True,  # when true, zooms to bounds of feature (e.g. polygon) on click
                     hoverStyle=arrow_function(dict(weight=5, color='#666', dashArray='')),  # style applied on hover
                     hideout=dict(colorscale=colorscale, style=style, colorProp=f"Average_{element}"),
                     id="geojson")
# Create info control.
info = html.Div(children=get_info(), id="info", className="info",
                style={"position": "absolute", "top": "10px", "right": "10px", "zIndex": "1000"})
# Create app.
app = Dash(prevent_initial_callbacks=True)
app.layout = dl.Map(children=[
    dl.TileLayer(url = url, attribution='Moon',minZoom=2), geojson, colorbar, info
], style={'height': '100vh'}, center=[56, 10], zoom=6)

@app.callback(Output("info", "children"), Input("geojson", "hoverData"))
def info_hover(feature):
    return get_info(feature)

if __name__ == '__main__':
    app.run_server()