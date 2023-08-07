import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.colors import Normalize, LinearSegmentedColormap
from matplotlib.colors import to_rgba
import numpy as np
import matplotlib as mpl
from shapely.geometry import Point, LineString


# Modify the figure size and DPI to match the desired resolution (7680x4320)
mpl.rcParams['figure.figsize'] = (16, 9)  # 16 inches by 9 inches aspect ratio
mpl.rcParams['savefig.dpi'] = 432  # 432 DPI for 7680x4320 resolution

# Read the CSV data
data = pd.read_csv('./1950-2022_actual_tornadoes.csv')
data['date'] = pd.to_datetime(data['date'])  # Convert 'date' column to datetime format

# Create a list of geometries (either points or lines)
geometries = []
for i in range(len(data)):
    slon, slat = data['slon'][i], data['slat'][i]
    elon, elat = data['elon'][i], data['elat'][i]
    
    if elon == 0 and elat == 0:
        # If elon and elat are both zero, treat it as a point geometry
        geometries.append(Point(slon, slat))
    else:
        # Otherwise, treat it as a line geometry
        geometries.append(LineString([Point(slon, slat), Point(elon, elat)]))

# Create a GeoDataFrame with the new geometries
lines_gdf = gpd.GeoDataFrame(data, geometry=geometries, crs="EPSG:4326")

# Get unique dates from the data and count the number of tornadoes for each date
date_counts = data['date'].value_counts().sort_index()
unique_dates = pd.date_range(start=date_counts.index.min(), end=date_counts.index.max(), freq='D')
num_frames = len(unique_dates)

# Create the colormap for fading effect (121 colors for 120 days + 1 for the current day)
fade_colors = ['black'] + ['red', 'black']
cmap = LinearSegmentedColormap.from_list('fade_colormap', fade_colors, N=121)

# Create the base map
fig, (ax_map, ax_line) = plt.subplots(nrows=2, figsize=(10, 10), gridspec_kw={'height_ratios': [5, 2]})
ax_map.set_xlim(-125, -65)  # Adjust the extent based on your data.
ax_map.set_ylim(25, 50)

# Initialize the line graph
line, = ax_line.plot([], [], color='blue')
ax_line.set_xlim(0, len(unique_dates) - 1)
ax_line.set_xlabel('Days')
ax_line.set_ylabel('Tornado Count')

# Create a list to track tornado counts for each day
tornado_count_list = [0] * len(unique_dates)

# Create a list to track the Line2D objects for LineStrings
existing_lines = []

def interpolate_color(color1, color2, t):
    # Linear interpolation for each component of the RGBA tuple
    return tuple(c1 * (1 - t) + c2 * t for c1, c2 in zip(color1, color2))
    

def animate(i):
    current_date = unique_dates[i]
    print('Processing day:', current_date)  # Print the current day being processed
    ax_map.set_title('Tornado Tracks - Date: {}'.format(current_date))
    
    # print('Clear the previous scatter plots')
    # ax_map.collections.clear()
    
    # print('Get tornado tracks for the current day')
    current_date_data = lines_gdf[lines_gdf['date'] <= current_date]

    # print('Calculate days since tornado occurrence for the new tornadoes on the current day')
    new_tornadoes = current_date_data[current_date_data['date'] == current_date]
    days_since_tornado = (current_date - new_tornadoes['date']).dt.days

    # print('Normalize days_since_tornado to [0, 1] for linear interpolation')
    days_normalized = days_since_tornado / 120

    # print('Interpolate color values from red to black based on days_since_tornado')
    interpolated_colors = [interpolate_color((1.0, 0.0, 0.0, 1.0), (0.0, 0.0, 0.0, 1.0), t)
                           for t in days_normalized]

    # print('Plot the existing tornado tracks (all tornadoes up to the current day) with fading effect')
    existing_tornadoes_data = lines_gdf[lines_gdf['date'] < current_date]
    existing_tornadoes_days_since = (current_date - existing_tornadoes_data['date']).dt.days
    existing_tornadoes_days_normalized = existing_tornadoes_days_since / 120
    existing_tornadoes_colors = [interpolate_color((1.0, 0.0, 0.0, 1.0), (0.0, 0.0, 0.0, 1.0), t)
                                 for t in existing_tornadoes_days_normalized]

    # print('Ensure the colors are within the valid range [0, 1]')
    existing_tornadoes_colors = np.clip(existing_tornadoes_colors, 0, 1)

    # print('Plot the existing tornado tracks as lines or points based on geometries')
    for idx, row in existing_tornadoes_data.iterrows():
        color_idx = min(idx, len(existing_tornadoes_colors) - 1)  # Handle the case when idx exceeds colors length
        if isinstance(row['geometry'], LineString):
            ax_map.plot(*row['geometry'].xy,
                        color=existing_tornadoes_colors[color_idx], linewidth=row['mag'] * .5)
        elif isinstance(row['geometry'], Point):
            ax_map.scatter(*row['geometry'].xy,
                           color=existing_tornadoes_colors[color_idx], s=row['mag'] * 5)

    # print('Plot the new tornado tracks (tornadoes on the current day) with fading effect')
    if len(interpolated_colors) > 0:
        new_tornadoes_points = new_tornadoes[new_tornadoes.geometry.type == 'Point']
        new_tornadoes_lines = new_tornadoes[new_tornadoes.geometry.type == 'LineString']

        # print('Plot the new tornado points')
        ax_map.scatter(new_tornadoes_points.geometry.x, new_tornadoes_points.geometry.y,
                       color=interpolated_colors[-1], s=new_tornadoes_points['mag'] * 5)

        # print('Plot the new tornado lines')
        for idx, row in new_tornadoes_lines.iterrows():
            color_idx = min(i, len(interpolated_colors) - 1)  # Handle the case when i exceeds colors length
            ax_map.plot(*row['geometry'].xy, color=interpolated_colors[color_idx], linewidth=row['mag'] * .5)

    # print('Update tornado count list')
    tornado_count_list[i] = date_counts.get(current_date, 0)

    # print('Plot the line graph showing the number of tornadoes for each day')
    line.set_data(range(len(unique_dates)), tornado_count_list)
    line.set_linewidth(0.1666) 
    ax_line.relim()
    ax_line.autoscale_view()

# Create the animation
animation = FuncAnimation(fig=fig, func=animate, frames=num_frames, interval=1, repeat=False) 


# Save the animation as an mp4 video
output_file = "tornado_tracks.mp4"
animation.save(output_file, fps=60, extra_args=['-vcodec', 'libx265'])

# Show the animation in the interactive window (optional)
#plt.show()

