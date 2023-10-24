# this code inputs csv_files and assigns colors, then visualizes results with umap, find umap.UMAP( for adjustable parameters
# multiple generations are usually neccesary for optimal visual appeal and informative reference
# also outputs csv file containing all dot coordinates, for occasional requirement
# this code needs manual tweakage for optimal output, see #s below
 
import os
import csv
import umap.umap_ as umap
import matplotlib.pyplot as plt
import numpy as np
# Define the list of CSV files, any number's fine # the legend in the file's name, so change your file name into the legend you want to visualize
csv_files = ['AA-usage_perced.csv','another_AA-usage_perced.csv','even_more_AA-usage_perced.csv']
# Define the color map for each file with transparency
file_colors = {}
# Assign colors to each file based on the order of input # could use hexa... number like #9EA839 # match the number of csv inputs!
colors = ['blue', 'red', '#61FFE7']
# Define transparency for the dots # match the number of csv inputs!( i can't see ----- i can see goooood   :   0 - 1 )
dot_alphas = ['0.2','0.6','1.0']

for i, file in enumerate(csv_files):
    file_colors[os.path.splitext(file)[0]] = (colors[i], dot_alphas[i])   # Use base file name as key instead of full file name
# Define dot size, font size, and text color
dot_size = 1
font_size = 5
text_color = 'black'
# Define whether to include text from the CSV file onto the visualized figure. if dot much, False is recommended cuz i can't see
include_text = False
#
data = []
labels = []
dot_names = []
dot_coordinates = []
dot_csv_files = []
label_index= []
for file in csv_files:
    with open(file, 'r') as csv_file:
        csv_reader = csv.reader(csv_file)
      # next(csv_reader)  # Skip header row if neccesary
        for row in csv_reader:
            if row:  
                labels.append(os.path.splitext(file)[0])  
                data.append([float(x) for x in row[1:]]) 
                dot_names.append(row[0])  
                dot_csv_files.append(file)  
data = np.array(data)
# Perform UMAP dimensionality reduction with metaparameters, this is where you need to know how UMAP works, see 10.1038/nbt.4314 for introduction, or just any youtuber's video
# Feel free to adjust the parameters or add more for your need: https://umap-learn.readthedocs.io/en/latest/parameters.html
reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, n_components=2,, spread=1.0)
embedding = reducer.fit_transform(data)
# Create legend
legends = []
for label in set(labels):
    if label in file_colors:
        color, dot_alpha = file_colors[label]
        legends.append(plt.scatter([], [], color=color, alpha= 0.9, label=label)) # alpha=float(dot_alpha) if want the legend to fit the dot alpha as well
# Plot the UMAP visualization
grid_interval = 1  # Interval between grid lines
grid_alpha = 0.6  # Transparency value for grid lines
line_width = 0.5  # line width
x_tick_fontsize = 3  # font size in the x and y axis
y_tick_fontsize = 5
plt.figure(figsize=(10, 8))
plt.grid(True, which='both', linestyle='-', linewidth=line_width, alpha=grid_alpha)  # Add grids
plt.xticks(np.arange(min(embedding[:,0]), max(embedding[:,0]), grid_interval), fontsize=x_tick_fontsize)
plt.yticks(np.arange(min(embedding[:,1]), max(embedding[:,1]), grid_interval), fontsize=y_tick_fontsize)
for i, label in enumerate(set(labels)):
    indices = [j for j, x in enumerate(labels) if x == label]
    if label in file_colors:  
        color, dot_alpha = file_colors[label]
    else:  
        color = 'gray'
    plt.scatter(embedding[indices, 0], embedding[indices, 1], color=color, alpha=float(dot_alpha), s=dot_size)
    if include_text:
        for index in indices:
            plt.text(embedding[index, 0], embedding[index, 1], dot_names[index], fontsize=font_size, color=text_color)
    for index in indices:
        dot_coordinates.append([dot_names[index], embedding[index, 0], embedding[index, 1]])  # Add dot coordinates
        label_index.append(label)
plt.title('UMAP Visualization')
plt.xlabel('UMAP Dimension 1')
plt.ylabel('UMAP Dimension 2')
# show the legends. if not needed, delete
plt.legend(handles=legends)
# save the plot # feel free to change format into svg, png, pdf dadada... # pump those dpi up man 
plt.savefig('Visualized_Result.jpg', format='jpg', dpi=2000)
plt.show()
# export
output_file = 'Visualized_Result_coordinates.csv'
with open(output_file, 'w', newline='') as csv_file:
    writer = csv.writer(csv_file)
    writer.writerow(['Dot Name', 'Original CSV File', 'UMAP Dimension 1', 'UMAP Dimension 2']) 
    for i in range(len(dot_names)):
        file_name = os.path.basename(dot_csv_files[i])  
        writer.writerow([dot_coordinates[i][0], label_index[i], dot_coordinates[i][1], dot_coordinates[i][2]])
