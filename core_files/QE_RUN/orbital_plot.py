import numpy as np
import pandas as pd
import os
import cv2
import matplotlib.pyplot as plt
import matplotlib as mpl
import file_analyzer as fa

from functools import cmp_to_key
from PIL import Image

csfont = {'fontname': 'Arial Narrow'}

fig = plt.figure(figsize=(30, 19))
plt.rcParams.update({'font.size': 50})
plt.rcParams["font.family"] = "DejaVu Sans"
mpl.rcParams['axes.linewidth'] = 0.7  # set the value globally
plt.rcParams['text.antialiased'] = True

################### 1st Plot#########

ax1 = plt.subplot(121)  # No. of subplotscolo: Read this as 1 (rows) by 2 (coulomns) plot position number 1.

# Taking required inputs from users
# Set the path to the directory where all your pdos, .gnu, .out and bands.in file are stored
root_dir = './orbital_data'
#print('root_dir = ', root_dir)
#print("Please enter the energy range in which you want to analyze the orbitals contribution")




LCB =  5	#LCB = float(input("Enter Lower Conduction Band energy in float datatype:"))

UVB = -5	#UVB = float(input("Enter Upper Valence Band energy in float datatype:"))

# root_dir = input()

# Set the value of coarse
coarse = 4.5
print('coarse =', coarse)

# Set the value of fine
fine = 100
print('fine =', fine)

# Storing all the files of the dir
all_files = os.listdir(root_dir)

dat_gnu_file = ''
nscf_out_file = ''
bandx_out_file = ''

for file in all_files:
    if file.endswith('.gnu'):
        dat_gnu_file = file

    if file.endswith('nscf.out'):
        nscf_out_file = file

    if file.endswith('bandx.out'):
        bandx_out_file = file

# creating csv file for bandx.dat.gnu
csv_file_handler = fa.file_analyzer()  # object to create csv file
gnu_file_path = root_dir + '/' + dat_gnu_file  # path of .gnu file
gnu_csv_file_path = gnu_file_path + '.Bandx.csv'  # path of .csv file
csv_file_handler.write_csv(gnu_file_path, gnu_csv_file_path)  # creating csv file

nscf_file_path = root_dir + "/" + nscf_out_file  # path of nscf.out file
bandx_out_path = root_dir + "/" + bandx_out_file  # path of bands.in file

# opening csv file
table = pd.read_csv(gnu_csv_file_path,
                    header=None)  # this command is used to read csv files store the read value in 'table' variable
table_values = table.values  # .values will change the vlaues of 'table' in matrix form on which all the calculation can be perfomed

r = table.shape[0]  # number of rows in file
distinct_x = csv_file_handler.get_distinct_x(table_values)  # getting no of points on x-axis

# Getting offset value to make centre at 0
min_max_values = csv_file_handler.get_offset(nscf_file_path)
offset = abs(min_max_values[0])

# Getting k_values of the plot
k_values = csv_file_handler.get_kpoints(bandx_out_path)
print('k_values =', k_values)
tot_k_values = len(k_values)

# Enter the names of points on x-axis (K_names)
K_names = ['$\Gamma$', 'M', 'K', '$\Gamma$']

all_elements_map = {}
img_names = {}
# Plotting the 2-d plot
for i in range(0, r, distinct_x):  # for each repeatition represents a new line in the plot
    X = table_values[i:i + distinct_x - 1,
        0]  # store column 0 and column 1 in x and y variable for each set and plot it,
    Y = table_values[i:i + distinct_x - 1,
        1] + offset  # 'count' value is used to get the size,adding 1.91011 to all y values to make the centre at 0
    plt.plot(X, Y, lw=0.7, antialiased=True, color='w')  # ploting the x vs y graph as amny times as many repeation

plt.margins(0)  # to remove extra spaces from boundary

# creating csv file for -ve spin As p Orbital almost same as above

all_pdos_file_paths = []  # Fetching all pdos files
for file in all_files:
    if file.endswith(')') and file.find('pdos') != -1:
        all_pdos_file_paths.append(file)

# Printing different plots for different files
for file in all_pdos_file_paths:
    # Writing csv for pdos file
    pdos_file_path = root_dir + '/' + file
    pdos_csv_file_path = pdos_file_path + '.csv'
    csv_file_handler.write_csv(pdos_file_path, pdos_csv_file_path)

    current_element = ""
    flag = 0
    for it in file:
        if it == ')':
            break
        if flag == 1:
            current_element = current_element + it
        if it == '(':
            flag = 1

    # open csv file
    table1 = pd.read_csv(pdos_csv_file_path, sep=",",
                         header=None)  # pd.read_csv return file in dataframe type. .values is used to change in array form
    table_values1 = table1.values
    distinct_y = len(table1) // distinct_x

    # select the pdos column
    c = table1.iloc[:, 2]  # iloc is used to select particular column from the csv file

    c1 = np.array(
        c)  # work same as 'c.values'.since the all data is present in column .all the y values at same x is not bundled together

    c1 = c1.reshape(distinct_x,
                    distinct_y)  # changing the linear column into matrix form,now every row represents all the y values at a paricular x value.and different row means different x value
    c1 = c1.T  # transposing the matrix to fit the plot for bandx

    x = table.iloc[:distinct_x, 0]  # all the x value from bandx plot

    min_value = table_values1[0, 1]  # minimum value of energy
    max_value = table_values1[distinct_y - 1, 1]  # maximum value of energy
    diff = round(table_values1[1, 1] - table_values1[0, 1], 4)  # difference between energy levels

    # y = np.arange(min_value + offset, max_value + offset + diff,
    #               diff)  # all the y values at which pdosvalue is given,2nd column,added 2.1998 to coincide with bandx plot

    y = []
    r1 = min_value + offset
    r2 = max_value + offset + diff
    while r1 <= r2:
        y.append(r1)
        r1 += diff

    if len(y) != distinct_y:
        if len(y) < distinct_y:
            while len(y) != distinct_y:
                y.append(y[-1] + diff)
        else:
            while len(y) != distinct_y:
                y.pop()

    d1 = y[0]
    d2 = y[-1]
    xx, yy = np.meshgrid(x, y)  # creating x ,y grid for contour plot

    extreme_sum = 0
    for i in range(len(yy)):
        for j in range(len(yy[0])):
            if LCB >= yy[i][j] >= UVB:
                extreme_sum += c1[i][j]

    breaks = np.linspace(0, coarse, fine)
    # print(current_element)
    if current_element not in all_elements_map:
        all_elements_map[current_element] = extreme_sum
    else:
        all_elements_map[current_element] = max(extreme_sum, all_elements_map[current_element])

    p1 = plt.contourf(xx, yy, c1, breaks, cmap='hot')  # creating contour plot
    plt.grid(visible=None)
    # ax = plt.axes()
    # ax.set_facecolor('#fdfcf6')

    ax1.tick_params(axis='both', colors='black', length=3, width=0.3, direction='in')  # axes tick properties

    plt.ylim(-8, 8)
    plt.xlim(0, k_values[len(k_values) - 1])  # limiting the plot to -4 to 4
    plt.grid(visible=None)

    # creating vertical and horizontal lines for conduction band and valence band
    plt.axhline(y=0, linestyle=':', lw=0.7, color='w')  # horizontal lines
    plt.axhline(y=min_max_values[1] - min_max_values[0], linestyle=':', lw=0.7, color='w')

    for i in range(0, len(k_values)):
        plt.axvline(x=k_values[i], linestyle='--', lw=0.5, color='w')

    # labeling tau M and K
    a = np.array(k_values)  # labelling particular points on x axis

    label3 = plt.xticks(a, K_names, fontsize=50)
    plt.ylabel(r'E-E$_\mathrm{F}$ (eV)', fontsize=50)

    # ax1.text(0.02, -5.7, r'(a):Ge 4$\mathrm{p}_{\mathrm{m}_{\mathrm{s}=+\frac{1}{2}}}$', fontsize=50)

    name = file + '.png'

    img_names[extreme_sum] = name
    plt.title(name, pad=50)
    fig.savefig(root_dir + '/' + name, bbox_inches='tight', format='png', dpi=300)
    print("Generated image for", file)
    # exit()

font = {'size': 5}

mpl.rc('font', **font)
# create figure
fig = plt.figure(figsize=(15, 15))

for [key, value] in all_elements_map.items():
    best_img = Image.open(root_dir + '/' + img_names[value])
    best_img.copy().save('best(' + key + ').png')
    print("Generated best image for " + key)

# setting values to rows and column variables

all_images = []
all_images_name = []
all_files = os.listdir(root_dir)
for file in all_files:
    if file.endswith('.png'):
        if file.startswith('all'):
            continue
        all_images.append(cv2.imread(root_dir + '/' + file))
        all_images_name.append(file)

# Adds a subplot at the 1st position
tot_image_size = len(all_images)

hash_map = [[0, 0], [1, 1], [1, 2], [1, 3], [2, 2], [2, 3], [2, 3], [3, 3], [3, 3], [3, 3]]

rows = 0
columns = 0
cnt = 0

if tot_image_size > len(hash_map):
    rows = tot_image_size // 4

    if tot_image_size % 4 != 0:
        rows = rows + 1

    columns = 4
else:
    rows = hash_map[tot_image_size][0]
    columns = hash_map[tot_image_size][1]

for image in all_images:
    fig.add_subplot(rows, columns, cnt + 1)
    plt.imshow(image)
    plt.axis('off')
    cnt = cnt + 1

# fig.tight_layout()
fig.savefig(root_dir + '/all.png', bbox_inches='tight', format='png', dpi=300)
print('Generated combined image. Check orbital_plot_data folder and orbital_plot_data/orbital_data folder.')
