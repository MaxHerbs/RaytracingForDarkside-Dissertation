'''
Produces a heat map for the intensity across the surface of the PDU tile 
Data is produced by pdu_sim_master, stored in output_data/pdu_intensities.csv
Stores heatmap for spacial distribution of intensities in figures/pdu_heatmap.png

'''



from matplotlib import gridspec
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cbook, cm
from matplotlib.colors import LightSource
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import pandas as pd
import sys
from math import sqrt
from scipy.optimize import curve_fit
from scipy.stats import norm
import math 
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.optimize import curve_fit
import sys


pdu_file_path = "objects/10x10new.obj"

intensities_file_path = "output_data/pdu_intensities.csv"
PI = 3.14159265359

vertex_coordinate = np.zeros(shape=(59,3))
face_centre_coord = np.zeros(shape=(108,3))
object_file = open(pdu_file_path)
object_file.readline()
object_file.readline()
object_file.readline()
object_file.readline()

vertex_index = 1
face_index = 0


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def load_obj(file_path):
    vertices = []
    faces = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('v '):
                vertex = list(map(float, line.strip().split()[1:]))
                vertices.append(vertex)
            elif line.startswith('f '):
                face = list(map(int, line.strip().split()[1:]))
                face = [index - 1 for index in face]
                faces.append(face)
    return np.array(vertices), np.array(faces)


# Load .obj file
vertices, faces = load_obj('objects/10x10new.obj')

# Plot the object
x,y,z = [],[],[]
for face in faces:
        face_vertices = vertices[face, :]
        face_center = np.mean(face_vertices, axis=0)
        x.append(face_center[0])
        y.append(face_center[1])
        z.append(face_center[2])

z = np.array(z)
        

intensities = []
err_intensities = []
intensities_file = open(intensities_file_path, "r")
intensities_file.readline()

for line in intensities_file.readlines():
    data = line.replace("\n","").split(",")
    intensities.append(float(data[1]))
    err_intensities.append(float(data[2]))

intensities_array = np.array(intensities)
#intensities_array[np.where(intensities_array < 0)] = 0
err_intens_array = np.array(err_intensities)






bad_coords = np.where(z != 0)
x = np.delete(x, bad_coords)
y = np.delete(y,bad_coords)
z = np.delete(z, bad_coords)
intensities = np.delete(intensities_array, bad_coords)
err_intensities = np.delete(err_intens_array, bad_coords)
perc_err = np.abs(err_intensities/intensities * 100)







def get_index(x, xmin, xmax):
    coord_range = xmax - xmin
    relative_position = (x - xmin) / coord_range
    if relative_position < 0.1:
        return 9
    elif 0.1 < relative_position <= 0.2:
        return 8 
    elif 0.2 < relative_position <= 0.3:
        return 7
    elif 0.3 < relative_position <= 0.4:
        return 6
    elif 0.4 < relative_position <= 0.5:
        return 5
    elif 0.5 < relative_position <= 0.6:
        return 4
    elif 0.6 < relative_position <= 0.7:
        return 3
    elif 0.7 < relative_position <= 0.8:
        return 2
    elif 0.8 < relative_position <= 0.9:
        return 1
    elif relative_position > 0.9:
        return 0
    else:
        print("Problematic value: {}\nRelative Position: {}".format(x, relative_position))
        return "err"


def sort_into_grid(intensity_values, x_coordinates, y_coordinates, normalise_face_num=True):
    # Create a 4x4 grid to store the sum of intensities for each cell
    grid = np.zeros((10, 10), dtype=float)
    # Create a grid to store the number of points contributing to each cell
    count_grid = np.zeros((10, 10), dtype=int)
    print(np.mean(y_coordinates))
    # Normalize x and y coordinates to fit within the 4x4 grid
    max_x = 0.1
    max_y = 0.1
    min_x = -0.1
    min_y = -0.1

    for i in range(len(x_coordinates)):
        x_index = get_index(x_coordinates[i], min_x, max_x)
        y_index = get_index(y_coordinates[i], min_y, max_y)
        grid[y_index, x_index] += intensity_values[i]
        count_grid[y_index, x_index] += 1
    print(count_grid)
    # Normalize the grid by dividing each cell by the number of contributing points
    if(normalise_face_num):
        for y in range(10):
            for x in range(10):
                if count_grid[y, x] != 0:
                    grid[y, x] /= count_grid[y, x]

    return grid








#Upper Grid
upper_index = np.where(y < 0)
upper_x = np.delete(x, upper_index)
upper_y = np.delete(y, upper_index)
upper_intensities = np.delete(intensities, upper_index)
upper_intensities_error = np.delete(perc_err, upper_index)

upper_grid = sort_into_grid(upper_intensities, upper_x, upper_y-0.11)
upper_error_grid = sort_into_grid(upper_intensities_error, upper_x, upper_y-0.11)

#Lower Grid
lower_index = np.where(y > 0)
lower_x = np.delete(x, lower_index)
lower_y = np.delete(y, lower_index)
lower_intensities = np.delete(intensities, lower_index)
lower_intensities_error = np.delete(perc_err, lower_index)

lower_grid = sort_into_grid(lower_intensities, lower_x, lower_y+0.11)
lower_error_grid = sort_into_grid(lower_intensities_error, lower_x, lower_y+0.11)

plt.rcParams.update({'font.size': 16})  

whole_grid = np.vstack((upper_grid, lower_grid))
whole_grid_error = np.vstack((upper_error_grid, lower_error_grid))
whole_grid = whole_grid/np.max(whole_grid)

bar = plt.imshow(whole_grid, extent=[-10,10,-20,20])
plt.xlabel("Horizontal Position on PDUs/cm")
plt.ylabel("Vertical Position on PDUs/cm")
plt.colorbar(bar, label="Relative Intensity")
plt.savefig("figures/10x10_variation.png", bbox_inches='tight')
plt.clf()

bar = plt.imshow(whole_grid_error,extent=[-10,10,-20,20])
plt.colorbar(bar, label="Percentage Error of Intensity")
plt.xlabel("Horizontal Position on PDUs/cm")
plt.ylabel("Vertical Position on PDUs/cm")

plt.savefig("figures/10x10_variation_err.png", bbox_inches='tight')

plt.clf()






flat_array = np.ndarray.flatten(whole_grid)
flat_error = np.ndarray.flatten(whole_grid_error)


max_val = np.max(flat_array)
min_val = np.min(flat_array)

max_index = np.where(flat_array == max_val)
min_index = np.where(flat_array == min_val)

power_range = (max_val - min_val) * 100
power_range_error = flat_error[max_index] + flat_error[min_index]
print("Total Variation = {:0.2f}% +/- {:0.2f}%".format(power_range, power_range_error[0]))













cum_intensity = []
cum_intensity_error = []
for i in range(10):
    cum_total = 0
    cum_total_error = 0
    for j in range(10):
        cum_total += whole_grid[i,j]
        cum_total_error += whole_grid_error[i,j]
    cum_intensity.append(cum_total/10)
    cum_intensity_error.append(cum_total_error/1000)

cum_x_axis = np.linspace(-0.175,0.175,10)
cum_x_axis = np.flip(cum_x_axis)
granular_x = np.linspace(-0.2, 0.2, 200)


def double_norm(x,a,sig,mu,c, b,sig2,mu2):
    return a*np.exp(-0.5*((x-mu)/sig)**2) + b*np.exp(-0.5*((x-mu2)/sig2)**2) + c


bad_coords = (0,1)
cum_intensity = np.delete(cum_intensity, bad_coords)
cum_x_axis = np.delete(cum_x_axis, bad_coords)
cum_intensity_error = np.delete(cum_intensity_error, bad_coords)

plt.errorbar(cum_x_axis, cum_intensity, cum_intensity_error, fmt="o", capsize=5, markersize=5)
try:
    dn_parmas, dn_cov = curve_fit(double_norm, cum_x_axis, cum_intensity)
    plt.plot(granular_x, double_norm(granular_x, *dn_parmas))
    print("Fit one:\namp={}\nsig={}\nmu={}\nFit two:\namp={}\nsig={}\nmu={}\n+c={}".format(dn_parmas[0],dn_parmas[1],dn_parmas[2],dn_parmas[4],dn_parmas[5],dn_parmas[6], dn_parmas[3]))

except:
    print("Variation with verticle position fit failed")


plt.savefig("figures/verticle_intensity_change.png")
plt.clf()












def norm_fit(x, a, sig, mu):
    return a * np.exp(-0.5 * ((x-mu)/sig)**2)
grid_intensities = np.ndarray.flatten(whole_grid)
plt.rcParams.update({'font.size': 14})  
plt.xlabel("Relative Intensity")
plt.ylabel("Frequency")

bin_num = 12
plt.hist(grid_intensities, bins=bin_num)
hist, bins = np.histogram(grid_intensities, bin_num)


try:
    bins = bins + 0.5*(bins[1]-bins[0])
    p0 = [30,0.02,0.9]
    normal_params, cov = curve_fit(norm_fit, bins[:-1], hist, p0=p0)
    hist_x = np.linspace(np.min(bins[:-1]), np.max(bins[:-1]), 200)
    plt.plot(hist_x, norm_fit(hist_x, *normal_params), "r")
    txt = "amp = {:0.2f} +/- {:0.2f}\nsigma = {:0.3f} +/- {:0.3f}\nmu = {:0.2f} +/- {:0.2f}".format(normal_params[0], sqrt(cov[0][0]),normal_params[1],sqrt(cov[1][1]),normal_params[2],sqrt(cov[2][2]) )
    print(txt)
except:
    print("Failed to fit hisogram")

plt.savefig("figures/10x10_histogram.png")
plt.clf()


