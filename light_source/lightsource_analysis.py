'''
Produces scatter plot of the geometry of the imported arc.obj
Then produces a scatter of the angular intensity shown against a plot of the real life data

'''




import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cbook, cm
from matplotlib.colors import LightSource
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import pandas as pd
PI = 3.14159265359
import matplotlib.gridspec as gridspec

plt.rcParams.update({'font.size': 22})

pdu_file_path = "objects/arc.obj"
intensities_file_path = "output_data/lightsource_arc_intensities.csv"
fig_path = "figures/lightsource_"
data_path = "output_data/"

vertex_coordinate = np.zeros(shape=(97,3))
face_centre_coord = np.zeros(shape=(188,3))
object_file = open(pdu_file_path)
object_file.readline()
object_file.readline()
object_file.readline()
object_file.readline()

vertex_index = 1
face_index = 0


def get_position(row):
    return vertex_coordinate[row]

def get_point(x, row):
    if x == "x":
        return vertex_coordinate[row][0]
    if x == "y":
        return vertex_coordinate[row][1]
    if x == "z":
        return vertex_coordinate[row][2]





object_data = object_file.readlines()
object_file.close()
for line in object_data:
    data = line.replace("\n", "").split(" ")
    if data[0] == "v":
        vertex_coordinate[vertex_index][0] = float(data[1])
        vertex_coordinate[vertex_index][1] = float(data[2])
        vertex_coordinate[vertex_index][2] = float(data[3])
        vertex_index += 1
    if data[0] == "f":
        face_centre_coord[face_index][0] = get_point("x", int(data[1])) + get_point("x", int(data[2])) + get_point("x", int(data[3]))
        face_centre_coord[face_index][1] = get_point("y", int(data[1])) + get_point("y", int(data[2])) + get_point("y", int(data[3]))
        face_centre_coord[face_index][2] = get_point("z", int(data[1])) + get_point("z", int(data[2])) + get_point("z", int(data[3]))
        face_index += 1
face_centre_coord = face_centre_coord/3





intensities = []
intensities_err = []
intensities_file = open(intensities_file_path, "r")
intensities_file.readline()

for line in intensities_file.readlines():
    data = line.replace("\n","").split(",")
    intensities.append(float(data[1]))
    intensities_err.append(float(data[2]))

intensities_array = np.array(intensities)
intensities_err_array = np.array(intensities_err)

all_radii = []
radii = []
for i in range(len(intensities_array)):
    all_radii.append(round(np.sqrt(face_centre_coord[i][0]**2 + face_centre_coord[i][1]**2 + face_centre_coord[i][2]**2),3))
    if round(np.sqrt(face_centre_coord[i][0]**2 + face_centre_coord[i][1]**2 + face_centre_coord[i][2]**2),3) != 0.075:
        radii.append(i)
trans = face_centre_coord.transpose()


fig = plt.figure(figsize = (10,10))
ax = plt.axes(projection='3d')
ax.grid()
ax.scatter(trans[0], trans[2], trans[1], c = 'r', s = 50)
ax.set_title('Geomety of Arc Faces')
# Set axes label
ax.set_xlabel('x', labelpad=20)
ax.set_ylabel('z', labelpad=20)
ax.set_zlabel('y', labelpad=20)

plt.savefig(fig_path + "geometry.png")
plt.clf()


x = np.delete(trans[0], radii)
y = np.delete(trans[1], radii)
z = np.delete(trans[2], radii)
intensities_array = np.delete(intensities_array, radii)
intensities_err_array = np.delete(intensities_err_array, radii)


norm_x = x/np.linalg.norm(x)
norm_y = y/np.linalg.norm(y)
norm_z = z/np.linalg.norm(z)
norm_intensities = intensities_array/np.linalg.norm(intensities_array)
norm_errors = intensities_err_array/np.linalg.norm(intensities_array)

plt.plot(norm_x,norm_y)
plt.scatter(norm_x,norm_intensities, marker="+")
plt.xlabel("HorizontalPosition")
plt.ylabel("Verticle Position")
plt.title("Flattened Arc Geometry of Inside Face")
plt.savefig(fig_path + "2d_geometry.png")
plt.clf()

def cos2(x, a, b, c, amp):
    return amp*(np.cos(a*x + b)**2) + c

def fn(x):
    return 0.173557784474499*(np.cos(0.6599849191038133*x)**2)

w = 0.5704976547575429
b = -1.586601617597982e-12
c = -0.14338489427162548
amp = 1.0880679340886772


angles = np.zeros(len(intensities_array))
angles = np.arcsin(x/np.max(x))
#angles = np.arcsin()
print(np.sort(angles))


theory_y = fn(np.sort(angles))
theory_y = theory_y/np.linalg.norm(theory_y)

extreme_angles = np.where(np.abs(angles) > np.deg2rad(81))
print(extreme_angles)
angles = np.delete(angles,extreme_angles)
theory_y = np.delete(theory_y,extreme_angles)
norm_intensities = np.delete(norm_intensities, extreme_angles)
norm_errors = np.delete(norm_errors, extreme_angles)

gs = gridspec.GridSpec(2,1, height_ratios=[4, 1])
fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(gs[0])
#ax.set_title("Normalised Intensity vs Parameter Plot")
ax.plot(np.sort(angles), theory_y, label="Input fit")
#plt.scatter(angles,norm_intensities, marker="+",)

ax.errorbar(angles,norm_intensities, yerr=norm_errors, fmt="o", label="Simulation Output")
ax.legend()
ax.set_ylabel("Normalised Intensity/W")


f = open("output_data/lightsource_simulation__by_angles.txt", "w")
for i in range(len(angles)):
    f.write("{} {}\n".format(angles[i], norm_intensities[i]/np.max(norm_intensities)))
f.close()

ax = fig.add_subplot(gs[1], sharex=ax)
resid = (theory_y-norm_intensities)/norm_errors
ax.scatter(angles,resid)
ax.set_ylim((-3.5,3.5))
ax.set_ylabel("Normalised \nResiduals")
ax.axhspan(-1, 1, color='gray', alpha=0.3)
proportion_in_1 = len(np.where(resid <=1)[0])/len(resid)

proportion_in_2 = len(np.where(resid <=2)[0])/len(resid)
proportion_in_3 = len(np.where(resid <=3)[0])/len(resid)
print("Percentage within 1σ: {}% \n Percentage within 2σ: {}%\n Percentage within 3σ: {}%".format(proportion_in_1*100,proportion_in_2*100,proportion_in_3*100))
plt.legend()
plt.xlabel("Angular position on arc/radians")


#plt.xlim((-1.4,1.4))
plt.savefig(fig_path + "angular_intensity.png", bbox_inches='tight') 
plt.clf()



# #Side by side

# filename = "parameters/intensity_25mm.txt"
# f = open(filename)

# a = []
# b = []
# c = []
# for line in f.readlines():
#     data = line.replace("\n","").split()
#     a.append(float(data[0]))
#     b.append(float(data[1]))
#     c.append(float(data[2]))


# theta = np.array(a)
# phi = np.array(b)
# intensity = np.array(c)
# sparse_points = np.where(phi != 0.0)

# theta = np.delete(theta, sparse_points)
# inten = np.delete(intensity, sparse_points)
# theta = np.deg2rad(theta)
# inten = inten/np.linalg.norm(inten)

# theory_y = theory_y/np.linalg.norm(theory_y)
# inten = inten/np.linalg.norm(inten)
# norm_intensities = norm_intensities/np.linalg.norm(norm_intensities)



# plt.plot(np.sort(angles), theory_y)
# plt.scatter(theta, inten, marker="+")
# plt.scatter(angles,norm_intensities, marker="o")

# plt.savefig("figures/lightsource_sidebyside.png")