'''
Master simulation of PDU geometry and light source for real life model.

Contains:
Cryostat    Two light sources   16 PDU tiles

Outputs intensity data for faces of the simulated PDU

'''




# External imports
import csv
import os
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
import time
import datetime

# Internal imports
from raysect.optical import World, translate, rotate, Point3D, d65_white, ConstantSF, Node
from raysect.optical.observer import PinholeCamera, RGBPipeline2D, RGBAdaptiveSampler2D, MeshPixel, MeshCamera, PowerPipeline0D, PowerPipeline1D, MonoAdaptiveSampler1D
from raysect.optical.material.emitter import UniformVolumeEmitter, UniformSurfaceEmitter, UnityVolumeEmitter, InhomogeneousVolumeEmitter
from raysect.optical.material import Lambert, AbsorbingSurface, Roughen
from raysect.primitive import Box, Subtract, export_vtk
from raysect.primitive.mesh import import_obj
from raysect.optical.library import schott, Gold, Aluminium, Silicon
from raysect.primitive import Sphere, Box
from math import sin,cos
from raysect.optical.library.spectra.colours import *
import numpy as np


    
PI = 3.14159265359
cryostat_radius = 0.25
num_per_layer = 8
light_source_radius = 25 * 10**(-3)/2
pdu_model = "objects/single_slab.obj"
large_pdu_model = "objects/double_slab.obj"
bin_num = 500




light_x = np.linspace(350,750,500)
light_y = np.zeros(len(light_x))
light_y[np.where(np.round(light_x,0)==405.0)] = 100

light_x = np.linspace(350,750,500)
light_y = np.zeros(len(light_x))
light_y[np.where(np.round(light_x,0)==405.0)] = 100
class Cos2Glow(InhomogeneousVolumeEmitter): 
    def emission_function(self, point, direction, spectrum, world, ray, primitive, to_local, to_world): 
        w = 0.5704976547575429 
        b = -1.586601617597982e-12 
        c = -0.14338489427162548 
        amp = 1.0880679340886772 
        ang = np.arcsin(abs(point.y)/light_source_radius) 
        ang = PI/2 - ang 
        spectrum.samples[:] = (amp * np.cos(1.925*w*ang + b)**2 + c) * light_y
        return spectrum 




def write_results(basename, camera, mesh):
    # obtain frame from pipeline
    frame = camera.pipelines[0].frame

    # calculate power density
    np.savetxt(basename + "face_areas.txt", camera.collection_areas)

    power_density = frame.mean / camera.collection_areas
    error = frame.errors() / camera.collection_areas

    # write as csv
    with open(basename + 'pdu_intensities.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(('Triangle Index', 'Power Density / Wm-2', 'Error / Wm-2'))
        for index in range(frame.length):
            writer.writerow((index, power_density[index], error[index]))

    triangle_data = {'PowerDensity': power_density, 'PowerDensityError': error}
    export_vtk(mesh, basename + 'paraviewfile.vtk', triangle_data=triangle_data)









world = World()

sphere = Sphere(light_source_radius)
box = Box(Point3D(-light_source_radius,-light_source_radius,-light_source_radius), Point3D(light_source_radius,0,light_source_radius))
bottom_light = Subtract(sphere, box, parent=world, material=Cos2Glow(), transform=translate(0,-0.22,0))
top_light = Subtract(sphere, box, parent=world, material=Cos2Glow(), transform=translate(0,0.22,0)*rotate(0,180,0))

#floor = Box(Point3D(-100,-0.5,-100), Point3D(100,-0.4,100), parent=world, material=Lambert())

#region Build PDU Layout
lower_pdu_array = []
upper_pdu_array = []



large_pdu = import_obj(large_pdu_model, parent=world, transform=translate(cryostat_radius,0,0)*rotate(90,180,0), material=Roughen(Silicon(), 0.125))

for i in range(1,8):
    angle = 360*i/num_per_layer
    angle_rad = angle*(2*PI)/360
    x = cryostat_radius * cos(angle_rad)
    z = cryostat_radius * sin(angle_rad)

    pduU = import_obj(pdu_model, parent=world, transform=translate(x,0.11,z)*rotate(angle + 90,180,0), material=Roughen(Silicon(), 0.125))
    pduL = import_obj(pdu_model, parent=world, transform=translate(x,-0.11,z)*rotate(angle + 90,180,0), material=Roughen(Silicon(), 0.125))


    lower_pdu_array.append(pduL)
    upper_pdu_array.append(pduU)


#endregion

#region Add Cryostat Chamber
cryostat = import_obj("objects/cryostat.obj", parent=world, transform=translate(0,0,0)*rotate(0,0,0), material=Roughen(Aluminium(), 0.21875))
lid = Box(Point3D(-10,0.5,-10), Point3D(10,0.55,10),parent=world, material=Roughen(Aluminium(), 0.21875))
#center_light = Sphere(2*light_source_radius, parent=world, material=UniformSurfaceEmitter(d65_white), transform=translate(0,0,0))

#endregion



#region Mesh Observer


min_wl = 350
max_wl = 700
samples = 1000000

power = PowerPipeline1D()
sampler = MonoAdaptiveSampler1D(power, fraction=1, ratio=5, min_samples=1000, cutoff=0.01)
camera = MeshCamera(
    large_pdu,
    surface_offset=1e-6,  # launch rays 1mm off surface to avoid intersection with absorbing mesh
    pipelines=[power],
    frame_sampler=sampler,
    parent=world,
    spectral_bins=bin_num,
    min_wavelength=min_wl,
    max_wavelength=max_wl,
    pixel_samples=1000,
    transform=translate(cryostat_radius,0,0)*rotate(90,180,0),
    quiet=True
)

print('Observing the PDU mesh...')
output_basename = "output_data/"
render_pass = 0
while (not camera.render_complete):
    render_pass += 1
    print('Render pass {}:'.format(render_pass))
    camera.observe()
    write_results(output_basename, camera, large_pdu)

print('Observation complete!')
write_results(output_basename, camera, large_pdu)
# export final data as csv



#endregion








# # CAMERA
# rgb = RGBPipeline2D(display_unsaturated_fraction=0.96, name="sRGB")
# sampler = RGBAdaptiveSampler2D(rgb, ratio=10, fraction=0.2, min_samples=2000, cutoff=0.01)
# camera = PinholeCamera((128, 128), parent=world, transform=translate(0, 1, 0) * rotate(0,-90,0), pipelines=[rgb], frame_sampler=sampler, fov=45)
# camera.spectral_rays = 1
# camera.spectral_bins = 500
# camera.pixel_samples = 250
# camera.ray_max_depth = 500
# camera.ray_extinction_min_depth = 3
# camera.ray_extinction_prob = 0.01

# # RAY TRACE


# start_time = datetime.datetime.now()
# ion()
# name = 'single_pdu_render'
# timestamp = time.strftime("%Y-%m-%d_%H-%M-%S")
# render_pass = 1
# max_renders = 4
# while (not camera.render_complete) and render_pass < max_renders:
#     print("Rendering pass {}...".format(render_pass))
#     camera.observe()
#     rgb.save("{}_toppass_{}.png".format(name, render_pass))
#     print()
#     render_pass += 1

# ioff()
# rgb.display()

# print("Runtime: {}".format(datetime.datetime.now() - start_time))
