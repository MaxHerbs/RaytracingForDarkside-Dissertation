'''
Analysis of the lightsource we have designed.
Surrounds the hemisphere with an arc and records the incident intensities on every face which
can be analysed.
Also analysises the light spectrum of the diffuser to confirm continuity.

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
from raysect.optical.material.emitter import UniformVolumeEmitter, UniformSurfaceEmitter, UnityVolumeEmitter
from raysect.optical.material import Lambert, AbsorbingSurface
from raysect.primitive import Box, Subtract, export_vtk, Intersect
from raysect.primitive.mesh import import_obj
from raysect.optical.library import schott, Gold, Aluminium
from raysect.primitive import Sphere, Box, Union
from math import sin,cos
from raysect.optical.library.spectra.colours import *
import os
import time
import numpy as np
from matplotlib.pyplot import *
from numpy import sqrt, cos, pi, arcsin

from raysect.optical import World, translate, rotate, Point3D
from raysect.optical.library import RoughTitanium
from raysect.optical.material import InhomogeneousVolumeEmitter, AbsorbingSurface, Checkerboard, UniformSurfaceEmitter
from raysect.optical.library.spectra.colours import orange
from raysect.optical.observer import PinholeCamera, RGBPipeline2D, RGBAdaptiveSampler2D, Pixel, PowerPipeline0D
from raysect.primitive import Box, Sphere, Subtract, Cylinder
from raysect.core.workflow import MulticoreEngine

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
from raysect.optical.material import Lambert, AbsorbingSurface
from raysect.primitive import Box, Subtract, export_vtk
from raysect.primitive.mesh import import_obj
from raysect.optical.library import schott, Gold, Aluminium, Silicon
from raysect.primitive import Sphere, Box
from math import sin,cos
from raysect.optical.library.spectra.colours import *
from raysect.optical.observer import FibreOptic, PowerPipeline0D, RadiancePipeline0D, SpectralPowerPipeline0D, SpectralRadiancePipeline0D
#endregion


    
PI = 3.14159265359
cryostat_radius = 0.25
num_per_layer = 8
light_source_radius = 25 * 10**(-3)/2
arc = "objects/arc.obj"
output_basename = "output_data/lightsource_"
fig_basename = "figures/lightsource_"
bin_num = 500


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
        #spectrum.samples[:] = (np.cos(1.66*0.6599849191038133*ang)**2)
        return spectrum 







def write_results(basename, camera, mesh):
    # obtain frame from pipeline
    frame = camera.pipelines[0].frame

    power_density = frame.mean / camera.collection_areas
    error = frame.errors() / camera.collection_areas

    # write as csv
    with open(basename + 'arc_intensities.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(('Triangle Index', 'Power Density / Wm-2', 'Error / Wm-2'))
        for index in range(frame.length):
            writer.writerow((index, power_density[index], error[index]))

    triangle_data = {'PowerDensity': power_density, 'PowerDensityError': error}
    export_vtk(mesh, basename + 'out.vtk', triangle_data=triangle_data)









world = World()
sphere = Sphere(light_source_radius)
box = Box(Point3D(-light_source_radius,-light_source_radius,-light_source_radius), Point3D(light_source_radius,0,light_source_radius))
#dummy_light = Sphere(light_source_radius, parent=world, material=Cos2Glow(), transform=translate(0,0,0))
light = Subtract(sphere, box, parent=world, material=Cos2Glow(), transform=translate(0,0,0))
#shadow = Sphere(light_source_radius/10, parent=world, material=AbsorbingSurface(), transform=translate(cryostat_radius/2,-0.11,0))


#floor = Box(Point3D(-100, -0.1, -100), Point3D(100, -0.01, 100), world,transform=translate(0,-0.05,0), material=Lambert(ConstantSF(1.0)))

mesh = import_obj(arc, parent=world, transform=translate(0,0,0)*rotate(0,0,0), material=AbsorbingSurface())



#region Mesh Observer


min_wl = 350
max_wl = 700

power = PowerPipeline1D()
sampler = MonoAdaptiveSampler1D(power, fraction=0.2, ratio=25.0, min_samples=10000, cutoff=0.01)
camera = MeshCamera(
    mesh,
    surface_offset=1e-6,  # launch rays 1mm off surface to avoid intersection with absorbing mesh
    pipelines=[power],
    frame_sampler=sampler,
    parent=world,
    spectral_bins=500,
    min_wavelength=min_wl,
    max_wavelength=max_wl,
    pixel_samples=100000,
    transform=translate(0,0,0)*rotate(0,0,0)
)
print('Observing the arc mesh...')
render_pass = 0
while (not camera.render_complete) and (render_pass < 10):
    render_pass += 1
    print('Render pass {}:'.format(render_pass))
    camera.observe()
    write_results(output_basename, camera, mesh)

print('Observation complete!')
write_results(output_basename, camera, mesh)
# export final data as csv


#endregion


# #region Spectral Analyis
# plt.ion()
# spectral_power = SpectralPowerPipeline0D()
# spectral_radiance = SpectralRadiancePipeline0D()
# power = PowerPipeline0D()
# radiance = RadiancePipeline0D()
# fibre = FibreOptic([spectral_power, spectral_radiance, power, radiance], acceptance_angle=10, radius=1,
#                    spectral_bins=bin_num, spectral_rays=1, pixel_samples=1000000, transform=translate(0, 1, 0)*rotate(0,180,0), parent=world)
# fibre.observe()

# plt.ioff()
# plt.xlim((400,450))
# plt.savefig(fig_basename + "lightsource_spectrum.png")
# plt.clf()

# #endregion


# # region Camera
# # CAMERA
# rgb = RGBPipeline2D(display_unsaturated_fraction=0.96, name="sRGB")
# sampler = RGBAdaptiveSampler2D(rgb, ratio=10, fraction=0.2, min_samples=2000, cutoff=0.01)
# camera = PinholeCamera((512, 512), parent=world, transform=translate(-0.065,light_source_radius,-0.065) * rotate(-45,-15,0), pipelines=[rgb], frame_sampler=sampler, fov=45)
# camera.spectral_rays = 1
# camera.spectral_bins = bin_num
# camera.pixel_samples = 250
# camera.ray_max_depth = 500
# camera.ray_extinction_min_depth = 3
# camera.ray_extinction_prob = 0.01


# start_time = datetime.datetime.now()
# ion()
# timestamp = time.strftime("%Y-%m-%d_%H-%M-%S")
# render_pass = 1
# max_renders = 20
# while (not camera.render_complete) and render_pass < max_renders:
#     print("Rendering pass {}...".format(render_pass))
#     camera.observe()
#     rgb.save("{}Light_render_{}.png".format(fig_basename,render_pass))
#     print()
#     render_pass += 1

# ioff()
# rgb.display()
# #endregion
