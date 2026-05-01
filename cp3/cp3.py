# list of imports
import openmc
import os
from IPython.display import Image
import matplotlib.pyplot as plt

LAST_CROSS_SECTIONS = "/home/oss/ARFC/openmc_installation/XS_data_raw/jeff-3.3-hdf5/cross_sections.xml" # last place we installed jeff data
# os.environ["OPENMC_CROSS_SECTIONS"] = LAST_CROSS_SECTIONS # in case the environment variable ever gets unset


my_cell = openmc.Cell(name="box")  # creation of a cell
model = openmc.Model()  # creation of a model

zirconium = openmc.Material(1, "zirconium")  # create material (modify from XSection)
mat = openmc.Material()
zirconium.add_element("Zr", 1.0)
zirconium.set_density("g/cm3", 6.5)
# print(zirconium)

uo2 = openmc.Material(name="uo2")  # create UO2 material
uo2.add_element("U", 1.0, enrichment=3.5)
uo2.add_nuclide("O16", 2.0)
uo2.set_density("g/cm3", 11.0)
# print(uo2)

water = openmc.Material(name="water")  # create coolant material
water.add_element("H", 2.0)
water.add_nuclide("O16", 1.0)
water.set_density("g/cm3", 1.0)

water.add_s_alpha_beta("c_H_in_H2O")
# print(water)

model.materials = openmc.Materials(
    [uo2, zirconium, water]
)  # combine materials into the model
print(model.materials)

# create the geometry for the pin cell
fuel_outer_radius = openmc.ZCylinder(r=0.46955)
clad_inner_radius = openmc.ZCylinder(r=0.47910)
clad_outer_radius = openmc.ZCylinder(r=0.54640)

fuel_region = -fuel_outer_radius
gap_region = +fuel_outer_radius & -clad_inner_radius
clad_region = +clad_inner_radius & -clad_outer_radius

height = 300.0
top = openmc.ZPlane(z0=height/2, boundary_type="vacuum")
bot = openmc.ZPlane(z0=-height/2, boundary_type="vacuum")
layer = +bot & -top



fuel = openmc.Cell()
fuel.fill = uo2
fuel.region = fuel_region & layer

gap = openmc.Cell()
gap.region = gap_region & layer

clad = openmc.Cell()
clad.fill = zirconium
clad.region = clad_region & layer

pitch = 1.44270
left = openmc.XPlane(-pitch / 2, boundary_type="reflective")
right = openmc.XPlane(pitch / 2, boundary_type="reflective")
front = openmc.YPlane(-pitch / 2, boundary_type="reflective")
back = openmc.YPlane(pitch / 2, boundary_type="reflective")

water_region = +left & -right & +front & -back & +clad_outer_radius & layer

moderator = openmc.Cell()
moderator.fill = water
moderator.region = water_region

root_universe = openmc.Universe(cells=(fuel, gap, clad, moderator))
model.geometry = openmc.Geometry(root_universe)
# root_universe.plot(width=(pitch, pitch), basis='xy')


# Create a point source
point = openmc.stats.Point((0, 0, 0))
source = openmc.IndependentSource(space=point)

model.settings.source = source
model.settings.batches = 400
model.settings.inactive = 200
model.settings.particles = 1000

model.export_to_xml()
# plot = openmc.Plot.from_geometry(model.geometry)
# plot.pixels = (250, 250)
# plot.to_ipython_image()

#### 2.a)

# pcounts = [1000, 2000, 4000]
# statepoints = []
# kvals = []
# p2a_path = lambda pcount: f"p2a/statepoint_{pcount}.h5"
# for pcount in pcounts:
#     model.settings.particles = pcount
#     newpath = p2a_path(pcount)
#     if not os.path.exists(newpath):
#         statepoint = model.run()
#         os.rename(statepoint, newpath)
#     # model.settings.export_to_xml()

    
#     with openmc.StatePoint(newpath) as sp:
#         statepoints.append(sp)
#         kvals.append(sp.keff)
    
# import matplotlib.pyplot as plt

# knoms = [float(k.nominal_value) for k in kvals]
# kstds = [float(k.std_dev) for k in kvals]
# print(f"knoms: {knoms}")
# print(f"kdevs: {kstds}")
# # plt.plot(pcounts, knoms)
# # plt.plot([1, 2, 3], [1, 3, 2], label="test data")
# plt.errorbar(pcounts, knoms, yerr=kstds, label="400 batches (200 inactive)")
# plt.xlabel("Number of particles")
# plt.ylabel("k_eff")
# plt.legend()
# plt.show()


#### 2.b)

# model.settings.particles = 2000
# mesh = openmc.RegularMesh()
# mesh.lower_left = (-pitch/2, -pitch/2)
# mesh.upper_right = (pitch/2, pitch/2)
# mesh.dimension = (50, 50)
# mesh_filter = openmc.MeshFilter(mesh)
# heat = openmc.Tally()
# heat.scores = ['kappa-fission']
# heat.filters = [mesh_filter]
# model.tallies += [heat]

# p2b_path = "p2b/heat_tallied.h5"
# if not os.path.exists(p2b_path):
#     statepoint = model.run()
#     os.rename(statepoint, p2b_path)
# with openmc.StatePoint(p2b_path) as sp:
#     mesh_tally_out = sp.get_tally(id=heat.id)
# tally = mesh_tally_out.get_values().reshape(mesh.dimension)
# # print(tally.)
# plt.figure()
# hp = pitch/2
# img = plt.imshow(tally/1000, extent=[-hp, hp, -hp, hp])
# plt.colorbar(img).set_label("Recoverable Fission Energy (keV / source particle)")
# plt.xlabel("x (cm)")
# plt.ylabel("y (cm)")
# plt.show()


#### 2.c) 

# pitches = [1.10, 1.44720, 1.80, 2.00]
# kvals = []
# for pitch in pitches:
#     left = openmc.XPlane(-pitch / 2, boundary_type="reflective")
#     right = openmc.XPlane(pitch / 2, boundary_type="reflective")
#     front = openmc.YPlane(-pitch / 2, boundary_type="reflective")
#     back = openmc.YPlane(pitch / 2, boundary_type="reflective")

#     water_region = +left & -right & +front & -back & +clad_outer_radius & layer

#     moderator = openmc.Cell()
#     moderator.fill = water
#     moderator.region = water_region

#     root_universe = openmc.Universe(cells=(fuel, gap, clad, moderator))
#     model.geometry = openmc.Geometry(root_universe)
#     # root_universe.plot(width=(pitch, pitch), basis='xy')


#     # Create a point source
#     point = openmc.stats.Point((0, 0, 0))
#     source = openmc.IndependentSource(space=point)

#     model.settings.source = source
#     model.settings.batches = 400
#     model.settings.inactive = 200
#     model.settings.particles = 1000

#     model.export_to_xml()

#     p2c_path = f"p2c/sp_pitch_{str(pitch).replace(".", "p")}.h5"
#     if not os.path.exists(p2c_path):
#         statepoint = model.run()
#         os.rename(statepoint, p2c_path)
#     with openmc.StatePoint(p2c_path) as sp:
#         kvals.append(sp.keff)

# knoms = [float(k.nominal_value) for k in kvals]
# kstds = [float(k.std_dev) for k in kvals]
# print(f"knoms: {knoms}")
# print(f"kdevs: {kstds}")
# # plt.plot(pcounts, knoms)
# # plt.plot([1, 2, 3], [1, 3, 2], label="test data")
# plt.errorbar(pitches, knoms, yerr=kstds, label="400 batches (200 inactive), p=1000")
# plt.xlabel("pitch (cm)")
# plt.ylabel("k_eff")
# plt.legend()
# plt.show()


#### 2.d) 

water = openmc.Material(name="water")  # create coolant material
water.add_element("H", 2.0)
water.add_nuclide("O16", 1.0)
water.set_density("g/cm3", 1.0)


glucose = openmc.Material(name="glucose")
glucose.add_element("C", 6)
glucose.add_element("O", 6)
glucose.add_element("H", 12)
corn_syrup = openmc.Material.mix_materials([water, glucose], [0.2, 0.8], "wo")
corn_syrup.name = "corn syrup"
corn_syrup_densities = [0.5, 1, 1.33, 1.38, 1.42, 1.5, 2, 2.1, 2.2, 2.5, 2.55, 2.6, 2.65, 2.7, 3, 4, 6, 9] # g / cm3

quartz = openmc.Material(name="quartz")
quartz.add_element("Si",1)
quartz.add_element("O", 2)
# quartz.add_s_alpha_beta("c_Si_in_SiO2")
quartz_densities = [0.5, 0.7, 0.85, 1, 1.3, 1.4, 2, 3, 3.5, 4, 4.5, 5, 6, 8, 10, 12, 14, 15, 20, 30, 40, 50, 60, 80, 100, 120, 150, 200]


vanillin = openmc.Material(name="vanillin")
vanillin.add_element("C", 8)
vanillin.add_element("H", 8)
vanillin.add_element("O", 3)
# vanillin.set_density("g")
ethanol = openmc.Material(name="ethanol")
ethanol.add_element("C", 2)
ethanol.add_element("H", 6)
ethanol.add_element("O", 1)
vsyrup = openmc.Material.mix_materials([vanillin, ethanol, glucose], [0.3, 0.4, 0.3], "wo")
vsyrup.name = "vanilla syrup"
vsyrup_densities = [0.5, 0.7, 0.8, 0.9, 1, 1.2, 1.4, 1.6, 2, 2.1, 2.2, 2.4, 3, 4]


vanillin_densities = sorted(vsyrup_densities + [3.25, 3.5, 3.75, 5, 6])

kval_sets = []
moderators = [quartz, corn_syrup, vanillin]

density_sets = {
    quartz:quartz_densities,
    corn_syrup: corn_syrup_densities,
    vsyrup: vsyrup_densities,
    vanillin: vanillin_densities
}

special_pcounts = {
    quartz:4000
}

for mod in moderators:
    densities = density_sets[mod]
    print(mod)
    
    kvals = []
    for density in densities:

        # print(water)

        mod.set_density("g/cm3", density)

        model.materials = openmc.Materials(
            [uo2, zirconium, mod]
        )  # combine materials into the model
        print(model.materials)

        # # create the geometry for the pin cell
        # fuel_outer_radius = openmc.ZCylinder(r=0.46955)
        # clad_inner_radius = openmc.ZCylinder(r=0.47910)
        # clad_outer_radius = openmc.ZCylinder(r=0.54640)

        # fuel_region = -fuel_outer_radius
        # gap_region = +fuel_outer_radius & -clad_inner_radius
        # clad_region = +clad_inner_radius & -clad_outer_radius

        # height = 300.0
        # top = openmc.ZPlane(z0=height/2, boundary_type="vacuum")
        # bot = openmc.ZPlane(z0=-height/2, boundary_type="vacuum")
        # layer = +bot & -top


        # left = openmc.XPlane(-pitch / 2, boundary_type="reflective")
        # right = openmc.XPlane(pitch / 2, boundary_type="reflective")
        # front = openmc.YPlane(-pitch / 2, boundary_type="reflective")
        # back = openmc.YPlane(pitch / 2, boundary_type="reflective")

        water_region = +left & -right & +front & -back & +clad_outer_radius & layer

        moderator = openmc.Cell()
        moderator.fill = mod
        moderator.region = water_region

        root_universe = openmc.Universe(cells=(fuel, gap, clad, moderator))
        model.geometry = openmc.Geometry(root_universe)
        # root_universe.plot(width=(pitch, pitch), basis='xy')


        # Create a point source
        point = openmc.stats.Point((0, 0, 0))
        source = openmc.IndependentSource(space=point)

        model.settings.source = source
        model.settings.batches = 400
        model.settings.inactive = 200

        pcount = 1000
        if mod in special_pcounts:
            pcound = special_pcounts[mod]
        model.settings.particles = pcount

        model.export_to_xml()
        

        p2d_path = f"p2d/sp_{mod.name}_{str(density).replace(".", "p")}.h5"
        if not os.path.exists(p2d_path):
            statepoint = model.run()
            os.rename(statepoint, p2d_path)
        with openmc.StatePoint(p2d_path) as sp:
            kvals.append(sp.keff)
    kval_sets.append(kvals)

for i, mod in enumerate(moderators):
    kvals = kval_sets[i]
    densities = density_sets[mod]

    knoms = [float(k.nominal_value) for k in kvals]
    kstds = [float(k.std_dev) for k in kvals]
    print(f"knoms: {knoms}")
    print(f"kdevs: {kstds}")
    # plt.plot(pcounts, knoms)
    # plt.plot([1, 2, 3], [1, 3, 2], label="test data")
    plt.errorbar(densities, knoms, yerr=kstds, label=f"Moderated by {mod.name}")
plt.xlabel("density (g/cm3)")
plt.ylabel("k_eff")
plt.legend()
plt.show()
