# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------
from math import pi
import numpy as np
from fourier_prop.laser_input import (laser_parameters, propagation_parameters, utils)
from fourier_prop.read_laser import sim_grid_parameters as grid
from fourier_prop.read_laser import read_laser
from fourier_prop.sim_helpers import (foil_shapes, sim_parameters)


def microns_to_norm_units(l):
    return utils.microns_to_norm_units(l, laser_parameters.REF_FREQ)


def fs_to_norm_units(t):
    return utils.fs_to_norm_units(t, laser_parameters.REF_FREQ)


# SIMULATION DIMENSIONS
l0 = 2. * pi  # laser wavelength [in code units]
t0 = l0  # optical cycle
Lsim = [
    microns_to_norm_units(grid.X_LENGTH),
    microns_to_norm_units(grid.Y_HEIGHT),
    microns_to_norm_units(grid.Z_HEIGHT)
]  # length of the simulation
Tsim = fs_to_norm_units(grid.T_LENGTH)

dt = t0 * grid.DT_SIM

sim_grid_parameters = grid.compute_sim_grid(
    propagation_parameters.TIMES,
    propagation_parameters.Y_VALS_OUTPUT,
    propagation_parameters.Z_VALS_OUTPUT
)

by_func = read_laser.get_By_function(propagation_parameters.DATA_DIRECTORY_PATH, sim_grid_parameters)
bz_func = read_laser.get_Bz_function(propagation_parameters.DATA_DIRECTORY_PATH, sim_grid_parameters)

ppc = 16

Main(
    geometry="3Dcartesian",
    solve_poisson=True,

    interpolation_order=2,

    cell_length=[l0 * grid.DX_SIM, l0 * grid.DY_SIM, l0 * grid.DZ_SIM],
    grid_length=Lsim,

    number_of_patches=[16, 16, 16],

    timestep=dt,
    simulation_time=Tsim,
    reference_angular_frequency_SI=laser_parameters.REF_FREQ,  # for ionization

    EM_boundary_conditions=[
        ['silver-muller'],
        ['silver-muller'],
        ['silver-muller']
    ],
)


def foil_shape_profile(n0):
    return foil_shapes.circular_foil(
        n0,
        microns_to_norm_units(sim_parameters.FOIL_LEFT_X),
        microns_to_norm_units(sim_parameters.FOIL_RADIUS),
        microns_to_norm_units(sim_parameters.FOIL_THICKNESS),
        microns_to_norm_units(grid.Y_HEIGHT / 2.),
        microns_to_norm_units(grid.Z_HEIGHT / 2.),
        sim_parameters.PRE_PLASMA_PARAMS["PRE_PLASMA"],
        microns_to_norm_units(sim_parameters.PRE_PLASMA_PARAMS["CHAR_LENGTH"]),
        sim_parameters.PRE_PLASMA_PARAMS["CUT_OFF_DENSITY"]
    )


n0_h = 30
frozen_time = 0
cold_or_mj = 'cold'

Species(
    name='hydrogen_ions',
    position_initialization='random',
    momentum_initialization=cold_or_mj,
    particles_per_cell=ppc,
    mass=1836. * 1,
    charge=1.,
    number_density=foil_shape_profile(
        n0_h,
    ),
    boundary_conditions=[
        ["remove", "remove"],
        ["remove", "remove"],
        ["remove", "remove"],
    ],
    time_frozen=frozen_time,
    temperature=[0.],
)

Species(
    name='hydrogen_electrons',
    position_initialization='hydrogen_ions',
    momentum_initialization=cold_or_mj,
    particles_per_cell=ppc,
    mass=1.,
    charge=-1.,
    number_density=foil_shape_profile(
        n0_h,
    ),
    boundary_conditions=[
        ["remove", "remove"],
        ["remove", "remove"],
        ["remove", "remove"],
    ],
    time_frozen=frozen_time,
    temperature=[0.],
)

Laser(
    box_side="xmin",
    space_time_profile=[by_func, bz_func]
)

##### DIAGNOSTICS #####

period_timestep = t0 / dt
data_sample_rate = 1 * period_timestep

fields = ["Ey", "Ez", "Ex", "Rho_hydrogen_ions", "Rho_hydrogen_electrons", "By", "Bz", "Jx"]

DiagScalar(
    every=10,
    vars=["Utot", "Ukin", "Uelm"],
    precision=10
)

# XY Plane
DiagProbe(
    # name = "my_probe",
    every=5 * data_sample_rate,
    origin=[0., 0., Lsim[2] / 2.],
    corners=[
        [Lsim[0], 0., Lsim[2] / 2.],
        [0., Lsim[1], Lsim[2] / 2.],
    ],
    number=[1000, 1000],
    fields=fields
)

# XZ Plane
DiagProbe(
    # name = "my_probe",
    every=5 * data_sample_rate,
    origin=[0., Lsim[1] / 2., 0.],
    corners=[
        [Lsim[0], Lsim[1] / 2., 0.],
        [0., Lsim[1] / 2., Lsim[2]],
    ],
    number=[1000, 1000],
    fields=fields
)

# XY + 45 deg
DiagProbe(
    # name = "my_probe",
    every=5 * data_sample_rate,
    origin=[0., 0., 0.],
    corners=[
        [Lsim[0], 0., 0.],
        [0., Lsim[1], Lsim[2]],
    ],
    number=[1000, 1000],
    fields=fields
)


# YZ Plane 20um
DiagProbe(
    # name = "my_probe",
    every=5 * data_sample_rate,
    origin=[microns_to_norm_units(20), 0., 0.],
    corners=[
        [microns_to_norm_units(40), Lsim[1], 0.],
        [microns_to_norm_units(40), 0., Lsim[2]],
    ],
    number=[1000, 1000],
    fields=["Rho_hydrogen_ions", "Rho_hydrogen_electrons", "Bz", "By"]
)

# YZ Plane 40um
DiagProbe(
    # name = "my_probe",
    every=5 * data_sample_rate,
    origin=[microns_to_norm_units(40), 0., 0.],
    corners=[
        [microns_to_norm_units(40), Lsim[1], 0.],
        [microns_to_norm_units(40), 0., Lsim[2]],
    ],
    number=[1000, 1000],
    fields=["Rho_hydrogen_ions", "Rho_hydrogen_electrons", "Bz", "By"]
)

# NUMBER DENSITY SCREEN 40, and 47um
DiagScreen(
    # name = "my screen",
    shape="plane",
    point=[microns_to_norm_units(40), 0., 0.],
    vector=[1., 0., 0.],
    direction="both",
    deposited_quantity="weight",
    species=["hydrogen_ions"],
    axes=[["y", 0, Lsim[1], 400], ["z", 0, Lsim[2], 400]],
    every=(Tsim / dt) - 1  # the entire sim
)

DiagScreen(
    # name = "my screen",
    shape="plane",
    point=[microns_to_norm_units(47), 0., 0.],
    vector=[1., 0., 0.],
    direction="both",
    deposited_quantity="weight",
    species=["hydrogen_ions"],
    axes=[["y", 0, Lsim[1], 400], ["z", 0, Lsim[2], 400]],
    every=(Tsim / dt) - 1  # the entire sim
)

DiagScreen(
    # name = "my screen",
    shape="plane",
    point=[microns_to_norm_units(40), 0., 0.],
    vector=[1., 0., 0.],
    direction="both",
    deposited_quantity="weight",
    species=["hydrogen_electrons"],
    axes=[["y", 0, Lsim[1], 400], ["z", 0, Lsim[2], 400]],
    every=(Tsim / dt) - 1  # the entire sim
)

DiagScreen(
    # name = "my screen",
    shape="plane",
    point=[microns_to_norm_units(47), 0., 0.],
    vector=[1., 0., 0.],
    direction="both",
    deposited_quantity="weight",
    species=["hydrogen_electrons"],
    axes=[["y", 0, Lsim[1], 400], ["z", 0, Lsim[2], 400]],
    every=(Tsim / dt) - 1  # the entire sim
)

# 20 and 30
DiagScreen(
    # name = "my screen",
    shape="plane",
    point=[microns_to_norm_units(20), 0., 0.],
    vector=[1., 0., 0.],
    direction="both",
    deposited_quantity="weight",
    species=["hydrogen_ions"],
    axes=[["y", 0, Lsim[1], 400], ["z", 0, Lsim[2], 400]],
    every=(Tsim / dt) - 1  # the entire sim
)

DiagScreen(
    # name = "my screen",
    shape="plane",
    point=[microns_to_norm_units(30), 0., 0.],
    vector=[1., 0., 0.],
    direction="both",
    deposited_quantity="weight",
    species=["hydrogen_ions"],
    axes=[["y", 0, Lsim[1], 400], ["z", 0, Lsim[2], 400]],
    every=(Tsim / dt) - 1  # the entire sim
)

DiagScreen(
    # name = "my screen",
    shape="plane",
    point=[microns_to_norm_units(20), 0., 0.],
    vector=[1., 0., 0.],
    direction="both",
    deposited_quantity="weight",
    species=["hydrogen_electrons"],
    axes=[["y", 0, Lsim[1], 400], ["z", 0, Lsim[2], 400]],
    every=(Tsim / dt) - 1  # the entire sim
)

DiagScreen(
    # name = "my screen",
    shape="plane",
    point=[microns_to_norm_units(30), 0., 0.],
    vector=[1., 0., 0.],
    direction="both",
    deposited_quantity="weight",
    species=["hydrogen_electrons"],
    axes=[["y", 0, Lsim[1], 400], ["z", 0, Lsim[2], 400]],
    every=(Tsim / dt) - 1  # the entire sim
)

# DIVERGENCE SCREEN 30, 40, and 47um
def angle(p):
    r = np.sqrt(p.py ** 2 + p.pz ** 2)
    r = np.where((p.py < 0) ^ (p.pz < 0), -r, r)
    return (180 / np.pi) * np.arctan2(r, p.px)

def angle_y(p):
    return (180 / np.pi) * np.arctan2(p.py, p.px)

def angle_z(p):
    return (180 / np.pi) * np.arctan2(p.pz, p.px)

DiagScreen(
    # name = "my screen",
    shape="plane",
    point=[microns_to_norm_units(30), 0., 0.],
    vector=[1., 0., 0.],
    direction="both",
    deposited_quantity="weight_ekin",
    species=["hydrogen_ions"],
    axes=[[angle_y, -180, 180, 400]],
    every=(Tsim / dt) - 1  # the entire sim
)

DiagScreen(
    # name = "my screen",
    shape="plane",
    point=[microns_to_norm_units(40), 0., 0.],
    vector=[1., 0., 0.],
    direction="both",
    deposited_quantity="weight_ekin",
    species=["hydrogen_ions"],
    axes=[[angle_y, -180, 180, 400]],
    every=(Tsim / dt) - 1  # the entire sim
)

DiagScreen(
    # name = "my screen",
    shape="plane",
    point=[microns_to_norm_units(47), 0., 0.],
    vector=[1., 0., 0.],
    direction="both",
    deposited_quantity="weight_ekin",
    species=["hydrogen_ions"],
    axes=[[angle_y, -180, 180, 400]],
    every=(Tsim / dt) - 1  # the entire sim
)

# ENERGY DENSITY SCREEN 40, and 47um
DiagScreen(
    # name = "my screen",
    shape="plane",
    point=[microns_to_norm_units(40), 0., 0.],
    vector=[1., 0., 0.],
    direction="both",
    deposited_quantity="weight_ekin",
    species=["hydrogen_ions"],
    axes=[["y", 0, Lsim[1], 400], ["z", 0, Lsim[2], 400]],
    every=(Tsim / dt) - 1  # the entire sim
)

DiagScreen(
    # name = "my screen",
    shape="plane",
    point=[microns_to_norm_units(47), 0., 0.],
    vector=[1., 0., 0.],
    direction="both",
    deposited_quantity="weight_ekin",
    species=["hydrogen_ions"],
    axes=[["y", 0, Lsim[1], 400], ["z", 0, Lsim[2], 400]],
    every=(Tsim / dt) - 1  # the entire sim
)

DiagScreen(
    # name = "my screen",
    shape="plane",
    point=[microns_to_norm_units(40), 0., 0.],
    vector=[1., 0., 0.],
    direction="both",
    deposited_quantity="weight_ekin",
    species=["hydrogen_electrons"],
    axes=[["y", 0, Lsim[1], 400], ["z", 0, Lsim[2], 400]],
    every=(Tsim / dt) - 1  # the entire sim
)

DiagScreen(
    # name = "my screen",
    shape="plane",
    point=[microns_to_norm_units(47), 0., 0.],
    vector=[1., 0., 0.],
    direction="both",
    deposited_quantity="weight_ekin",
    species=["hydrogen_electrons"],
    axes=[["y", 0, Lsim[1], 400], ["z", 0, Lsim[2], 400]],
    every=(Tsim / dt) - 1  # the entire sim
)

# 20 and 30
DiagScreen(
    # name = "my screen",
    shape="plane",
    point=[microns_to_norm_units(20), 0., 0.],
    vector=[1., 0., 0.],
    direction="both",
    deposited_quantity="weight_ekin",
    species=["hydrogen_ions"],
    axes=[["y", 0, Lsim[1], 400], ["z", 0, Lsim[2], 400]],
    every=(Tsim / dt) - 1  # the entire sim
)

DiagScreen(
    # name = "my screen",
    shape="plane",
    point=[microns_to_norm_units(30), 0., 0.],
    vector=[1., 0., 0.],
    direction="both",
    deposited_quantity="weight_ekin",
    species=["hydrogen_ions"],
    axes=[["y", 0, Lsim[1], 400], ["z", 0, Lsim[2], 400]],
    every=(Tsim / dt) - 1  # the entire sim
)

DiagScreen(
    # name = "my screen",
    shape="plane",
    point=[microns_to_norm_units(20), 0., 0.],
    vector=[1., 0., 0.],
    direction="both",
    deposited_quantity="weight_ekin",
    species=["hydrogen_electrons"],
    axes=[["y", 0, Lsim[1], 400], ["z", 0, Lsim[2], 400]],
    every=(Tsim / dt) - 1  # the entire sim
)

DiagScreen(
    # name = "my screen",
    shape="plane",
    point=[microns_to_norm_units(30), 0., 0.],
    vector=[1., 0., 0.],
    direction="both",
    deposited_quantity="weight_ekin",
    species=["hydrogen_electrons"],
    axes=[["y", 0, Lsim[1], 400], ["z", 0, Lsim[2], 400]],
    every=(Tsim / dt) - 1  # the entire sim
)

# Spectrum
DiagParticleBinning(
    # name = "my binning",
    deposited_quantity="weight",
    every=1 * data_sample_rate,
    time_average=1,
    species=["hydrogen_ions"],
    axes=[
        ["ekin", "auto", "auto", 400]
    ]
)

DiagParticleBinning(
    # name = "my binning",
    deposited_quantity="weight",
    every=1 * data_sample_rate,
    time_average=1,
    species=["hydrogen_electrons"],
    axes=[
        ["ekin", "auto", "auto", 400]
    ]
)


# Cylinder Bin
def radius(p):
    centery = Lsim[1] / 2.
    centerz = Lsim[2] / 2.
    r = np.sqrt((p.y - centery) ** 2 + (p.z - centerz) ** 2)
    return r


def cylinder_filter(p, r_um):
    radius = microns_to_norm_units(r_um)
    start_x = microns_to_norm_units(20)
    radius_p = np.sqrt((p.y - (Lsim[1]/2.))**2 + (p.z - (Lsim[2]/2.))**2)
    
    return (p.x > start_x) & (radius_p < radius)

def divergence_y(p):
    divergence = np.arctan2(p.py, p.px)
    
    return np.where(cylinder_filter(p, 5), divergence, 1.2*np.pi + 0*p.x)

def divergence_z(p):
    divergence = np.arctan2(p.pz, p.px)
    
    return np.where(cylinder_filter(p, 5), divergence, 1.2*np.pi + 0*p.x)

def divergence_y_small(p):
    divergence = np.arctan2(p.py, p.px)

    return np.where(cylinder_filter(p, 2.5), divergence, 1.2*np.pi + 0*p.x)

def divergence_z_small(p):
    divergence = np.arctan2(p.pz, p.px)

    return np.where(cylinder_filter(p, 2.5), divergence, 1.2*np.pi + 0*p.x)

def divergence_y_tiny(p):
    divergence = np.arctan2(p.py, p.px)

    return np.where(cylinder_filter(p, 1), divergence, 1.2*np.pi + 0*p.x)

def divergence_z_tiny(p):
    divergence = np.arctan2(p.pz, p.px)

    return np.where(cylinder_filter(p, 1), divergence, 1.2*np.pi + 0*p.x)



DiagParticleBinning(
    deposited_quantity="weight",
    every=data_sample_rate,
    time_average=1,
    species=["hydrogen_ions"],
    axes=[
        ["y", 0, Lsim[1], 1000],
        [divergence_y, -np.pi, np.pi, 1000],
    ],
)

DiagParticleBinning(
    deposited_quantity="weight",
    every=data_sample_rate,
    time_average=1,
    species=["hydrogen_ions"],
    axes=[
        ["z", 0, Lsim[2], 1000],
        [divergence_z, -np.pi, np.pi, 1000],
    ],
)

DiagParticleBinning(
    deposited_quantity="weight",
    every=data_sample_rate,
    time_average=1,
    species=["hydrogen_ions"],
    axes=[
        ["y", 0, Lsim[1], 1000],
        [divergence_y_small, -np.pi, np.pi, 1000],
    ],
)

DiagParticleBinning(
    deposited_quantity="weight",
    every=data_sample_rate,
    time_average=1,
    species=["hydrogen_ions"],
    axes=[
        ["z", 0, Lsim[2], 1000],
        [divergence_z_small, -np.pi, np.pi, 1000],
    ],
)

DiagParticleBinning(
    deposited_quantity="weight",
    every=data_sample_rate,
    time_average=1,
    species=["hydrogen_ions"],
    axes=[
        ["y", 0, Lsim[1], 1000],
        [divergence_y_tiny, -np.pi, np.pi, 1000],
    ],
)

DiagParticleBinning(
    deposited_quantity="weight",
    every=data_sample_rate,
    time_average=1,
    species=["hydrogen_ions"],
    axes=[
        ["z", 0, Lsim[2], 1000],
        [divergence_z_tiny, -np.pi, np.pi, 1000],
    ],
)


# Spectrum of Cylinders
def energy_cylinder(p):
    energy = 1836 * (np.sqrt(1 + p.pz**2 + p.px**2 + p.py**2) - 1) * 0.511

    return np.where(cylinder_filter(p, 5), energy, -1 + 0*p.x)

def energy_cylinder_small(p):
    energy = 1836 * (np.sqrt(1 + p.pz**2 + p.px**2 + p.py**2) - 1) * 0.511

    return np.where(cylinder_filter(p, 2.5), energy, -1 + 0*p.x)

def energy_cylinder_tiny(p):
    energy = 1836 * (np.sqrt(1 + p.pz**2 + p.px**2 + p.py**2) - 1) * 0.511

    return np.where(cylinder_filter(p, 1), energy, -1 + 0*p.x)

def energy_cylinder_e(p):
    energy = (np.sqrt(1 + p.pz**2 + p.px**2 + p.py**2) - 1) * 0.511

    return np.where(cylinder_filter(p, 5), energy, -1 + 0*p.x)

def energy_cylinder_small_e(p):
    energy = (np.sqrt(1 + p.pz**2 + p.px**2 + p.py**2) - 1) * 0.511

    return np.where(cylinder_filter(p, 2.5), energy, -1 + 0*p.x)

def energy_cylinder_tiny_e(p):
    energy = (np.sqrt(1 + p.pz**2 + p.px**2 + p.py**2) - 1) * 0.511

    return np.where(cylinder_filter(p, 1), energy, -1 + 0*p.x)



DiagParticleBinning(
    deposited_quantity="weight",
    every=data_sample_rate,
    time_average=1,
    species=["hydrogen_ions"],
    axes=[
        [energy_cylinder, 0, 220, 3000],
    ],
)

DiagParticleBinning(
    deposited_quantity="weight",
    every=data_sample_rate,
    time_average=1,
    species=["hydrogen_ions"],
    axes=[
        [energy_cylinder_small, 0, 220, 3000],
    ],
)

DiagParticleBinning(
    deposited_quantity="weight",
    every=data_sample_rate,
    time_average=1,
    species=["hydrogen_ions"],
    axes=[
        [energy_cylinder_tiny, 0, 220, 3000],
    ],
)

DiagParticleBinning(
    deposited_quantity="weight",
    every=data_sample_rate,
    time_average=1,
    species=["hydrogen_electrons"],
    axes=[
        [energy_cylinder_e, 0, 220, 3000],
    ],
)

DiagParticleBinning(
    deposited_quantity="weight",
    every=data_sample_rate,
    time_average=1,
    species=["hydrogen_electrons"],
    axes=[
        [energy_cylinder_small_e, 0, 220, 3000],
    ],
)

DiagParticleBinning(
    deposited_quantity="weight",
    every=data_sample_rate,
    time_average=1,
    species=["hydrogen_electrons"],
    axes=[
        [energy_cylinder_tiny_e, 0, 220, 3000],
    ],
)

