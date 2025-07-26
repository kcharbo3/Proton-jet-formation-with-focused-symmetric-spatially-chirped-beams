from fourier_prop.laser_input import constants
from dataclasses import dataclass
import numpy as np

WAVELENGTH = 0.8  # um
REF_FREQ = (2*np.pi*constants.C_SPEED) / (WAVELENGTH * 1.e-6)
OMEGA0 = REF_FREQ * 1e-15  # rad / PHz

POLARIZATION = constants.LINEAR_Y

SPATIAL_SHAPE = constants.GAUSSIAN
SPATIAL_GAUSSIAN_ORDER = 1
TEMPORAL_SHAPE = constants.GAUSSIAN_T
TEMPORAL_GAUSSIAN_ORDER = 1
PHASE_OFFSET = 0.

WAIST_IN = 88.4374e4
DELTAX = 0. * WAIST_IN

PULSE_FWHM = 30.
SPOT_SIZE = 0.904593714644869
OUTPUT_DISTANCE_FROM_FOCUS = -37.25073157017444

NORMALIZE_TO_A0 = False
PEAK_A0 = 21.
TOTAL_ENERGY = 203144608.7520115

# Used for Special Laser Shapes
L = 1  # LG
# Petal Beam Parameters
NUM_PETALS = 4
WAIST_IN_RADIAL = 20e4
WAIST_IN_AZIMUTHAL = (20. / 1.438399455970175) * 1e4

@dataclass
class LaserParameters:
    wavelength: float
    ref_freq: float
    omega0: float
    polarization: str
    spatial_shape: str
    spatial_gaussian_order: int
    temporal_shape: str
    temporal_gaussian_order: int
    phase_offset: float
    use_grating_eq: bool
    alpha: float
    grating_separation: float
    deltax: float
    pulse_fwhm: float
    spot_size: float
    waist_in: float
    output_distance_from_focus: float
    normalize_to_a0: bool
    peak_a0: float
    total_energy: float
    l: int
    num_petals: int
    waist_in_radial: float
    waist_in_azimuthal: float


laser_parameters_obj = LaserParameters(
    wavelength=WAVELENGTH, ref_freq=REF_FREQ, omega0=OMEGA0, polarization=POLARIZATION,
    spatial_shape=SPATIAL_SHAPE, spatial_gaussian_order=SPATIAL_GAUSSIAN_ORDER, temporal_shape=TEMPORAL_SHAPE,
    temporal_gaussian_order=TEMPORAL_GAUSSIAN_ORDER, phase_offset=PHASE_OFFSET, deltax=DELTAX, pulse_fwhm=PULSE_FWHM, 
    spot_size=SPOT_SIZE, waist_in=WAIST_IN, output_distance_from_focus=OUTPUT_DISTANCE_FROM_FOCUS, 
    normalize_to_a0=NORMALIZE_TO_A0, peak_a0=PEAK_A0, total_energy=TOTAL_ENERGY, l=L, num_petals=NUM_PETALS, 
    waist_in_radial=WAIST_IN_RADIAL, waist_in_azimuthal=WAIST_IN_AZIMUTHAL
)
