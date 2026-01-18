"""
Physics constants for TauRunner.jl

Based on the PhysicsConstants class from the Python implementation.
All values in natural units (eV-based) unless otherwise noted.
"""

"""
    PhysicsConstants

Module containing physical constants used in neutrino propagation simulations.
"""
module PhysicsConstants

# Mathematical constants
const PI = 3.14159265
const PIby2 = 1.5707963268
const SQRT2 = 1.4142135624
const LN2 = log(2.0)

# Earth and Sun radii [km]
const EARTH_RADIUS_KM = 6371.0
const SUN_RADIUS_KM = 109 * EARTH_RADIUS_KM

# Fundamental constants
const GF = 1.16639e-23           # Fermi constant [eV^-2]
const NA = 6.0221415e23          # Avogadro number [mol^-1]
const SW_SQ = 0.2312             # sin^2(theta_weinberg)
const G = 6.67300e-11            # Gravitational constant [m^3 kg^-1 s^-2]
const ALPHA = 1.0/137.0          # Fine-structure constant

# Particle masses [GeV]
const MUON_MASS_GEV = 0.10565
const NEUTRON_MASS_GEV = 0.939565
const PROTON_MASS_GEV = 0.938272
const ELECTRON_MASS_GEV = 0.510998910e-3
const TAU_MASS_GEV = 1.776

# Particle lifetimes [seconds]
const TAU_LIFETIME_SEC = 2.9e-13
const MUON_LIFETIME_SEC = 2.2e-6

# Atomic mass unit [g]
const ATOMIC_MASS_UNIT = 1.660538e-24

# Isoscalar mass (average of proton and neutron) [GeV]
const ISOSCALAR_MASS_GEV = (PROTON_MASS_GEV + NEUTRON_MASS_GEV) / 2

# Neutrino oscillation parameters
# Mixing angles [radians]
const TH12 = 0.579639
const TH13 = 0.150098
const TH23 = PIby2 / 2.0

# CP phase
const DELTA_CP = 5.235987

# Squared mass differences [eV^2]
const DM21SQ = 7.50e-5
const DM31SQ = 2.47e-3
const DM32SQ = -2.43e-3

# Chemistry: element -> (atomic_number, atomic_mass_amu)
const CHEMISTRY = Dict(
    "H"  => (1,  1.00784),
    "Si" => (14, 28.0855),
    "Mg" => (12, 24.305),
    "Ni" => (28, 58.6934),
    "Fe" => (56, 55.845),
    "O"  => (8,  15.999),
    "Al" => (13, 26.9815),
    "Ca" => (20, 40.078),
    "Na" => (11, 22.9897),
    "K"  => (19, 39.0983),
    "Ti" => (22, 47.867),
)

end # module PhysicsConstants
