"""
Unit conversion factors for TauRunner.jl

All conversions are to/from natural units (eV-based system where c = hbar = 1).
The convention is: quantity_in_eV = quantity_in_unit * UNIT_FACTOR

For example:
    energy_eV = energy_GeV * GeV
    length_eV_inv = length_km * km
"""
module units

# Energy conversions [eV per unit]
const TeV = 1.0e12
const GeV = 1.0e9
const MeV = 1.0e6
const keV = 1.0e3
const eV = 1.0
const Joule = 1 / 1.60225e-19

# Mass conversions [eV per unit]
const kg = 5.62e35
const gr = 1e-3 * kg

# Time conversions [eV^-1 per unit]
const sec = 1.523e15
const hour = 3600.0 * sec
const day = 24.0 * hour
const year = 365.0 * day
const yearstosec = sec / year

# Distance conversions [eV^-1 per unit]
const meter = 5.06773093741e6
const cm = 1.0e-2 * meter
const km = 1.0e3 * meter
const fermi = 1.0e-15 * meter
const angstrom = 1.0e-10 * meter
const AU = 149.60e9 * meter
const parsec = 3.08568025e16 * meter

# Area/Cross-section conversions [eV^-2 per unit]
const picobarn = 1.0e-36 * cm^2
const femtobarn = 1.0e-39 * cm^2

# Pressure conversions [eV^4 per unit]
const Pascal = Joule / meter^3
const hPascal = 100.0 * Pascal
const atm = 101325.0 * Pascal
const psi = 6893.0 * Pascal

# Temperature conversions [eV per unit]
const kelvin = 1 / 1.1604505e4

# Angle conversions [radians per unit]
const degree = pi / 180.0

# Magnetic field [eV^2 per unit]
const T = 0.000692445  # Tesla

# Legacy notation (for compatibility)
const cm3toev3 = 7.68351405e-15
const KmtoEv = 5.0677288532e9

# Derived quantities commonly used
"""Density conversion factor: gr/cm^3 to eV^4"""
const DENSITY_CONV = gr / cm^3

end # module units
