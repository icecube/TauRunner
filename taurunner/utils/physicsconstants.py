"""
Author  : C.A. Arguelles
Date    : 10/MAY/2011

Contains Physics constants and global variables.

Log :
- Modified on 23/ABR/2012 by C.Arguelles
    + Changed the definition of PhysicsConstants to
    include an __init__ to separate the class global
    properties from its instances.
"""

# python standard modules
import numpy as np

class PhysicsConstants():
    
    def __init__(self):
        ## PHYSICS CONSTANTS
        #===========================================================================
        # NAME
        #===========================================================================
        
        self.name = "STD"                    # Default values
        self.linestyle = "solid"             # Default linestyle in plots
        self.markerstyle = "*"               # Default marker style
        self.colorstyle = "red"              # Default color style
        self.savefilename = "output.dat"     # Default color style
        
        #===============================================================================
        # ## MATH
        #===============================================================================
        self.PI=3.14159265		            # Pi
        self.PIby2=1.5707963268	            # Pi/2
        self.sqr2=1.4142135624	            # Sqrt[2]
        self.ln2 = np.log(2.0)
        
        #===============================================================================
        # ## EARTH 
        #===============================================================================
        self.EARTHRADIUS = 6371.0	        # [km] Earth radius
        #===============================================================================
        # ## SUN 
        #===============================================================================
        self.SUNRADIUS = 109*self.EARTHRADIUS     # [km] Sun radius 
        
        #===============================================================================
        # # PHYSICAL CONSTANTS
        #===============================================================================
        self.GF = 1.16639e-23	            # [eV^-2] Fermi Constant 
        self.Na = 6.0221415e+23 		        # [mol cm^-3] Avogadro Number
        self.sw_sq = 0.2312                  # [dimensionless] sin(th_weinberg) ^2
        self.G  = 6.67300e-11                # [m^3 kg^-1 s^-2]
        self.alpha = 1.0/137.0               # [dimensionless] fine-structure constant 
        
        #===============================================================================
        # ## UNIT CONVERSION FACTORS
        #===============================================================================
        # Energy
        self.TeV = 1.0e12                    # [eV/TeV]
        self.GeV = 1.0e9                     # [eV/GeV]
        self.MeV = 1.0e6                     # [eV/MeV]
        self.keV = 1.0e3                     # [eV/keV]
        self.Joule = 1/1.60225e-19           # [eV/J]
        # Mass
        self.kg = 5.62e35                    # [eV/kg]
        self.gr = 1e-3*self.kg               # [eV/g] 
        # Time
        self.sec = 1.523e15                  # [eV^-1/s]
        self.hour = 3600.0*self.sec          # [eV^-1/h]
        self.day = 24.0*self.hour            # [eV^-1/d]
        self.year = 365.0*self.day           # [eV^-1/yr]
        self.yearstosec = self.sec/self.year # [s/yr]
        # Distance
        self.meter = 5.06773093741e6         # [eV^-1/m]
        self.cm = 1.0e-2*self.meter          # [eV^-1/cm]
        self.km = 1.0e3*self.meter           # [eV^-1/km]
        self.fermi = 1.0e-15*self.meter      # [eV^-1/fm]
        self.angstrom = 1.0e-10*self.meter   # [eV^-1/A]
        self.AU = 149.60e9*self.meter        # [eV^-1/AU]
        self.parsec = 3.08568025e16*self.meter# [eV^-1/parsec]
        # Integrated Luminocity # review
        self.picobarn = 1.0e-36*self.cm**2   # [eV^-2/pb]
        self.femtobarn = 1.0e-39*self.cm**2  # [eV^-2/fb]
        # Presure
        self.Pascal = self.Joule/self.meter**3 # [eV^4/Pa]
        self.hPascal = 100.0*self.Pascal     # [eV^4/hPa]
        self.atm = 101325.0*self.Pascal      # [eV^4/atm]
        self.psi = 6893.0*self.Pascal        # [eV^4/psi]
        # Temperature
        self.kelvin = 1/1.1604505e4          # [eV/K]
        # Angle
        self.degree = self.PI/180.0          # [rad/degree]
        # magnetic field
        self.T = 0.000692445                 # [eV^2/T]
        
        # old notation
        self.cm3toev3 = 7.68351405e-15       # cm^3-> ev^3
        self.KmtoEv =5.0677288532e+9         # km -> eV
        self.yearstosec = 31536.0e3          # years -> sec
        
        #===============================================================================
        # ## NEUTRINO OSCILLATION PARAMETERS ##
        #===============================================================================
        
        self.numneu = 3                      # number of neutrinos
        self.numneumax = 6                   # maximum neutrino number
        self.neutype = 'neutrino'
        #neutype = 'antineutrino'
        
        # values updated according to 1209.3023 Table 1 FreeFluxes + RSBL
        
        # MIXING ANGLES
       
        self.th12 = 0.579639
        self.th13 = 0.150098
        self.th23 = self.PIby2/2.0
        self.th14 = 0.0
        self.th24 = 0.0
        self.th34 = 0.0
        self.th15 = 0.0
        self.th25 = 0.0
        self.th35 = 0.0
        self.th45 = 0.0
        self.th16 = 0.0
        self.th26 = 0.0
        self.th36 = 0.0
        self.th46 = 0.0
        self.th56 = 0.0
        
        # mixing angles matrix array
        self.th = np.zeros([self.numneumax+1,self.numneumax+1],float)
        self.th[1,2] = self.th12
        self.th[1,3] = self.th13
        self.th[2,3] = self.th23
        self.th[1,4] = self.th14
        self.th[2,4] = self.th24
        self.th[3,4] = self.th34
        self.th[1,5] = self.th15
        self.th[2,5] = self.th25
        self.th[3,5] = self.th35
        self.th[4,5] = self.th45
        self.th[1,6] = self.th16
        self.th[2,6] = self.th26
        self.th[3,6] = self.th36
        self.th[4,6] = self.th46
        self.th[5,6] = self.th56
        
        self.s12 = np.sin(self.th12)
        self.c12 = np.cos(self.th12)
        self.s13 = np.sin(self.th13)
        self.c13 = np.cos(self.th13)
        self.s23 = np.sin(self.th23)
        self.c23 = np.cos(self.th23)
        self.s14 = np.sin(self.th14)
        self.c14 = np.cos(self.th14)
        self.s24 = np.sin(self.th24)
        self.c24 = np.cos(self.th24)
        self.s34 = np.sin(self.th34)
        self.c34 = np.cos(self.th34)
        self.s15 = np.sin(self.th15)
        self.c15 = np.cos(self.th15)
        self.s25 = np.sin(self.th25)
        self.c25 = np.cos(self.th25)
        self.s35 = np.sin(self.th35)
        self.c35 = np.cos(self.th35)
        self.s45 = np.sin(self.th45)
        self.c45 = np.cos(self.th45)
        self.s16 = np.sin(self.th16)
        self.c16 = np.cos(self.th16)
        self.s26 = np.sin(self.th26)
        self.c26 = np.cos(self.th26)
        self.s36 = np.sin(self.th36)
        self.c36 = np.cos(self.th36)
        self.s46 = np.sin(self.th46)
        self.c46 = np.cos(self.th46)
        self.s56 = np.sin(self.th56)
        self.c56 = np.cos(self.th56)    
        
        # cos(th_ij) matrix array
        self.c = np.zeros([self.numneumax+1,self.numneumax+1],float)
        self.c[1,2] = self.c12
        self.c[1,3] = self.c13
        self.c[1,4] = self.c14
        self.c[2,3] = self.c23
        self.c[2,4] = self.c24
        self.c[3,4] = self.c34
        self.c[1,5] = self.c15
        self.c[2,5] = self.c25
        self.c[3,5] = self.c35
        self.c[4,5] = self.c45
        self.c[1,6] = self.c16
        self.c[2,6] = self.c26
        self.c[3,6] = self.c36
        self.c[4,6] = self.c46
        self.c[5,6] = self.c56    
        
        # sin(th_ij) matrix array
        self.s = np.zeros([self.numneumax+1,self.numneumax+1],float)
        self.s[1,2] = self.s12
        self.s[1,3] = self.s13
        self.s[1,4] = self.s14
        self.s[2,3] = self.s23
        self.s[2,4] = self.s24
        self.s[3,4] = self.s34
        self.s[1,5] = self.s15
        self.s[2,5] = self.s25
        self.s[3,5] = self.s35
        self.s[4,5] = self.s45
        self.s[1,6] = self.s16
        self.s[2,6] = self.s26
        self.s[3,6] = self.s36
        self.s[4,6] = self.s46
        self.s[5,6] = self.s56    
        
        # CP PHASES
        #self.delta21=3.3e-5
        #self.delta32=3.1e-3
        #self.delta31=self.delta32+self.delta21
        #self.deltaCP=self.PIby2
        
        # CP Phases
        self.deltaCP = 5.235987
        self.delta1 = self.deltaCP
        self.delta2 = 0.0
        self.delta3 = 0.0
        
        # d-CP phases
        self.dcp = np.zeros([self.numneumax-2+1],complex)
        self.dcp[0] = 1.0
        self.dcp[1] = self.delta1
        self.dcp[2] = self.delta2
        self.dcp[3] = self.delta3
        
        # SQUARED MASS DIFFERENCE
        self.dm21sq = 7.50e-5	             # [eV^2]
        self.dm31sq = 2.47e-3	             # [eV^2]
        self.dm32sq = -2.43e-3               # [eV^2]
        # STERILE 
        self.dm41sq = 0.0                    # [eV^2]
        self.dm51sq = 0.0                    # [eV^2]
        self.dm61sq = 0.0                    # [eV^2]
        # SQUARED MASS DIFFERENCE MATRIX
        self.dmsq = np.zeros([self.numneumax+2],float)
        self.dmsq[2] = self.dm21sq
        self.dmsq[3] = self.dm31sq
        self.dmsq[4] = self.dm41sq
        self.dmsq[5] = self.dm51sq
        self.dmsq[6] = self.dm61sq    
        
        self.dm2 = np.zeros([self.numneumax+1,self.numneumax+1],float)
        self.dm2[1,2] = self.dm21sq
        self.dm2[1,3] = self.dm31sq
        self.dm2[2,3] = self.dm32sq
        self.dm2[1,4] = self.dm41sq
        self.dm2[1,5] = self.dm51sq
        self.dm2[1,6] = self.dm61sq    
        
        # MIXING MATRIX
        self.U = None
        
        #===============================================================================
        # # PARTICLE MASSES
        #===============================================================================
        self.muon_mass = 0.10565 	     # [GeV]
        self.neutron_mass = 0.939565         # [GeV]
        self.proton_mass = 0.938272          # [GeV]
        self.electron_mass = 0.510998910e-3  # [GeV]
        
        self.atomic_mass_unit = 1.660538e-24 # [g]
        
        ## names
        self.electron = 0
        self.muon = 1
        self.tau = 2
        self.sterile1 = 3
        self.sterile2 = 4
        self.sterile3 = 5
    
    #===============================================================================
    # REFRESH
    #===============================================================================
    
    def Refresh(self):
        # Refresh angles
        self.s12 = np.sin(self.th12)
        self.c12 = np.cos(self.th12)
        self.s13 = np.sin(self.th13)
        self.c13 = np.cos(self.th13)
        self.s23 = np.sin(self.th23)
        self.c23 = np.cos(self.th23)
        self.s14 = np.sin(self.th14)
        self.c14 = np.cos(self.th14)
        self.s24 = np.sin(self.th24)
        self.c24 = np.cos(self.th24)
        self.s34 = np.sin(self.th34)
        self.c34 = np.cos(self.th34)
        self.s15 = np.sin(self.th15)
        self.c15 = np.cos(self.th15)
        self.s25 = np.sin(self.th25)
        self.c25 = np.cos(self.th25)
        self.s35 = np.sin(self.th35)
        self.c35 = np.cos(self.th35)
        self.s45 = np.sin(self.th45)
        self.c45 = np.cos(self.th45)
        self.s16 = np.sin(self.th16)
        self.c16 = np.cos(self.th16)
        self.s26 = np.sin(self.th26)
        self.c26 = np.cos(self.th26)
        self.s36 = np.sin(self.th36)
        self.c36 = np.cos(self.th36)
        self.s46 = np.sin(self.th46)
        self.c46 = np.cos(self.th46)
        self.s56 = np.sin(self.th56)
        self.c56 = np.cos(self.th56)                
        
        th = self.th
        th[1,2] = self.th12
        th[1,3] = self.th13
        th[2,3] = self.th23
        th[1,4] = self.th14
        th[2,4] = self.th24
        th[3,4] = self.th34
        th[1,5] = self.th15
        th[2,5] = self.th25
        th[3,5] = self.th35
        th[4,5] = self.th45
        th[1,6] = self.th16
        th[2,6] = self.th26
        th[3,6] = self.th36
        th[4,6] = self.th46
        th[5,6] = self.th56        
        # Refresh cos(th_ij)
        c = self.c
        c[1,2] = self.c12
        c[1,3] = self.c13
        c[1,4] = self.c14
        c[2,3] = self.c23
        c[2,4] = self.c24
        c[3,4] = self.c34
        c[1,5] = self.c15
        c[2,5] = self.c25
        c[3,5] = self.c35
        c[4,5] = self.c45
        c[1,6] = self.c16
        c[2,6] = self.c26
        c[3,6] = self.c36
        c[4,6] = self.c46
        c[5,6] = self.c56        
        # Refresh sin(th_ij)
        s = self.s
        self.s[1,2] = self.s12
        self.s[1,3] = self.s13
        self.s[1,4] = self.s14
        self.s[2,3] = self.s23
        self.s[2,4] = self.s24
        self.s[3,4] = self.s34
        self.s[1,5] = self.s15
        self.s[2,5] = self.s25
        self.s[3,5] = self.s35
        self.s[4,5] = self.s45
        self.s[1,6] = self.s16
        self.s[2,6] = self.s26
        self.s[3,6] = self.s36
        self.s[4,6] = self.s46
        self.s[5,6] = self.s56                
        # Refresh CP-Phases
        dcp = self.dcp
        dcp[0] = 1.0
        dcp[1] = self.delta1
        dcp[2] = self.delta2
        dcp[3] = self.delta3
        #dcp[4] = self.delta2        
        # Refresh Square mass differences
        dmsq = self.dmsq
        dmsq[2] = self.dm21sq
        dmsq[3] = self.dm31sq
        dmsq[4] = self.dm41sq
        dmsq[5] = self.dm51sq
        dmsq[6] = self.dm61sq        
        
        dm2 = self.dm2
        dm2[1,2] = self.dm21sq
        dm2[1,3] = self.dm31sq
        dm2[2,3] = self.dm32sq
        dm2[1,4] = self.dm41sq
        dm2[1,5] = self.dm51sq
        dm2[1,6] = self.dm61sq
        
    class MemoryLoadObjects():
        """ Stores : variables, arrays, etc, on memory for faster use.
        """
        body_name     = ""
        body_rdensity = []
        body_ry       = []
        body_rnc      = []
        
#B temporary memory variables #
# USED IN DM.py
act_channel     = 0.0
act_DM_mass     = 0.0
act_neuflavor   = 0.0
act_inter       = 0.0
flag_inter      = False
intershadow_neu  = 0.0
intershadow_aneu = 0.0
# USED IN XSECTIONS.py
E_act           = 0.0
E_NC_act           = 0.0
E_CC_act           = 0.0
act_dsde_NC_inter  = 0.0
act_dsde_NC_n_inter  = 0.0
act_dsde_NC_a_inter  = 0.0
act_sig_NC_n_inter  = 0.0
act_sig_NC_a_inter  = 0.0
act_dsde_CC_inter  = 0.0

act_dsde_CCe_inter  = 0.0
act_dsde_CCm_inter  = 0.0
act_dsde_CCt_inter  = 0.0

#neutrino
act_dsde_CCe_n_inter  = 0.0
act_dsde_CCm_n_inter  = 0.0
act_dsde_CCt_n_inter  = 0.0
#antineutrino
act_dsde_CCe_a_inter  = 0.0
act_dsde_CCm_a_inter  = 0.0
act_dsde_CCt_a_inter  = 0.0

#neutrino
act_sig_CCe_n_inter  = 0.0
act_sig_CCm_n_inter  = 0.0
act_sig_CCt_n_inter  = 0.0
#antineutrino
act_sig_CCe_a_inter  = 0.0
act_sig_CCm_a_inter  = 0.0
act_sig_CCt_a_inter  = 0.0

# USED IN NEUOSC
tauola_init = False
pabs_n_inter = 0.0
pabs_a_inter = 0.0
# neutrino
p_sur_inter_array_en = 0.0
p_sur_inter_array_mn = 0.0
p_sur_inter_array_tn = 0.0
# antineutrino
p_sur_inter_array_ea = 0.0
p_sur_inter_array_ma = 0.0
p_sur_inter_array_ta = 0.0
#E temporary memory variables # 
        
        
