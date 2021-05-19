#cross section module
import os
import sys
info = sys.version_info
pyv  = int(info.major)
import numpy as np
from taurunner.modules import units


class CrossSections(object):

    def __init__(self, model):

        self.model = model
        #cross section tables
        ######################################
        self.cross_section_path = os.path.dirname(os.path.realpath(__file__))+'/cross_section_tables/'
        if(self.model=='dipole'):
            if(pyv==3):
                f_NC = np.load(cross_section_path+'NC_table_py3.npy', allow_pickle=True).item()
                f_CC = np.load(cross_section_path+'CC_table_py3.npy', allow_pickle=True).item()

                dsdy_spline_CC = np.load(cross_section_path + 'dsigma_dy_CC_py3.npy', allow_pickle=True).item()
                dsdy_spline_CC_lowe = np.load(cross_section_path + 'dsigma_dy_CC_lowE_py3.npy', allow_pickle=True).item()
                dsdy_spline_NC = np.load(cross_section_path + 'dsigma_dy_NC_py3.npy', allow_pickle=True).item()
                dsdy_spline_NC_lowe = np.load(cross_section_path + 'dsigma_dy_NC_lowE_py3.npy', allow_pickle=True).item()
            else:
                f_NC = np.load(cross_section_path+'NC_table.npy', allow_pickle=True).item()
                f_CC = np.load(cross_section_path+'CC_table.npy', allow_pickle=True).item()

                dsdy_spline_CC = np.load(cross_section_path + 'dsigma_dy_CC.npy', allow_pickle=True).item()
                dsdy_spline_CC_lowe = np.load(cross_section_path + 'dsigma_dy_CC_lowE.npy', allow_pickle=True).item()
                dsdy_spline_NC = np.load(cross_section_path + 'dsigma_dy_NC.npy', allow_pickle=True).item()
                dsdy_spline_NC_lowe = np.load(cross_section_path + 'dsigma_dy_NC_lowE.npy', allow_pickle=True).item()
        elif(self.model=='CSMS'):
                f_NC = 
                f_CC = 

                dsdy_spline_CC = 
                dsdy_spline_CC_lowe = dsdy_spline_CC
                dsdy_spline_NC =
                dsdy_spline_NC_lowe = dsdy_spline_NC
         ##########################
         ### Load your tables ####

    def TotalNeutrinoCrossSection(self, enu, interaction = 'NC'):
        r'''
        Calculates total neutrino cross section. returns the value of sigma_CC (or NC)
        in natural units.
        Parameters
        ----------
        enu:         float
            neutrino energy in eV
        interaction: str
            string defining the interaction type (CC or NC). default is NC
        Returns
        -------
        TotalCrossSection: float
            Total neutrino cross section at the given energy in natural units.
        '''
        if(np.log10(enu) < 0.):
            raise ValueError("Going below a GeV. this region is not supported.")
        if(interaction == 'NC'):
            return((10**f_NC(np.log10(enu/1e9)))*(units.cm)**2)
        else:
            return((10**f_CC(np.log10(enu/1e9)))*(units.cm)**2)

    def DifferentialOutGoingLeptonDistribution(self, ein, eout, interaction):
        r'''
        Calculates Differential neutrino cross section. returns the value of d$\sigma$/dy
        in natural units.
        Parameters
        ----------
        ein:         float
            incoming lepton energy in GeV
        eout:         float
            outgoing lepton energy in GeV
        interaction: nusquids obj
            string defining the interaction type (CC or NC).
        Returns
        -------
        diff: float
            d$\sigma$/d(1-y) at the given energies in natural units where y is the bjorken-y.
        '''
        if(np.log10(ein) < 0):
            diff = 0
            return diff
        if(interaction=='CC'):
            if (np.log10(eout) >= np.log10(ein)):
                diff = 0.
                return diff
            elif(np.log10(eout) < 3.):
                diff = 10**dsdy_spline_CC_lowe(np.log10(ein), np.log10(eout))[0][0]/ein 
            else:
                diff = 10**dsdy_spline_CC(np.log10(ein), np.log10(eout))[0][0]/ein
        elif(interaction=='NC'):
            if (np.log10(eout) >= np.log10(ein)):
                diff = 0.
                return diff
            elif( np.log10(eout) <= 3.):
                diff = 10**dsdy_spline_NC_lowe(np.log10(ein), np.log10(eout))[0][0]/ein
            else:
                diff = 10**dsdy_spline_NC(np.log10(ein), np.log10(eout))[0][0]/ein
#        elif(self.model=='CSMS'):
#            diff = dis.SingleDifferentialCrossSection(ein*units.GeV, eout*units.GeV, 
#                        nsq.NeutrinoCrossSections_NeutrinoFlavor.tau, 
#                        nsq.NeutrinoCrossSections_NeutrinoType.neutrino,
#                        getattr(nsq.NeutrinoCrossSections_Current, interaction) 
#                       )
        return diff
