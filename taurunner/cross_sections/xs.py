import os, sys, pickle
import numpy as np
from importlib.resources import path

from taurunner.modules import units

def hima_tot_xs(E, spl): # pragma: no cover
    pass

def jeff_tot_xs(E, spl):
    return np.exp(spl(np.log(E/1e9)))

def hima_diff_xs(E_in, E_out, spl): # pragma: no cover
    return(10**spl(np.log10(E_in), np.log10(E_out))[0][0]/E_in)

def jeff_diff_xs(E_in, E_out, spl):
    #E_min = np.power(10, spl.extents[0][0])
    E_min = 1 # Lowest knot on spline in GeV
    z     = (E_out-E_min)/(E_in-E_min)
    res = np.power(10,spl(np.log10(E_in),z)[0])/E_in
    return res

def get_file_path(name):
    with path('taurunner.resources.cross_section_tables', name) as p:
        return str(p)

class CrossSections(object):

    def __init__(self, model):

        self.model = model

        #self.cross_section_path = os.path.dirname(os.path.realpath(__file__))+'/cross_section_tables/'
        #cross_section_path      = self.cross_section_path
        if(self.model=='dipole'):
            self.f_NC = np.load(get_file_path('NC_table_py3.npy'), allow_pickle=True).item()
            self.f_CC = np.load(get_file_path('CC_table_py3.npy'), allow_pickle=True).item()

            self.dsdy_spline_CC = np.load(get_file_path('dsigma_dy_CC_py3.npy'), allow_pickle=True).item()
            self.dsdy_spline_CC_lowe = np.load(get_file_path('dsigma_dy_CC_lowE_py3.npy'), allow_pickle=True).item()
            self.dsdy_spline_NC = np.load(get_file_path('dsigma_dy_NC_py3.npy'), allow_pickle=True).item()
            self.dsdy_spline_NC_lowe = np.load(get_file_path('dsigma_dy_NC_lowE_py3.npy'), allow_pickle=True).item()
        elif(self.model=='CSMS'):
            with open(get_file_path('nu_n_dsde_CC.pkl'), 'rb') as pkl_f:
                diff_nu_n_CC = pickle.load(pkl_f)
            with open(get_file_path('nu_p_dsde_CC.pkl'), 'rb') as pkl_f:
                diff_nu_p_CC = pickle.load(pkl_f)
            with open(get_file_path('nu_n_dsde_NC.pkl'), 'rb') as pkl_f:
                diff_nu_n_NC = pickle.load(pkl_f)
            with open(get_file_path('nu_n_dsde_NC.pkl'), 'rb') as pkl_f:
                diff_nu_p_NC = pickle.load(pkl_f)
            with open(get_file_path('nu_n_sigma_CC.pkl'), 'rb') as pkl_f:
                nu_n_CC = pickle.load(pkl_f)
            with open(get_file_path('nu_p_sigma_CC.pkl'), 'rb') as pkl_f:
                nu_p_CC = pickle.load(pkl_f)
            with open(get_file_path('nu_n_sigma_NC.pkl'), 'rb') as pkl_f:
                nu_n_NC = pickle.load(pkl_f)
            with open(get_file_path('nu_p_sigma_NC.pkl'), 'rb') as pkl_f:
                nu_p_NC = pickle.load(pkl_f)
            self._nu_p_NC = lambda E:jeff_tot_xs(E,nu_p_NC)
            self._nu_n_NC = lambda E:jeff_tot_xs(E,nu_n_NC)
            self._nu_p_CC = lambda E:jeff_tot_xs(E,nu_p_CC)
            self._nu_n_CC = lambda E:jeff_tot_xs(E,nu_n_CC)
            self.f_NC = lambda E: (jeff_tot_xs(E, nu_p_NC)+jeff_tot_xs(E, nu_n_NC))/2.
            self.f_CC = lambda E: (jeff_tot_xs(E, nu_p_CC)+jeff_tot_xs(E, nu_n_CC))/2.
                
            self.dsdy_spline_CC = lambda Ein, Eout: (jeff_diff_xs(Ein,Eout,diff_nu_p_CC)+jeff_diff_xs(Ein,Eout,diff_nu_n_CC))/2
            self.dsdy_spline_CC_lowe = self.dsdy_spline_CC
            self.dsdy_spline_NC = lambda Ein, Eout: (jeff_diff_xs(Ein,Eout,diff_nu_p_NC)+jeff_diff_xs(Ein,Eout,diff_nu_n_NC))/2
            self.dsdy_spline_NC_lowe = self.dsdy_spline_NC

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
        if self.model=='CSMS':
            if interaction=='NC':
                return self.f_NC(enu)
            elif interaction=='CC':
                return self.f_CC(enu)
        elif self.model=='dipole':
            if(interaction == 'NC'):
                return((10**self.f_NC(np.log10(enu/1e9)))*(units.cm)**2)
            else:
                return((10**self.f_CC(np.log10(enu/1e9)))*(units.cm)**2)

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
        eout = np.atleast_1d(eout)
        diff = np.zeros_like(eout)
        if np.log10(ein) < 0:
            return diff
        if self.model=='dipole':
            if(interaction=='CC'):
                greater_out_msk = (np.log10(eout) >= np.log10(ein))
                diff[greater_out_msk] = 0.
                lowe_msk = (np.log10(eout) < 3.) & ~greater_out_msk
                if np.count_nonzero(lowe_msk) > 0:
                    diff[lowe_msk] = 10**self.dsdy_spline_CC_lowe(np.log10(ein), 
                        np.log10(eout[lowe_msk])
                        )[0][:]/ein
                other_msk = ~(greater_out_msk | lowe_msk)
                if np.count_nonzero(other_msk) > 0:
                    diff[other_msk] = 10**self.dsdy_spline_CC(np.log10(ein), 
                        np.log10(eout[other_msk])
                        )[0][:]/ein
            elif(interaction=='NC'):
                greater_out_msk = np.log10(eout) >= np.log10(ein)
                diff[greater_out_msk] = 0.
                lowe_msk = (np.log10(eout) <= 3.) & ~greater_out_msk
                if np.count_nonzero(lowe_msk) > 0:
                    diff[lowe_msk] = 10**self.dsdy_spline_NC_lowe(np.log10(ein), 
                        np.log10(eout[lowe_msk])
                        )[0][:]/ein
                other_msk = ~(greater_out_msk | lowe_msk)
                if np.count_nonzero(other_msk) > 0:
                    diff[other_msk] = 10**self.dsdy_spline_NC(np.log10(ein), 
                        np.log10(eout[other_msk])
                        )[0][:]/ein
        elif(self.model=='CSMS'):
            if interaction=='NC':
                diff = self.dsdy_spline_NC(ein, eout)
            elif interaction=='CC':
                diff = self.dsdy_spline_CC(ein, eout)
            else:
                raise ValueError('Interaction %s not allowed. Must be "CC" or "NC"')
#            diff = dis.SingleDifferentialCrossSection(ein*units.GeV, eout*units.GeV, 
#                        nsq.NeutrinoCrossSections_NeutrinoFlavor.tau, 
#                        nsq.NeutrinoCrossSections_NeutrinoType.neutrino,
#                        getattr(nsq.NeutrinoCrossSections_Current, interaction) 
#                       )
        return diff
