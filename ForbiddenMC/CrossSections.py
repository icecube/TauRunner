from Casino import *
units = nsq.Const()
dis = nsq.NeutrinoDISCrossSectionsFromTables()
tds = nsq.TauDecaySpectra()


class CrossSections(object):
    def __init__(self, cross_section_path, seed):
        #cross sections from 1e8 - 1e16 GeV patched with nuSQUIDS 
	#temporary until EHE cross sections are added to nuSQUIDS
        
	self.rand = np.random.RandomState(seed=seed)

	######################################
	charged = np.load(cross_section_path + 'nuXS_CC_8-16.npy')
	neutral = np.load(cross_section_path + 'nuXS_NC_8-16.npy')
	energies = np.logspace(17, 25, 500)
	energies = (energies[:-1] + energies[1:])/2

	self.log_e = np.log10(energies)
	self.log_XS_CC = np.log10(charged)
	self.log_XS_NC = np.log10(neutral)

	self.f_NC = interp1d(self.log_e, self.log_XS_NC)
	self.f_CC = interp1d(self.log_e, self.log_XS_CC)

	self.dsdy_spline_CC = np.load(cross_section_path + 'dsigma_dy_CC.npy').item()
	self.dsdy_spline_CC_lowe = np.load(cross_section_path + 'dsigma_dy_CC_lowE.npy').item()
	self.dsdy_spline_NC = np.load(cross_section_path + 'dsigma_dy_NC.npy').item()
	self.dsdy_spline_NC_lowe = np.load(cross_section_path + 'dsigma_dy_NC_lowE.npy').item()

	#####################################

    def TotalNeutrinoCrossSection(self, enu,
                              flavor = nsq.NeutrinoCrossSections_NeutrinoFlavor.tau,
                              neutype = nsq.NeutrinoCrossSections_NeutrinoType.neutrino,
                              interaction = nsq.NeutrinoCrossSections_Current.NC):
        r'''
        Calculates total neutrino cross section. returns the value of sigma_CC (or NC)
	in natural units.
        Parameters
        ----------
        enu:         float
	    neutrino energy in eV
        flavor:      nusquids obj
	    nusquids object defining neutrino flavor. default is tau
        interaction: nusquids obj
	    nusquids object defining the interaction type (CC or NC). default is NC
	Returns
        -------
        TotalCrossSection: float
	    Total neutrino cross section at the given energy in natural units.
        '''

	if(np.log10(enu) > self.log_e[0]):

            if(interaction == nsq.NeutrinoCrossSections_Current.NC):
                return((10**self.f_NC(np.log10(enu)))*(units.cm)**2)
            else:
                return((10**self.f_CC(np.log10(enu)))*(units.cm)**2)
	else:

	    return dis.TotalCrossSection(enu,flavor,neutype,interaction)*(units.cm)**2

    def DifferentialOutGoingLeptonDistribution(self, ein, eout, interaction):
	if(np.log10(ein) < 0):
	    diff = 0
	    return diff
        if(interaction==nsq.NeutrinoCrossSections_Current.CC):
            if (np.log10(eout) >= np.log10(ein)):
                diff = 0.
		return diff
            elif(np.log10(eout) < 3.):
		diff = 10**self.dsdy_spline_CC_lowe(np.log10(ein), np.log10(eout))[0][0]/ein 
	    else:
                diff = 10**self.dsdy_spline_CC(np.log10(ein), np.log10(eout))[0][0]/ein
        elif(interaction==nsq.NeutrinoCrossSections_Current.NC):
            if (np.log10(eout) >= np.log10(ein)):
                diff = 0.
		return diff
            elif( np.log10(eout) <= 3.):
	        diff = 10**self.dsdy_spline_NC_lowe(np.log10(ein), np.log10(eout))[0][0]/ein
	    else:
                diff = 10**self.dsdy_spline_NC(np.log10(ein), np.log10(eout))[0][0]/ein

        return diff








