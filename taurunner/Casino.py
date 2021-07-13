import os
os.environ['HDF5_DISABLE_VERSION_CHECK']='2'
import numpy as np
import subprocess
from scipy.interpolate import InterpolatedUnivariateSpline as iuvs
from taurunner.modules import units

TOL  = 0.0

# Auxiliary function for secondary antinu sampling

def chunks(lst, n): # pragma: no cover
    for i in range(0, len(lst), n):
        yield lst[i:i+n]

def DoAllCCThings(objects, xs, losses=True):
    r'''
    Calling MMC requires overhead, so handle all MMC calls per iterations
    over the injected events at once
    Parameters
    ----------
    objects: list
        List of CasinoEvents that need to have tau losses sampled stochastically.
    xs: str 
        Cross section model to use for the photohadronic losses
    losses: bool
        This can be set to False to turn off energy losses. In this case, the particle decays at rest.
    Returns
    -------
    objects: list
        List of CasinoEvents after losses are calculated
    '''
    final_values= []
    efinal, distance = [], []
    pid = int(objects[0][-1])
    if pid in [13,14]:
      flavor='mu'
    elif pid in [15,16]: # pragma: no cover
      flavor='tau'
    e      = [obj[0]/units.GeV for obj in objects]                    #MMC takes initial energy in GeV 
    dists  = [1e3*(obj[6] - obj[2])/units.km for obj in objects]      #distance to propagate is total distance minus the current position in m 
    mult   = [obj[-2]*(units.cm**3)/units.gr/2.7 for obj in objects]  #convert density back to normal (not natural) units
    sort         = sorted(list(zip(mult, e, dists, objects)))
    sorted_mult  = np.asarray(list(zip(*sort))[0])
    sorted_e     = np.asarray(list(zip(*sort))[1])
    sorted_dists = np.asarray(list(zip(*sort))[2])
    sorted_obj   = np.asarray(list(zip(*sort))[3])

    if(not losses):
        final_energies = sorted_e
        final_distances = np.zeros(len(sorted_e))
        for i, obj in enumerate(sorted_obj):
            obj[0] = final_energies[i]*units.GeV
            obj[5] = final_distances[i]
        return(sorted_obj)
    else: # pragma: no cover
        split = np.append(np.append([-1], np.where(sorted_mult[:-1] != sorted_mult[1:])[0]), len(sorted_mult))
        propagate_path = os.path.dirname(os.path.realpath(__file__))
        if(xs.model=='dipole'):
            propagate_path+='/propagate_{}s.sh'.format(flavor)
        elif(xs.model=='CSMS'):
            propagate_path+='/propagate_{}s_ALLM.sh'.format(flavor)
        else:
            raise ValueError("Cross section model error. Cross section model is %s" % xs.model)
        for i in range(len(split)-1):
            multis = sorted_mult[split[i]+1:split[i+1]+1]
            eni = sorted_e[split[i]+1:split[i+1]+1]
            din = sorted_dists[split[i]+1:split[i+1]+1]
            max_arg = 500
            eni_str = ['{} {}'.format(e, d) for e,d in list(zip(eni, din))]
            eni_str = list(chunks(eni_str, max_arg))
            for kk in range(len(eni_str)):
                eni_str[kk].append(str(multis[0]))
                eni_str[kk].insert(0, propagate_path)
                process = subprocess.check_output(eni_str[kk])
                for line in process.split(b'\n')[:-1]:
                    final_values.append(float(line.replace(b'\n',b'')))

        final_energies = np.asarray(final_values)[::2]
        final_distances = np.abs(np.asarray(final_values)[1::2])/1e3
        final_distances += final_distances*0.05
        for i, obj in enumerate(sorted_obj):
            obj[0] = final_energies[i]*units.GeV/1e3
            obj[5] = final_distances[i]
        return(sorted_obj)



##########################################################
##########################################################


#This is the propagation algorithm. The MCmeat, if you will.
def Propagate(particle, track, body):
    total_column_depth = track.total_column_depth(body)
    #keep iterating until final column depth is reached or a charged lepton is made
    while(not np.any((particle.position >= 1.) or (particle.isCC))):
        if(particle.ID in [12, 14, 16]):
            #Determine how far you're going
            p1 = particle.rand.random_sample()
            DepthStep = particle.GetProposedDepthStep(p1)
            CurrentDepth=track.x_to_X(body, particle.position)
            if(CurrentDepth+DepthStep >= total_column_depth):
                particle.position=1.
                return particle
            else:
                particle.position=track.X_to_x(body, CurrentDepth+DepthStep)
            #now pick an interaction
            p2 = particle.rand.random_sample()
            CC_lint = particle.GetInteractionDepth(interaction='CC')
            p_int_CC = particle.GetTotalInteractionDepth() / CC_lint
            if(p2 <= p_int_CC):
                particle.Interact('CC')
            else:
                particle.Interact('NC')
            if(particle.isCC):
                continue
        elif(np.logical_or(particle.ID == 15, particle.ID == 13)):
            current_distance=track.x_to_d(particle.position)
            charged_distance = particle.chargedposition*units.km/body.radius
            if(track.d_to_x(current_distance+charged_distance) >=1.): # pragma: no cover
                particle.position=1.
            else:
                current_distance+=charged_distance
                particle.position=track.d_to_x(current_distance)
            if(particle.position >= 1-TOL): # pragma: no cover
                return particle
                continue
            else:
                particle.Decay()
    return particle
