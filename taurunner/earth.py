import numpy as np

#########################################################
#####  Earth Model Handling. Here we define PREM ########
####   params and handle ice layers/detector     ########
#########################################################

prem_params = {0: [13.0885, 0, -8.8381, 0], 1: [12.5815, -1.2638, -3.6426, -5.5281], 
          2: [7.9565, -6.4761, 5.5283, -3.0807], 3: [5.3197, -1.4836, 0, 0], 
          4: [11.2494, -8.0298, 0, 0], 5: [7.1089, -3.8045, 0, 0],
          6: [2.6910, 0.6924, 0, 0], 7: [2.9, 0, 0, 0], 8:[2.6, 0, 0, 0],
          9: [0.92, 0, 0, 0]}

def rho_earth(theta, x, d = 0, water_layer = 0):
    """ Returns the Earth density in gr/cm^3.

    Args:
        theta: zenith angle in radians.
        x: position along the trajectory in km.
        d: depth under the Earth's surface. (default set to 0)

    Returns:
        rho: density in gr/cm^3
    """
    r = np.sqrt((REarth - d)**2 + x**2 + 2. * (REarth - d) * x * np.cos(theta))

    if r <= 1221.:
        param = prem_params[0]
    elif r <= 3480:
        param = prem_params[1]
    elif r <= 5701:
        param = prem_params[2]
    elif r <= 5771:
        param = prem_params[3]
    elif r <= 5971:
        param = prem_params[4]
    elif r <= 6151:
        param = prem_params[5]
    elif r <= 6346.6:
        param = prem_params[6]
    elif r <= 6356:
        param = prem_params[7]
    elif r <= 6368:
        param = prem_params[8]
    elif(r <= REarth):
        param = prem_params[9]
    else:
        raise RuntimeError("sampling densities outside earth boundaries")

    rho = param[0] + (param[1] * r/REarth) + (param[2] * (r/REarth)**2) + (param[3] * (r/REarth)**3)
    return rho

def get_params(layers, aradius):
    layer = np.min(np.where(layers > aradius))
    return prem_params[layer]

def get_average(param, a, b, REarth):
    piece1 = (a/REarth)*param[0] + (param[1]*(a/REarth)**2)/2. + (param[2]*(a/REarth)**3)/3. + (param[3]*(a/REarth)**4)/4.
    piece2 = (b/REarth)*param[0] + (param[1]*(b/REarth)**2)/2. + (param[2]*(b/REarth)**3)/3. + (param[3]*(b/REarth)**4)/4.
    piece3 = REarth/(b - a)
    return(piece3*(piece2 - piece1))

def get_radii_densities(water_layer):
    REarth = 6368.
    
    layers = [1221., 3480., 5701., 5771., 5971., 6151., 6346.6, 6356., 6368., REarth]
    jump_depths = REarth - np.asarray(layers)
    jump_depths = np.sort(np.append(jump_depths[:-1], np.linspace(REarth - layers[2], REarth, 13)[1:]))
    avg_rho = []
    jump_radius = REarth - jump_depths
    jump_radius = np.append(jump_radius, REarth)
    for left_jump, right_jump in zip(jump_radius[:-1], jump_radius[1:]):
        aradius = (left_jump + right_jump)/2.
        param = get_params(layers, aradius)
        avg_rho.append(get_average(param, left_jump, right_jump, REarth))

    radii, densities = np.flip(jump_radius, axis=0)[1:], np.flip(avg_rho, axis=0)
    if(water_layer > 0):
        radii = np.append(radii, REarth+water_layer)
        densities = np.append(densities, 0.92)
    return radii, densities

def GetDistancesPerSection(theta, radii, water_layer=0):
        r'''
        Calculates distance traveled in each uniform density region,
        given in order with index of density region
        Parameters
        ----------
        radii: list
            radii of earth regions of different uniform density
        Returns
        -------
        dists_by_sect: list
            distance traveled in each section of uniform density,
            along with the index of the region in an additional dimension
        '''

        closest_approach = radii[-1] * np.sin(theta)
        smallest_sect = np.min(np.where(closest_approach < radii))
        dists = []
        for i in range(smallest_sect, len(radii)):
            side_length = np.sqrt(radii[i]**2 - closest_approach**2)
            dists.append(side_length - np.sum(dists))
        dists = np.flip(dists, axis=0)
        dists[-1] = 2 * dists[-1]
        dists = np.append(dists, np.flip(dists[:-1], axis=0))
        sects = list(range(len(radii) - 1, smallest_sect, -1)) + \
                    list(range(smallest_sect, len(radii)))
        
        return dists, sects

###########################################
##### Earth Model handling ends here ######
###########################################
