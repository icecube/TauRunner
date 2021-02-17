import numpy as np
from scipy.integrate import quad
from physicsconstants import PhysicsConstants
units = PhysicsConstants()

class Track(object):

    def __init__(self, depth=0.0):
        self.depth = depth
      
    def x_to_d(self, x):
        '''
        Convert from track parameter to distance traveled
        params
        ______
        x (float): Affine parameter between 0 and 1 which parametrizes the track
    
        returns
        _______
        d (float) Distance traveled along the track in units of radius of the body.
        '''
        pass

    def d_to_x(self, d) :
        '''
        Convert from distance traveled to track parameter
        params
        ______
        d (float) Distance traveled along the track in units of radius of the body.
    
        returns
        _______
        x (float): Affine parameter between 0 and 1 which parametrizes the track
        '''
        pass 

    def r_to_x(self, r):
        '''
        Convert from radiu to track parameter
        params
        ______
        r (float) Radius. Must be between 0 and 1.
    
        returns
        _______
        x (float): Affine parameter between 0 and 1 which parametrizes the track
        '''
        pass

    def x_to_r(self, x):
        '''
        Convert from radius to track parameter
        params
        ______
        x (float): Affine parameter between 0 and 1 which parametrizes the track
    
        returns
        _______
        r (float) Radius. Must be between 0 and 1.
        '''
        pass

    def x_to_r_prime(self, x):
        '''
        Derivative of self.x_to_r with respect to the affine track parameter
        params
        ______
        x (float): Affine parameter between 0 and 1 which parametrizes the track
    
        returns
        _______
        dr/dx (float) Derivative of self.x_to_r with respect to the affine track parameter
        '''
        pass

    def column_depth(self, body, xi, xf, safe_mode=True):
        '''
        params
        ______
        body (Body)  : TauRunner Body object for which you want the column depth
        xi   (float) : Affine track parameter at which to start the integeration.
        xf   (float) : Affine track parameter at which to end the integeration.

        returns
        _______
        column_depth (float) : column depth on the portion of the track from xi to xf [natural units]
        '''
        if not (xi<=xf):
            raise RuntimeError('xi must be less than or equal to xf')
        integrand = lambda x: body.get_density(self.x_to_r(x))*np.abs(self.x_to_r_prime(x))*body.radius
        # find where the path intersects layer boundaries
        xx        = np.append(sorted(self.r_to_x(body.layer_boundaries)), 1)
        # NaNs means it does not intersect the layer
        xx        = xx[np.where(~np.isnan(xx))[0]] # non-nanize
        # Remove xs before and after the integration limits
        mask      = np.where(np.logical_and(xi<xx, xx<xf))[0]
        xx        = np.hstack([[xi], xx[mask], [xf]])
        II        = []
        for xi, xf in zip(xx[:-1], xx[1:]):
            I = quad(integrand, xi, xf)
            if safe_mode:
               if I[1]/I[0] > 1e-4:
                    raise RuntimeError('Error too large')
            II.append(np.abs(I[0]))
        return np.sum(II)
    
    def total_column_depth(self, body, safe_mode=False):
        '''
        params
        ______
        body (Body) : TauRunner Body object for which you want the column depth

        returns
        _______
        column_depth (float) : total column depth for the entire track [natural units]
        '''
        return self.column_depth(body, 0, 1)
