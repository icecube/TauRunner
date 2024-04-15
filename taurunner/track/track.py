import numpy as np
from taurunner.body import Body
from taurunner.track.utils import get_hash, set_spline

class Track(object):

    def __init__(self, depth=0.0, theta=None):
        self.depth                = depth
        self._column_depth_lookup = {}
        self.desc                 = 'generic_track'
      
    def d_to_x(self, d: float) -> float:
        '''
        Convert from distance traveled to track parameter

        params
        ______
        d : Distance traveled along the track in units of radius of the body.
    
        returns
        _______
        x : Affine parameter between 0 and 1 which parametrizes the track
        '''
        pass 

    def r_to_x(self, r: float) -> float:
        '''
        Convert from radius to track parameter

        params
        ______
        r : Radius. Must be between 0 and 1.
    
        returns
        _______
        x : Affine parameter between 0 and 1 which parametrizes the track
        '''
        pass

    def total_column_depth(self, body: Body, safe_mode=True) -> float:
        '''
        params
        ______
        body      : TauRunner Body object for which you want the column depth
        safe_mode : If True make sure the error on the integral is small

        returns
        _______
        column_depth : total column depth for the entire track [natural units]
        '''
        hash_s = get_hash(self, body)
        if hash_s not in self._column_depth_lookup.keys():
            set_spline(self, body)
        return self._column_depth_lookup[hash_s][0]

    def x_to_cartesian_direction(x):
        '''
        Convert from track parameter to direction in cartesian space

        params
        ______
        x (float): Affine parameter between 0 and 1 which parametrizes the track

        returns
        _______
        cart_dir (tuple): x-, y-, and z-directions in a cartesian system whose origin coincides
                          with the center of the body

        '''
        pass

    def x_to_d(self, x: float) -> float:
        '''
        Convert from track parameter to distance traveled

        params
        ______
        x : Affine parameter between 0 and 1 which parametrizes the track
    
        returns
        _______
        d : Distance traveled along the track in units of radius of the body.
        '''
        pass

    def x_to_d_prime(self, x: float) -> float:
        '''
        Derivative of self.x_to_d with respect to the affine track parameter

        params
        ______
        x : Affine parameter between 0 and 1 which parametrizes the track
    
        returns
        _______
        dd/dx : Derivative of self.x_to_r with respect to the affine track parameter
        '''
        pass

    def x_to_r(self, x:float) -> float:
        '''
        Convert from radius to track parameter

        params
        ______
        x : Affine parameter between 0 and 1 which parametrizes the track
    
        returns
        _______
        r : Radius. Must be between 0 and 1.
        '''
        pass

    def X_to_x(self, body : Body, X: float) -> float:
        '''
        Returns the affine track parameter at which the track has traversed a given column depth
        in a given body. If column depth is greater than the total column depth it raises an error
        
        params
        ______
        body : TauRunner body object
        X    : Column depth [natural units]

        returns
        _______
        x : Affine track parameter
        '''
        hash_s = get_hash(self, body)
        if hash_s not in self._column_depth_lookup.keys():
            set_spline(self, body)
        max_X = self._column_depth_lookup[hash_s][0]
        if X<=max_X:
            return self._column_depth_lookup[hash_s][2](X)
        else:
            raise ValueError('Column depth was greater than total')

    def x_to_X(self, body: Body, x: float) -> float:
        r'''
        Returns the total column depth the particle has crossed after traversing up 
        to a given affine track parameter

        Params
        ______
        body : TauRunner body object
        x    : Affine track parameter

        Returns
        _______
        X : Column depth [natural units]
        '''
        hash_s = get_hash(self, body)
        if hash_s not in self._column_depth_lookup.keys():
            set_spline(self, body)
        return self._column_depth_lookup[hash_s][1](x)

    def x_to_pp_dir(self, x):
        pass

    def x_to_pp_pos(self, x):
        pass
