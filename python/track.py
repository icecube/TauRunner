class Track():

    def __init__(self, depth=0.0):
        self.depth = depth
        

    def r_to_x(self, r):
        pass

    def x_to_r(self, x):
        pass

    def r_to_x_prime(self, x):
        pass

    def column_depth(self, body, xi, xf, safe_mode=True):
        r'''
        params
        ______
        body (Body) :

        returns
        _______
        column_depth (float) :
        '''
        if not (xi<=xf):
            raise RuntimeError('xi must be less than or equal to xf')
        integrand = lambda x: body.get_density(self.x_to_r(x))*self.x_to_r_prime(x)*body.radius
        # find where the path intersects layer boundaries
        xx        = np.append(sorted(self.r_to_x(body.layer_boundaries)), 1)
        # NaNs means it does not intersect the layer
        xx        = xx[np.where(~np.isnan(xx))[0]] # non-nanize
        # Remove xs before and after the integrastion limits
        mask      = np.where(np.logical_and(xi<xx, xx<xf))[0]
        xx        = np.hstack([[xi], xx[mask], [xf]])
        II        = []
        for xi, xf in zip(xx[:-1], xx[1:]):
            I = quad(integrand, xi, xf)
            if safe_mode:
               if I[1]/I[0] > 1e-4:
                    raise RuntimeError('Error too large')
            II.append(np.abs(I[0]))
        return np.sum(II)/units.gr*units.cm**2
    
    def total_column_depth(self, body, safe_mode=False):
        r'''
        params
        ______
        track (Track) :

        returns
        _______
        column_depth (float) :
        '''
        return self.column_depth(body, 0, 1)
