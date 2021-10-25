import warnings
import numpy as np
from scipy.interpolate import RectBivariateSpline, splev, splrep
from scipy.integrate import quad
from importlib.resources import path as imppath
import os
import pickle as pkl

from taurunner.utils import FileLock
import taurunner as tr


with imppath('taurunner.resources.column_depth_splines', '__init__.py') as p:
    SPLINE_PATH = str(p).split('__init__.py')[0]


def column_depth(track,
                 body,
                 xi: float,
                 xf: float,
                 safe_mode=True
    ) -> float:
    '''
    params
    ______
    body      : TauRunner Body object for which you want the column depth
    xi        : Affine track parameter at which to start the integeration.
    xf        : Affine track parameter at which to end the integeration.
    safe_mode : If True make sure the error on the integral is small

    returns
    _______
    column_depth : column depth on the portion of the track from xi to xf [natural units]
    '''
    if not (xi<=xf):
        raise RuntimeError('xi must be less than or equal to xf')
    integrand = lambda x: body.get_density(track.x_to_r(x))*track.x_to_d_prime(x)*body.radius
    # find where the path intersects layer boundaries
    xx        = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for r in body.layer_boundaries:
            xx = np.append(xx, track.r_to_x(r))
    # NaNs means it does not intersect the layer
    xx        = xx[np.where(~np.isnan(xx))[0]] # non-nanize
    # Remove xs before and after the integration limits
    mask      = np.where(np.logical_and(xi<xx, xx<xf))[0]
    xx        = np.hstack([[xi], xx[mask], [xf]])
    II        = []
    for xi, xf in zip(xx[:-1], xx[1:]):
        I = quad(integrand, xi, xf, full_output=1)
        if safe_mode:
            if I[0]!=0:
                if I[1]/I[0]>1e-3:
                    raise RuntimeError('Error too large')
        II.append(I[0])
    return np.sum(II)

def round_sig(x, sig=3):
    if x==0:
        return 0.
    else:
        return round(x, sig-int(np.floor(np.log10(abs(x))))-1)


# Precompute column depths when the object is initialized so that we are computing integrals
# every time
def column_depth_helper(track, body):
    xx            = np.linspace(0, 1, 301)
    column_depths = np.append(0, np.cumsum([column_depth(track, body, x[0], x[1]) for x in zip(xx[:-1], xx[1:])]))
    # Pad arrays to help spline stability
    npad          = 5
    xpad          = np.linspace(0.01, 0.05, npad)
    padded_xx     = np.hstack([-xpad[::-1], xx, 1+xpad])
    padded_cds    = np.hstack([np.full(npad, column_depths[0]), column_depths, np.full(npad, column_depths[-1])])
    # Make the splines
    x_to_X_tck    = splrep(padded_xx, padded_cds)
    x_to_X        = lambda x: splev(x, x_to_X_tck)
    
    return column_depths[-1], x_to_X

def get_hash(track, body):
    s = f'{track.desc}_{body.radius}{track.depth}'
    for bd in body.layer_boundaries:
        s += str(round_sig(bd, sig=3))
        s += str(round_sig(body.get_density(bd), sig=3))
    s = s.replace('.', '').replace('+', '')
    hash_s = s
    return hash_s

def spline_fname(hash_s):    
    x2X_fname = f'{SPLINE_PATH}/{hash_s}_x_to_X.pkl'
    spline_exists = os.path.isfile(x2X_fname)
    return x2X_fname, spline_exists

def construct_X2x(x2X):
    if x2X(0)==x2X(1): # there is no bijection since derivative of this function must be non-neg
        if x2X(0.5)!=0: # This should only happen in the case that there is no column depth
            raise ValueError('Constant X2x is only allow for tracks with no column depth')
        else:
            def X2x(x):
                raise Exception('You should not be calling this for null propagations')
    else:
        xx         = np.linspace(0, 1, 301)
        cds        = [float(x2X(x)) for x in xx]
        cds[0]     = 0. # Sometimes spline support can be imperfect
        X2x_tck    = splrep(cds, xx)
        def X2x(X):
            return float(splev(X, X2x_tck))
    return X2x

def set_spline(track, body, npad=5):
    hash_s = get_hash(track, body)
    # We do some fancy precomputing for this stuff
    if track.desc=='chord':
        x2X_fname, spline_exists = spline_fname(hash_s)
        if spline_exists:
            with FileLock(x2X_fname):
                with open(x2X_fname, 'rb') as pkl_file:
                    x2X_spline =  pkl.load(pkl_file)
        else:
            from taurunner.track import chord
            if track.depth==0:
                skip = True
            else:
                skip = False
            ths = np.arccos(np.linspace(-1, 1, 500))
            chs = [chord(theta=th, depth=track.depth) for th in ths]
            xx  = np.linspace(0, 1, 301)
            res = np.zeros((ths.shape[0], xx.shape[0]))
            for i, ch in enumerate(chs):
                if skip and np.cos(ch.theta)<=0:
                    continue
                tot_X, x2X = column_depth_helper(ch, body)
                for j, x in enumerate(xx):
                    res[i,j] = x2X(x)
            res[np.where(res<=0)] = np.min(res[np.where(res>0)])
            x2X_spline = RectBivariateSpline(np.cos(ths), xx, res, kx=2, ky=2)
            with FileLock(x2X_fname):
                with open(x2X_fname, 'wb') as pkl_f:
                    pkl.dump(x2X_spline, pkl_f)
    
        # reduce dimensionality of the spline to the desired theta    
        x2X = lambda x: x2X_spline(np.cos(track.theta), x)[0]
    else: # we're not gonna do anything smart
        xx  = np.linspace(0, 1, 301)
        tot_X, x2X = column_depth_helper(track, body)
        
    tot_X = x2X(1)
    
    # Construct the X2x on the fly
    X2x = construct_X2x(x2X)
    
    track._column_depth_lookup[hash_s] = (tot_X, x2X, X2x)
