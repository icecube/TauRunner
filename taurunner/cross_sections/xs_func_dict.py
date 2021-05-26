with open('wrkfh.pkl', 'rb') as pkl_f:
    f1 = pkl.load(f)

tot_dict = {
            ('NC',):lambda e: (f1(e)+f2(e)+f3(e)+f4(e)/4.),
            ('CC',),
            ('NC', 'p'),
            ('CC', 'p'),
            ('NC', 'n'),
            ('CC', 'n'),
            ('NC', 'p', 'nu'),
            ('CC', 'p', 'nu'),
            ('NC', 'n', 'nu'),
            ('CC', 'n', 'nu'),
            ('NC', 'p', 'nubar'),
            ('CC', 'p', 'nubar'),
            ('NC', 'n', 'nubar'),
            ('CC', 'n', 'nubar'),
            ('my_func',):their_func,

           }

class XS(object):
    
    def __init__(self, xs_model, iso=True):
        if xs_model not in ['dipole', 'CSMS']:
            raise ValueError('Cross section model %s not recognized' % str(xs_model))
        self.xs_model    = xs_model
        self.spline_dict = {}
        self.key    
        
    def total_xs(self, e, key,):
        return tot_dict[(adkjfhdsf)](e)
