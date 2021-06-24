import os, json
from .make_outdir import make_outdir
def setup_outdir(TR_specs):
    if not TR_specs['prefix']:
        outdir = make_outdir(TR_specs['base_savedir'])
        TR_specs['prefix'] = outdir.split('/')[-1]
    else:
        outdir = '%s/%s' % (TR_specs['base_savedir'], TR_specs['prefix'])
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    if not os.path.isdir(outdir):
        raise ValueError('Specified outdir is a file already')
    if not TR_specs['seed']: # we must have a seed if saving
        seed = hash(outdir) % (2**32) # Make seed from hashing outdir
        TR_specs['seed'] = seed
    return TR_specs
