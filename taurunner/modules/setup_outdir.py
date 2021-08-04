import os, json
from datetime import date

def make_dirname(savedir, i):
    today    = date.today()
    TODAYSTR = today.strftime("%Y%m%d")
    stri     = ('000000'+str(i))[-6:]
    dir_name = '%s/%s_%s' % (savedir, TODAYSTR, stri)
    return dir_name

def make_outdir(savedir):
    i            = 0
    proposed_dir = make_dirname(savedir, i)
    while os.path.isdir(proposed_dir):
        i+=1
        proposed_dir = make_dirname(savedir, i)
    return proposed_dir

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
