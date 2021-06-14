import os, json
from .make_outdir import todaystr, make_outdir
def setup_outdir(args):
    if args.seed is None:
        seed = int(float(savedir.split('/')[-1].replace('_', ''))) % 2**32
    else:
        seed = args.seed
    if args.save is not None:
        d = vars(args)
        d['seed'] = args.seed
        savedir = os.path.join(args.save, '')
        if not os.path.isdir(savedir):
            raise RuntimeError("Directory to save output is not a valid directory")
        savedir = make_outdir(savedir, todaystr)
        os.mkdir(savedir)
        params_file = savedir+"/params.json"
        output_file = savedir+'/output.npy'
        d = vars(args)
        # Check this
        d['seed'] = args.seed
        j = json.dumps(d)
        f = open(params_file,"w")
        f.write(j)
        f.close()
    else:
        savedir=None
        params_file=None
        output_file=None
    if args.seed is None:
        seed = int(float(savedir.split('/')[-1].replace('_', ''))) % 2**32
    else:
        seed = args.seed
    
    return seed, savedir, params_file, output_file
