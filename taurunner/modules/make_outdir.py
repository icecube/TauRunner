from datetime import date
import os
today = date.today()
TODAYSTR = today.strftime("%Y%m%d")
def make_dirname(savedir, i):
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
