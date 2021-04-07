from datetime import date
import os
today = date.today()
todaystr = today.strftime("%Y%m%d")
def make_dirname(savedir, todaystr, i):
    stri     = ('000000'+str(i))[-6:]
    dir_name = '%s/%s_%s' % (savedir, todaystr, stri)
    return dir_name

def make_outdir(savedir, todaystr):
    i            = 0
    proposed_dir = make_dirname(savedir, todaystr, i)
    while os.path.isdir(proposed_dir):
        i+=1
        proposed_dir = make_dirname(savedir, todaystr, i)
    return proposed_dir
