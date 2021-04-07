import os
def cleanup_outdir(savedir, output_file, params_file):
    if os.path.exists(output_file):
        os.remove(output_file)
    if os.path.exists(params_file):
        os.remove(params_file)
    if os.path.isdir(savedir):
        os.rmdir(savedir)
