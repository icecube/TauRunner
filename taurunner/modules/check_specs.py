def check_specs(TR_specs):
    
    if TR_specs['nevents'] <= 0:
        raise RuntimeError('You must specify more than 0 events to simulate (-n)') 
    # Clean up args so that there is angular arg and energy args
    if (not gzk and theta is None) or (energy is None and spectrum is None):
        raise RuntimeError('You must either pick an energy and theta, use a spectrum, or use the GZK flux')
    
