from taurunner import track

def make_track(theta, depth=0, track_type='chord'):
    track_style = getattr(track, track_type.lower())
    my_track       = track_style(theta=theta, depth=depth)
    return my_track
