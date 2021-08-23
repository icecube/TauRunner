from taurunner import track

def make_tracks(thetas, depth=0, track_type='chord'):
    track_style = getattr(track, track_type.lower())
    # Premake all necessary tracks in case of redundancies
    tracks      = {theta : track_style(theta=theta, depth=depth) for theta in set(thetas)}
    return tracks
