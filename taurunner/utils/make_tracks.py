from taurunner import track

def make_tracks(TR_specs, thetas, radius):
    track_style = getattr(track, TR_specs['track'])
    # Premake all necessary tracks in case of redundancies
    tracks      = {theta : track_style(theta=theta, depth=TR_specs['depth']/radius) for theta in set(thetas)}
    return tracks
