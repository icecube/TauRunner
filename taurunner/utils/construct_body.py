def construct_body(body_name, water_layer, radius=None):
    if body_name=='earth':
        from taurunner.body import construct_earth
        if water_layer > 0:
            body = construct_earth([water_layer, 1])
        else:
            body = construct_earth()

    elif 'sun' in body_name.lower(): # TODO make this less stupid
        from taurunner.body import construct_sun
        body = construct_sun(body_name)
    else:
        density = float(body_name)
        if density <= 0:
            raise ValueError('Density must be strictly positive.')
        if radius is None or radius<=0:
            raise ValueError('Radius must be strictly positive')
        from taurunner.body import Body
        body = Body(density, radius)
    return body
