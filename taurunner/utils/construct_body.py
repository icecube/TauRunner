def construct_body(TR_specs):
    if TR_specs['body']=='earth':
        from taurunner.body import lumen_sit
        if TR_specs['water']>0:
            body = lumen_sit([(TR_specs['water'], 1)])
        else:
            body = lumen_sit()
    elif 'sun' in TR_specs['body'].lower():
        body = getattr(taurunner.body, TR_specs['body'])
    else:
        try:
            density = float(TR_specs['body'])
            if density<=0:
                raise ValueError('Density must be strictly positive.')
            if TR_specs['radius']<=0:
                raise ValueError('Radius must be strictly positive')
            from taurunner.body import Body
            radius_arg = float(TR_specs['radius'])
            body = Body(density, radius_arg)
        except ValueError as e:
            print('Not known how to handle body arg, %s' % args.body)
            raise e
    return body
