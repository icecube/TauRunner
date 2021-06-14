def construct_body(body_arg, radius_arg, **kwargs):
    if body_arg in ['earth', 'sun']:
        if(body_arg=='earth'):
            from taurunner.body import create_earth
            body = create_earth(kwargs['layer'], kwargs['density'])
        elif(body_arg=='sun'):
            from taurunner.body import HZ_Sun
            body = HZ_Sun
    else:
        try:
            density = float(body_arg)
            if density<=0:
                raise ValueError('Density must be strictly positive.')
            if radius_arg is None:
                raise ValueError('Must specify a radius if you give a numerical value for body')
            from taurunner.body import Body
            body = Body(density, radius_arg)
        except ValueError as e:
            print('Not known how to handle body arg, %s' % args.body)
            raise e
    return body
