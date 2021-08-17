def is_floatable(element):
    try:
        float(element)
        return True
    except:
        return False
