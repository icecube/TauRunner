class Callable():

    def __init__(self, obj, scale=1.0):
        self.obj = obj
        self.scale = scale
        if hasattr(obj, '__call__'):
            self.iscallable = True
        else:
            self.iscallable = False

    def __mul__(self, other):
        return Callable(self.obj, scale=other)

    def __rmul__(self, other):
        return Callable(self.obj, scale=other)

    def __call__(self, x):
        if self.iscallable:
            return self.obj(x)*self.scale
        else:
            return self.obj*self.scale
