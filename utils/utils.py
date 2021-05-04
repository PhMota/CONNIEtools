
class Namespace:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
    
    def __contains__(self, key):
        return key in self.__dict__.keys()
