
class Namespace:
    def __init__(self, **kwargs):
        self.__dict__.update( {key: value for key, value in kwargs.items() if not '_' in key} )
    
    def __contains__(self, key):
        return key in self.__dict__.keys()
