from abc import ABCMeta

class SingletonMetaclass(ABCMeta):
    """
    SingletonMetaclass is a metaclass that ensures a class only has one instance.

    This metaclass overrides the __call__ method to check if an instance of the class
    already exists. If it does, it returns the existing instance. If not, it creates
    a new instance, stores it in the _instances dictionary, and then returns it.

    Attributes:
        _instances (dict): A dictionary that stores the single instances of classes
                           that use this metaclass.

    Methods:
        __call__(cls, *args, **kwargs): Checks if an instance of the class exists.
                                        If not, creates and stores a new instance.
                                        Returns the single instance of the class.
    """
    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            instance = super().__call__(*args, **kwargs)
            cls._instances[cls] = instance
        return cls._instances[cls]