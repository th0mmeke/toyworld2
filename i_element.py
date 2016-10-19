from abc import ABCMeta, abstractmethod


class IElement(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def get_symbol(self):
        pass
