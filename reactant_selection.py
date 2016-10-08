from abc import ABCMeta, abstractmethod


class IReactantSelection(object):
    __metaclass__ = ABCMeta

    '''
    Abstract class (interface), defining the interface for reactant selection.
    '''

    @abstractmethod
    def get_reactants(self):
        '''

        :return: list of Element
        '''

        pass

    @abstractmethod
    def react(self, reaction):
        pass
