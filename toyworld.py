from reaction import Reaction


class ToyWorld:

    def __init__(self, reactor, product_selection):
        self.product_selection = product_selection
        self.reactor = reactor

    def run(self, generations):

        for i in range(generations):
            reactants = self.reactor.get_reactants()
            products = self.product_selection(reactants)
            self.reactor.react(Reaction(reactants=reactants, products=products))

            print(self.reactor)
