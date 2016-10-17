from collections import Counter


class ToyWorld:

    def __init__(self, reactor, chemistry, product_selection):
        self.product_selection = product_selection
        self.chemistry = chemistry
        self.reactor = reactor

    def run(self, generations, state):

        print(Counter([str(x) for x in self.reactor.get_population()]))

        for i in range(generations):
            reactants = self.reactor.get_reactants()
            if len(reactants) == 0:
                break

            reactions = self.chemistry.enumerate(reactants)
            if len(reactions) == 0:
                break

            reaction = self.product_selection(reactions)
            self.reactor.react(reaction)

            state.add(reaction.as_dict())

        print(Counter([str(x) for x in self.reactor.get_population()]))

        return state
