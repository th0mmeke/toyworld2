from collections import Counter
import logging


class ToyWorld2:

    def __init__(self, reactor, chemistry, product_selection):
        self.product_selection = product_selection
        self.chemistry = chemistry
        self.reactor = reactor

    def run(self, generations, state):

        for i in range(generations):
            partial_reaction = self.reactor.get_reactants()
            reactions = self.chemistry.enumerate(partial_reaction)
            reaction = self.product_selection(reactions)

            if reaction is not None:
                self.reactor.react(reaction)
                logging.info("Reaction between " + str([r.get_symbol() for r in reaction.reactants]) + " giving " + str([p.get_symbol() for p in reaction.products]))
                state.add(reaction.as_dict())

        return state
