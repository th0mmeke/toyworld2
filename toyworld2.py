import logging
from weighted_selection import weighted_selection


class ToyWorld2:

    def __init__(self, reactor, chemistry, product_selection):
        self.chemistry = chemistry
        self.reactor = reactor
        self.product_selection = product_selection

    def run(self, generations, state, environment):

        generation = 1
        non_reaction = 0
        while generation <= generations:
            partial_reaction = self.reactor.get_reactants()
            reactions = self.chemistry.enumerate(partial_reaction)
            reaction = weighted_selection(reactions, self.product_selection)

            if reaction is None:
                non_reaction += 1
                if non_reaction > 200:
                    logging.info("Stopping, lack of reactions")
                    break
            else:
                self.reactor.react(reaction)
                logging.info("{}: reaction between {} giving {}".format(generation, str([r.get_symbol() for r in reaction.reactants]),
                                                                        str([p.get_symbol() for p in reaction.products])))
                state.add(reaction.as_dict())
                generation += 1
                non_reaction = 0

        return state
