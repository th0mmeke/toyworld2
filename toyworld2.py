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
        e = environment.pop()
        while generation <= generations:
            try:
                partial_reaction = self.reactor.get_reactants()
            except ValueError:
                logging.info("Stopping, lack of energy")
                break
            reactions = self.chemistry.enumerate(partial_reaction)
            reaction = weighted_selection(reactions, self.product_selection, e)

            if reaction is None:
                non_reaction += 1
                if non_reaction > 200:
                    logging.info("Stopping, lack of reactions")
                    break
            else:
                self.reactor.react(reaction)
                logging.info("{}: bias={}, reaction between {} giving {}".format(generation, e, str([r.get_symbol() for r in reaction.reactants]),
                                                                        str([p.get_symbol() for p in reaction.products])))
                state.add(reaction.as_dict())
                # self.reactor.update_environment(environment.pop())
                generation += 1
                non_reaction = 0
                e = environment.pop()

        return state
