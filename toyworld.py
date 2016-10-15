class ToyWorld:

    def __init__(self, reactor, chemistry, product_selection):
        self.product_selection = product_selection
        self.chemistry = chemistry
        self.reactor = reactor

    def run(self, generations, state):

        for i in range(generations):
            reactants = self.reactor.get_reactants()
            reactions = self.chemistry.enumerate(reactants)
            if len(reactions) == 0:
                break

            reaction = self.product_selection(reactions)
            self.reactor.react(reaction)

            state.add(reaction.as_dict())

        return state
