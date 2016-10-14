class ToyWorld:

    def __init__(self, reactor, generate_reaction_set, product_selection):
        self.product_selection = product_selection
        self.generate_reaction_set = generate_reaction_set
        self.reactor = reactor

    def run(self, generations, state):

        for i in range(generations):
            reactants = self.reactor.get_reactants()
            reactions = self.generate_reaction_set(reactants)
            reaction = self.product_selection(reactions)

            self.reactor.react(reaction)

            state.add(reaction.as_dict())

        return state
