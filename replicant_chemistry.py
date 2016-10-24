from semi_realistic_chemistry import SemiRealisticChemistry
import distance
import itertools
import sys


class ReplicantChemistry(SemiRealisticChemistry):

    def enumerate(self, partial_reaction):

        """
        Selection pressure towards A+X -> 2A + Y

        Weight reaction options by:

        - 3 or more products
        - One product contains a section of a reactant, the more the better
        - Longer A's are better

        :param partial_reaction:
        :return:
        """

        reaction_options = super(ReplicantChemistry, self).enumerate(partial_reaction=partial_reaction)

        # Now weight options

        for reaction in reaction_options:
            reaction.product_value = 0

            # Only consider reactions with 3 products - this also eliminates unproductive null "reactions" where the products == reactants
            if len(reaction.products) > 2:

                # Similarity between reactants and products
                maximum = 0
                for reactant in reaction.reactants:

                    # Now find the best combination of products - 2 products with the smallest combined distance from this reactant
                    for combination in itertools.combinations(reaction.products, 2):
                        weight = len(reactant.get_symbol())  # start with length of A
                        weight -= sum([distance.levenshtein(reactant.get_symbol(), product.get_symbol()) for product in combination])
                        if weight > maximum:
                            maximum = weight

                reaction.product_value = maximum  # add to weight for multiple products

        return reaction_options
