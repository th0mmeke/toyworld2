class State(object):

    def __init__(self, persistence):
        self.persistence = persistence

    def add(self, state_entry):
        return self.persistence(state_entry)
