class State:
    def __init__(self, emissions):
        self.emissions = emissions
        
    def __getitem__(self,i):
        return self.emissions[i]
        
    def __len__(self):
        return len(self.emissions)
        
    def __str__(self):
        return self.emissions
        
    def is_match_state(self):
        return self.emissions.count('-') < len(self.emissions) / 2.0
        
    def get_emissions(self):
        return self.emissions