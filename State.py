class State:
    def __init__(self, position, label, emissions, prev_state, prev_mat_state, prev_ins_state):
        self.position = position
        self.label = label
        self.emissions = emissions
        #the general case 
        self.prev_state = prev_state
        #the specific states 
        self.prev_mat_state = prev_mat_state
        self.prev_ins_state = prev_ins_state
        
    def __getitem__(self,i):
        return self.emissions[i]
        
    def __len__(self):
        return len(self.emissions)
        
    def __str__(self):
        if self.label == 'mat' or self.label == 'ins' or self.label == 'beg' or self.label == 'end':
            return self.label
        return 'Invalid'
        
    def get_emissions(self):
        return self.emissions