class HMMModel:
    def __init__(self,alphabet,match_state_count,psuedocount):
        self.alphabet = [letter for letter in alphabet]
        width = match_state_count + 1 # accounts for position 0
        every_state = ['mat','ins']
        every_trans = [('mat','mat'),('mat','ins'),('mat','del'),('ins','mat'),('ins','ins'),('ins','del'),('del','mat'),('del','ins'),('del','del')]
        
        self.state = dict.fromkeys(set(every_state))
        for key in every_state:
            self.state[key] = dict.fromkeys(self.alphabet)
        for key in every_state:
            for letter in self.alphabet:
                if key == 'mat':
                    self.state[key][letter] = [None] + [float(psuedocount) for i in range(width-1)]
                else:
                    self.state[key][letter] = [float(psuedocount) for i in range(width)]
        
        self.trans = dict.fromkeys(set(every_trans))
        for key in every_trans:
            if key[0] == 'del':
                self.trans[key] = [None] + [float(psuedocount) for i in range(width-1)]
            else:
                self.trans[key] = [float(psuedocount) for i in range(width)]
        
    def __getitem__(self,item):
        if self.state.has_key(item): 
            return self.state[item]
        elif self.trans.has_key(item): 
            return self.trans[item]
        return None
        
    def normalize(self):
        for key in self.state:
            for letter in self.state[key]:
                total_count = float(0)
                
                for i in range(len(self.state[key][letter])):
                    if self.state[key][letter][i] != None:
                        total_count += self.state[key][letter][i]
                for i in range(len(self.state[key][letter])):
                    if self.state[key][letter][i] != None:
                        self.state[key][letter][i] /= total_count
        
        for key in self.trans:
            total_count = float(0)
            
            for i in range(len(self.trans[key])):
                if self.trans[key][i] != None:
                    total_count += self.trans[key][i]
            for i in range(len(self.trans[key])):
                if self.trans[key][i] != None:
                    self.trans[key][i] /= total_count
                
    def print_model(self):
        print 'MAT Emissions'
        for letter in self.state['mat']:
            print letter, ':      ', self.state['mat'][letter]
        print
        print 'INS Emissions:'
        for letter in self.state['ins']:
            print letter, ':      ', self.state['ins'][letter]
        print
        print 'STATE Transitions'
        for s1 in ['mat','ins','del']:
            print '-----'
            for s2 in ['mat','del','ins']:
                print '(%s,%s)' % (s1,s2), self.trans[(s1,s2)]
        print '\n---------------------------------\n'

M = HMMModel('ACGT',3,0)

M['mat']['A'][1] += 4
M['mat']['C'][3] += 4
M['mat']['G'][2] += 3

M['ins']['A'][2] += 6
M['ins']['G'][2] += 1

M[('mat','mat')][0] += 4
M[('mat','mat')][1] += 3
M[('mat','mat')][2] += 2
M[('mat','mat')][3] += 4

M[('mat','del')][0] += 1
M[('mat','del')][1] += 1

M[('mat','ins')][2] += 1

M[('ins','mat')][2] += 2
M[('ins','del')][2] += 1
M[('ins','ins')][2] += 4

M[('del','mat')][3] += 1
M[('del','del')][1] += 1
M[('del','ins')][2] += 2

M.print_model()

def viterbi(model):
    pass
    
viterbi(M)