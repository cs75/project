from collections import deque
from time import sleep


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

        letters = ['A', 'C', 'G', 'T']

        #for each key in the current dictionary
        for key in self.state:

            #for i in the length of the number of states 
            for i in range(len(self.state[key][letters[0]])):

                total_count = float(0)        

                #for each character in the amino acid alphabet 
                for letter in letters:

                    if self.state[key][letter][i] != None:
                        total_count += self.state[key][letter][i]
                        # print total_count

                for letter in letters:
                    if self.state[key][letter][i] != None:
                        self.state[key][letter][i] /= total_count
                    # print "after normalizing: ", self.state[key][letter][i]

            
            init_trans = self.trans.keys()[0]


            for i in range(len(self.trans[init_trans])):

                total_count = float(0)
                
                for key in self.trans:
                    if self.trans[key][i] != None:
                        total_count += self.trans[key][i]

                print "total: ", total_count

                for key in self.trans:
                    if self.trans[key][i] != None:
                        self.trans[key][i] /= total_count
                        print "normalized: ", self.trans[key][i]

                # sleep(1)
                
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
