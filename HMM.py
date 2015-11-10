from HMMModel import *
from State import *

def Align(data):
    '''Description
    '''    
    training_data = open(data, 'r')
    aligned_sequences = []
    
    curr_seq = None
    curr_line = training_data.readline()
    
    while curr_line != '': # populate aligned_sequences from training_data
        if curr_line[0] != '>':
            curr_seq += curr_line[:-1] # remove newline character
        elif curr_line[0] == '>' and curr_seq != None:
            aligned_sequences.append(curr_seq)
            curr_seq = ''
        else:
            curr_seq = ''
        curr_line = training_data.readline()
    
    training_data.close()
    
    return aligned_sequences
    
# aligned_sequences = Align('training_data.txt') # test
aligned_sequences = [   'AG---C',
                        'A-AG-C',
                        'AG-AA-',
                        '--AAAC',
                        'AG---C'
                    ]
                    
def count_symbols(string, symbols, psuedocount=float(0)):
    symbol_counts = dict.fromkeys(symbols,psuedocount)
    for symbol in string:
        symbol_counts[symbol] += 1
    return symbol_counts
    

def build_model(aligned_sequences, psuedo_count=1, nucleotides=True):
    '''Description
    '''
    alphabet = 'ACGT' if nucleotides == True else 'GPAVLIMCFYWHKRQNEDST'
    
    states = []
    match_state_count = 3 # change to 0 after check
    for pos in range(len(aligned_sequences[0])):
        emissions = ''
        for seq in aligned_sequences:
            emissions += seq[pos]
        state = State(emissions)
        states.append(state)
        # if state.is_match_state():
            # match_state_count += 1
    # print 'STATES', ':', states, '\n'
    
    model = HMMModel(alphabet, match_state_count, psuedo_count)
    
    prev_mat_states = []
    prev_ins_states = ''
    for i in range(len(states)):
        # if states[i].is_match_state()
        if i == 0 or i == 1 or i == 5: # match state
            for emission in states[i]: # update mat emissions
                if emission != '-':
                    model['mat'][emission][len(prev_mat_states)+1] += 1
            
            for emission in prev_ins_states: # update ins emissions
                if emission != '-':
                    model['ins'][emission][len(prev_mat_states)] += 1
                    
            if len(prev_mat_states) > 0: # update state transitions
                for j in range(len(states[i])):
                    e1 = '%s' % states[prev_mat_states[-1]][j]
                    e2 = '%s' % states[i][j]
                    # print e1, '->', e2
                    
                    if e1 != '-' and e2 != '-':
                        model[('mat','mat')][prev_mat_states[-1]+1] += 1
                    elif e1 != '-' and e2 == '-':
                        model[('mat','del')][prev_mat_states[-1]] += 1
                    elif False:
                        model[('del','mat')][prev_mat_states[-1]] += 1
                    elif False:
                        model[('del','del')][prev_mat_states[-1]] += 1
                    
                # print 
                
                
                
                
            
            # cleanup
            prev_mat_states.append(i)
            prev_ins_states = ''
            
        else: # insert state
            prev_ins_states += states[i].get_emissions() # keep track until next match state
    
    return model
    
    
    
# model = build_model(aligned_sequences,0) # test
# model.print_model()





M = HMMModel('ACGT',3,1)

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

# M.print_model()


class Sequence(object):

    def __init__(self, letter):
        self.letter = letter
        self.state = None


def viterbi(model,sequence):
    M.print_model()
    M.normalize()
    M.print_model()
    # M.print_model()


    queue = deque

    for pos in range(len(sequence)):

        if pos == 0:
            lastState = "mat"

        if sequence[pos] == "-":
            currState = "del"


            # clarify: convergent or modify input sequence?
        # a -> t
        # a -> a





    pass

viterbi(M,"ACTGA")