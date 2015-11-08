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
                    print e1, '->', e2
                    
                    if e1 != '-' and e2 != '-':
                        model[('mat','mat')][prev_mat_states[-1]+1] += 1
                    elif e1 != '-' and e2 == '-':
                        model[('mat','del')][prev_mat_states[-1]] += 1
                    elif False:
                        model[('del','mat')][prev_mat_states[-1]] += 1
                    elif False:
                        model[('del','del')][prev_mat_states[-1]] += 1
                    
                print 
                
                
                
                
            
            # cleanup
            prev_mat_states.append(i)
            prev_ins_states = ''
            
        else: # insert state
            prev_ins_states += states[i].get_emissions() # keep track until next match state
    
    return model
    
    
    
model = build_model(aligned_sequences,0) # test
model.print_model()


def build_model1(aligned_sequences, psuedocount, nucleotides=True):
    '''Description
    '''
    
    symbols = 'ACGT' if nucleotides == True else 'GPAVLIMCFYWHKRQNEDST'
    states = ['match','insert','delete']
    
    total_count = len(aligned_sequences)
    length = len(aligned_sequences[0])
    
    psuedocount = float(psuedocount)
    
    emissions = dict.fromkeys(symbols, psuedocount) 
    transitions = dict.fromkeys(states, dict.fromkeys(states, psuedocount))   
     
    match_scores = [None]
    insert_scores = []
    
    transition_scores = []
    
    columns = ''
    
    prev_match_states = {}
    prev_insert_states = {}
    
    for pos in range(length):
        gap_count = 0
        match_count = dict(emissions)
        column = ''
        curr_state = ''
        for seq in aligned_sequences:
            emission = seq[pos]
            if emission == '-':
                pass
            else:
                match_count[emission] += 1
                column += seq[pos]
            curr_state += emission
        
        gap_count = curr_state.count('-')
        # match_count = total_count - gap_count
        
        if pos == 0: # initial column
            pass
        
        if match_count[max(match_count, key=match_count.get)] > total_count / 2.0: # match state
            match_scores.append(dict(emissions))
            match_scores[len(match_scores)-1] = count_symbols(column, symbols, psuedocount)
            
            insert_scores.append(dict(emissions))
            insert_scores[len(insert_scores)-1] = count_symbols(columns, symbols, psuedocount)
            
            transition_scores.append(dict(transitions))
            
            columns = ''
            prev_match_states[pos] = curr_state
            
        else: # insert state
            columns += column
            prev_insert_states[pos] = curr_state
            
        if pos == length-1: # final column 
            insert_scores.append(dict(emissions))
            insert_scores[len(match_scores)-1] = count_symbols(columns, symbols, psuedocount)
            
            transition_scores.append(dict(transitions))
            transition_scores[len(transition_scores)-1]['insert']['insert'] = 0.0
                            
    print 'match emissions', '\n'
    for pos in range(len(match_scores)):
        print pos, ':', match_scores[pos], '\n'
        
    print 'insert emissions', '\n'
    for pos in range(len(insert_scores)):
        print pos, ':', insert_scores[pos], '\n'
    
    print 'state transitions', '\n'
    for pos in range(len(transition_scores)):
        print pos, ':', transition_scores[pos], '\n'
        
    print 'match states: ', prev_match_states
    print 'insert states:', prev_insert_states