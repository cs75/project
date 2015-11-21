from HMMModel import *
from State import *
import math
import sys
import numpy as np  

def Align(data):
    '''Description
    '''
    training_data = open(data, 'r')
    aligned_sequences = []

    curr_seq = None
    curr_line = training_data.readline()

    while curr_line != '':  # populate aligned_sequences from training_data
        if curr_line[0] != '>':
            curr_seq += curr_line[:-1]  # remove newline character
        elif curr_line[0] == '>' and curr_seq != None:
            aligned_sequences.append(curr_seq)
            curr_seq = ''
        else:
            curr_seq = ''
        curr_line = training_data.readline()

    training_data.close()

    return aligned_sequences


def build_model(aligned_sequences, alphabet, psuedo_count):

    states = [State(0, 'beg', None, None, None, None)]
    # intialize the states
    prev_state = states[-1]
    prev_mat_state = None
    prev_ins_state = None

    curr_model_position = 0

    prev_ins_emissions = ''

    # WIDTH = len(sorted(aligned_sequences, key=len)[-1])


    count = 0
    WIDTH = len(aligned_sequences[0])

    for col in range(WIDTH + 1):  # +1 accounts for end state
        if col < len(aligned_sequences[0]):
            emissions = ''
            for sequence in aligned_sequences:
                # for each sequence, input the emissions
                emissions += sequence[col]

            if emissions.count('-') / float(len(aligned_sequences)) < 0.50:
                # is there an insert state?
                if len(prev_ins_emissions) > 0:
                    # instantiate
                    state = State(curr_model_position, 'ins', str(
                        prev_ins_emissions), prev_state, prev_mat_state, prev_ins_state)
                    prev_ins_emissions = ''
                    states.append(state)
                    # update the general and specific "previous" states
                    prev_state = state
                    prev_ins_state = state
                # increment the position, since the position corresponds to a
                # match state... that is you can only have a position in the
                # model if it is a match state
                curr_model_position += 1
                # if there was not an insert state before it, instantiate
                state = State(curr_model_position, 'mat', str(
                    emissions), prev_state, prev_mat_state, prev_ins_state)
                states.append(state)

                prev_state = state
                prev_mat_state = state

            else:                                      # insert state; wait until next mat state unless last col
                prev_ins_emissions += emissions
                if col == WIDTH - 1:  # last col
                    state = State(curr_model_position, 'ins', prev_ins_emissions,
                                  prev_state, prev_mat_state, prev_ins_state)
                    prev_ins_emissions = ''
                    states.append(state)

                    prev_state = state
                    prev_ins_state = state

        else:  # end state
            states.append(State(curr_model_position + 1, 'end',
                                None, prev_state, prev_mat_state, prev_ins_state))

    model = HMMModel(alphabet, curr_model_position, psuedo_count)

    # for state in states:
    #     print state.label,
    # print

    for curr_state in states:
        prev_state = curr_state.prev_state  # mat or ins
        prev_mat_state = curr_state.prev_mat_state  # prior to or equal to prev_state
        prev_ins_state = curr_state.prev_ins_state  # prior to or equal to prev_state

        if curr_state.label == 'mat':
            for emission in curr_state.emissions:           # mat emission counts
                if emission != '-':
                    model['mat'][emission][curr_state.position] += 1

            if prev_mat_state != None:                      # transition counts only after first mat state
                # curr_state will be a match state due to prior check
                for i in range(len(curr_state.emissions)):
                    prev_emission = prev_mat_state.emissions[i]
                    curr_emission = curr_state.emissions[i]

                    if prev_emission != '-' and curr_emission != '-':  # possibly gapped mat to mat for mat counts
                        model[('mat', 'mat')][prev_mat_state.position] += 1

                    if prev_state == prev_mat_state:  # immediate mat to mat for del counts
                        if prev_emission != '-' and curr_emission == '-':
                            model[('mat', 'del')][prev_mat_state.position] += 1
                        elif prev_emission == '-' and curr_emission == '-':
                            model[('del', 'del')][prev_mat_state.position] += 1
                        elif prev_emission == '-' and curr_emission != '-':
                            model[('del', 'mat')][prev_mat_state.position] += 1

                if prev_state == prev_ins_state:  # ins -> mat

                    split_emissions = [prev_state.emissions[i:i + len(curr_state.emissions)] for i in range(
                        0, len(prev_state.emissions), len(curr_state.emissions))]

                    for i in range(len(curr_state.emissions)):
                        curr_emission = curr_state.emissions[i]
                        ins_to_mat = False
                        ins_to_del = False

                        for prev_emissions in split_emissions:
                            prev_emission = prev_emissions[i]
                            if prev_emission != '-' and curr_emission != '-':
                                ins_to_mat = True
                            elif prev_emission != '-' and curr_emission == '-':
                                ins_to_del = True

                        if ins_to_mat:
                            model[('ins', 'mat')][prev_ins_state.position] += 1
                        elif ins_to_del:
                            model[('ins', 'del')][prev_ins_state.position] += 1

            else:                                   # beg state since there was no previous match state
                if prev_ins_state == None:             # -> mat state

                    for emission in curr_state.emissions:
                        if emission != '-':
                            model[('mat', 'mat')][0] += 1
                        else:
                            model[('mat', 'del')][0] += 1
                else:                                  # -> ins state
                    for emission in prev_ins_state.emissions:
                        if emission != '-':
                            model[('mat', 'ins')][0] += 1
                        else:
                            model[('mat', 'del')][0] += 1

        # ins to ins
        elif curr_state.label == 'ins':
            split_emissions = [curr_state.emissions[i:i + WIDTH]
                               for i in range(0, len(curr_state.emissions), WIDTH)]
            for i in range(len(split_emissions) - 2):
                for j in range(len(split_emissions[i])):
                    if split_emissions[i][j] != '-' and split_emissions[i + 1][j] != '-':
                        model[('ins', 'ins')][curr_state.position] += 1
            # print split_emissions

            # print split_emissions

            # match to insert
            if prev_state.label == "mat":
                model[('mat', 'ins')][prev_mat_state.position] += 1

            for emission in curr_state.emissions:           # ins emission counts
                if emission != '-':
                    model['ins'][emission][curr_state.position] += 1

            for i in range(len(prev_mat_state.emissions)):
                prev_emission = prev_mat_state.emissions[i]
                curr_emission = curr_state.emissions[i]

                if prev_emission == '-' and curr_emission != '-':
                    model[('del', 'ins')][prev_mat_state.position] += 1

        elif curr_state.label == 'end':  # -> end state
            if prev_state == prev_mat_state:        # mat state ->
                for emission in prev_mat_state.emissions:
                    if emission != '-':
                        model[('mat', 'mat')][prev_mat_state.position] += 1
                    else:
                        model[('del', 'mat')][prev_mat_state.position] += 1
            else:                                   # ins state ->
                for emission in prev_ins_state.emissions:
                    if emission != '-':
                        model[('ins', 'mat')][prev_ins_state.position] += 1

    return model

class WavefrontCell(object):
    """
    Represents a cell in the wavefront matrix

    row
    col
    state = {"mat", "del", "ins"} -> val
    back = cell we came from (best cell)
    """

    def __init__(self, row, col):
        self.row = row
        self.col = col
        self.states = None
        self.backs = None

    def __str__(self):
        return "row: " + str(self.row) + " col: " + str(self.col) + " state: "

        #, "back row: ", self.back.row, "back col: ", self.back.col, "value: ", self.value

    def __repr__(self):
        return str('self.states')
        

def standard_viterbi(model,sequence,alphabet):
    #populating the matrix here 
    matrix = [[WavefrontCell(row,col) for col in range(len(sequence)+2)] for row in range(len(model)+2)]
    bg_distribution = 1.0 / (len(alphabet) + 1)
    model.normalize()

    curr = None # eventually set to final cell for backtracking

    for row in range(len(matrix)):
        for col in range(len(matrix[0])):
            if (row,col) == (0,0):                                      # start case
                matrix[row][col].states = {'mat': float(0)}
            elif row == 0 and col < len(matrix[0])-1:                   # border cases: horizontal
                matrix[row][col].states = {}
                for prev_kind in matrix[row][col-1].states:
                    matrix[row][col].states['ins'] = matrix[row][col-1].states[prev_kind] + math.log(model[(prev_kind,'ins')][row])
                    matrix[row][col].backs = {'ins': prev_kind}
            elif col == 0 and row < len(matrix)-1:                      # border cases: vertical
                matrix[row][col].states = {}
                for prev_kind in matrix[row-1][col].states:
                    matrix[row][col].states['del'] = matrix[row-1][col].states[prev_kind] + math.log(model[(prev_kind,'del')][row-1])
                    matrix[row][col].backs = {'del': prev_kind}
            elif row < len(matrix)-1 and col < len(matrix[0])-1:        # normal cases: everything but last
                curr_cell = matrix[row][col]
                curr_cell.states = {}
                curr_cell.backs = {}
                
                from_mat_cell = matrix[row-1][col-1]
                from_ins_cell = matrix[row][col-1]
                from_del_cell = matrix[row-1][col]
        
                from_cells = {'mat': from_mat_cell, 'ins': from_ins_cell, 'del': from_del_cell}

                for curr_kind in ['mat','ins','del']:
                    from_cell = from_cells[curr_kind]
                    best_score = None
                    best_back = None
                    for prev_kind in from_cell.states:
                        score = from_cell.states[prev_kind] + math.log(model[(prev_kind,curr_kind)][from_cell.row])
                        if curr_kind == 'mat' and row < len(matrix)-1:
                            score += math.log(model[curr_kind][sequence[col-1]][row] / bg_distribution)
                        
                        if best_score == None or score > best_score:
                            best_score = score
                            best_back = prev_kind
                    
                    curr_cell.states[curr_kind] = best_score
                    curr_cell.backs[curr_kind] = best_back
            elif (row,col) == (len(matrix)-1,len(matrix[0])-1): # final cell
                final_cell = matrix[row][col]
                best_score = None
                best_back = None

                for prev_kind in matrix[row-1][col-1].states:
                    score = matrix[row-1][col-1].states[prev_kind] + math.log(model[(prev_kind,'mat')][row-1])
                    if best_score == None or score > best_score:
                        best_score = score
                        best_back = prev_kind

                final_cell.states = {'mat': best_score}
                final_cell.backs = {'mat': best_back}
                curr = final_cell
                
    best_path = []
    print 'FINAL SCORE =', curr.states['mat']
    while curr.backs != None: # back tracking
        best_kind = max(curr.states, key=curr.states.get)
        best_back = curr.backs[best_kind]

        if curr.col != len(matrix[0])-1 and (curr.row,curr.col) != (0,0):
            best_path.insert(0,best_kind)

        back_mat_cell = matrix[max(0,curr.row-1)][max(0,curr.col-1)]
        back_ins_cell = matrix[max(0,curr.row)][max(0,curr.col-1)]
        back_del_cell = matrix[max(0,curr.row-1)][max(0,curr.col)]
        back_cells = {'mat': back_mat_cell, 'ins': back_ins_cell, 'del': back_del_cell}
        
        curr = back_cells[best_kind]

    print best_path

def diagonal_viterbi(model,sequence):
    matrix = np.array([[WavefrontCell(row,col) for col in range(len(sequence)+2)] for row in range(len(model)+2)])
    bg_distribution = 1.0 / (len(alphabet) + 1)
    model.normalize()

    wavefront = [matrix[::-1,:].diagonal(i) for i in range(-matrix.shape[0]+1, matrix.shape[1])]
    
    curr = None  # eventually set to final cell for backtracking
    
    # print len(matrix), 'by', len(matrix[0])
    # for i in range(len(wavefront)):
    #     for j in range(len(wavefront[i])):
    #         print (wavefront[i][j].row, wavefront[i][j].col),
    #     print

    for i in range(len(wavefront)):
        for j in range(len(wavefront[i])):
            if i == 0:                                          # beg state
                wavefront[i][j].states = {'mat':float(0)}
        
            elif i == len(wavefront)-1:                            # end state
                final_cell = matrix[row][col]
                best_score = None
                best_back = None

                for prev_kind in matrix[row-1][col-1].states:
                    score = matrix[row-1][col-1].states[prev_kind] + math.log(model[(prev_kind,'mat')][row-1])
                    if best_score == None or score > best_score:
                        best_score = score
                        best_back = prev_kind

                final_cell.states = {'mat': best_score}
                final_cell.backs = {'mat': best_back}

                curr = final_cell

            elif wavefront[i][j].col == 0:
                if wavefront[i][j].row < len(matrix)-1:         # del state
                    wavefront[i][j].states = {}
                    for prev_kind in wavefront[i-1][0].states:
                        wavefront[i][j].states['del'] = wavefront[i-1][0].states[prev_kind] + math.log(model[(prev_kind,'del')][wavefront[i][j].row])
                        wavefront[i][j].backs = {'del': prev_kind}
            
            elif wavefront[i][j].row == 0:
                if wavefront[i][j].col < len(matrix[0])-1:      # ins state
                    wavefront[i][j].states = {}
                    for prev_kind in wavefront[i-1][-1].states:
                        wavefront[i][j].states['ins'] = wavefront[i-1][-1].states[prev_kind] + math.log(model[(prev_kind,'ins')][wavefront[i][j].row])
                        wavefront[i][j].backs = {'ins': prev_kind}

            elif wavefront[i][j].row < len(matrix)-1 and wavefront[i][j].col < len(matrix[0])-1:
                microwave = wavefront[i][j]
                microwave.states = {}
                microwave.backs = {}
                
                if len(wavefront[i-1]) < len(wavefront[i]): 
                    from_mat = wavefront[i-2][j-1]
                    from_ins = wavefront[i-1][j-1]
                    from_del = wavefront[i-1][j]
                elif len(wavefront[i-1]) > len(wavefront[i]): 
                    from_mat = wavefront[i-2][j]
                    from_ins = wavefront[i-1][j]
                    from_del = wavefront[i-1][j+1]
                else: 
                    from_mat = wavefront[i-2][j]
                    from_ins = wavefront[i-1][j]
                    from_del = wavefront[i-1][j+1]
                
                from_microwaves = {'mat': from_mat, 'ins': from_ins, 'del': from_del}

                for kind in ['mat','ins','del']:
                    from_microwave = from_microwaves[kind]
                    best_score = None
                    best_back = None
                    for prev_kind in from_microwave.states:
                        score = from_microwave.states[prev_kind] + math.log(model[(prev_kind,kind)][from_microwave.row])
                        if kind == 'mat':
                            score += math.log(model[kind][sequence[microwave.col-1]][microwave.row] / bg_distribution)
                        if best_score == None or score > best_score:
                            best_score = score
                            best_back = prev_kind

                    microwave.states[kind] = best_score
                    microwave.backs[kind] = best_back

    for i in range(len(wavefront)):
        for j in range(len(wavefront[i])):
            print (wavefront[i][j].row, wavefront[i][j].col, wavefront[i][j].states),
        print

    best_path = []
    print 'FINAL SCORE =', curr.states['mat']
    while curr.backs != None: # back tracking
        best_kind = max(curr.states, key=curr.states.get)
        best_back = curr.backs[best_kind]

        if curr.col != len(matrix[0])-1 and (curr.row,curr.col) != (0,0):
            best_path.insert(0,best_kind)

        back_mat_cell = matrix[max(0,curr.row-1)][max(0,curr.col-1)]
        back_ins_cell = matrix[max(0,curr.row)][max(0,curr.col-1)]
        back_del_cell = matrix[max(0,curr.row-1)][max(0,curr.col)]
        back_cells = {'mat': back_mat_cell, 'ins': back_ins_cell, 'del': back_del_cell}
        
        curr = back_cells[best_kind]

    print best_path
    


if __name__ == "__main__":

    # AAs
    # alphabet = 'GPAVLIMCFYWHKRQNEDST'
    # NTs
    alphabet = 'AGCT'
    # aligned_sequences = Align("./tests/PF00005_seed.txt")
    aligned_sequences = [   'AG---C',
                            'A-AG-C',
                            'AG-AA-',
                            '--AAAC',
                            'AG---C'
                        ]

    model = build_model(aligned_sequences, alphabet, 1)

    sequence = "AGAAACCCCCC"
    #sequence = "AGAAAAAAAAAAAAAAAAAAAAAAAAC"

    diagonal_viterbi(model, sequence)
    standard_viterbi(model, alphabet, sequence)