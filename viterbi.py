from HMMModel import *
from State import *
import math
import sys

NUM_STATES = 3

STATE_TABLE = {1: "mat", 2: "ins", 0: "del"}


class Cell(object):
    """
    Represents a cell in the matrix

    row
    col
    state = "mat" "del" "ins"
    back = cell we came from (best cell)
    """

    def __init__(self, row, col, value):
        self.row = row
        self.col = col
        self.state = None
        self.back = None
        self.value = value

    def __str__(self):
        return "row: " + str(self.row) + " col: " + str(self.col) + " state: " + str(self.state)

        #, "back row: ", self.back.row, "back col: ", self.back.col, "value: ", self.value

    def __repr__(self):
        return str(self.value)


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

    WIDTH = len(aligned_sequences)
    count = 0

    for col in range(WIDTH + 1):  # +1 accounts for end state
        if col < len(aligned_sequences[0]):
            emissions = ''
            for sequence in aligned_sequences:
                # for each sequence, input the emissions
                emissions += sequence[col]

            if emissions.count('-') / float(len(aligned_sequences)) < 0.50:

                # if col == 0 or col == 1 or col == 5:        # match state
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

                if col == len(aligned_sequences[0]) - 1:  # last col
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
            for i in range(len(split_emissions) - 1):
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


def count_symbols(string, symbols, psuedocount=float(0)):
    symbol_counts = dict.fromkeys(symbols, psuedocount)
    for symbol in string:
        symbol_counts[symbol] += 1
    return symbol_counts


# model = build_model(aligned_sequences,0) # test
# model.#print_model()


# M = HMMModel('ACGT', 3, 1)

# M['mat']['A'][1] += 4
# M['mat']['C'][3] += 4
# M['mat']['G'][2] += 3

# M['ins']['A'][2] += 6
# M['ins']['G'][2] += 1

# M[('mat', 'mat')][0] += 4
# M[('mat', 'mat')][1] += 3
# M[('mat', 'mat')][2] += 2
# M[('mat', 'mat')][3] += 4

# M[('mat', 'del')][0] += 1
# M[('mat', 'del')][1] += 1

# M[('mat', 'ins')][2] += 1

# M[('ins', 'mat')][2] += 2
# M[('ins', 'del')][2] += 1
# M[('ins', 'ins')][2] += 4

# M[('del', 'mat')][3] += 1
# M[('del', 'del')][1] += 1
# M[('del', 'ins')][2] += 2


# class Sequence(object):

#     def __init__(self, letter):
#         self.letter = letter
#         self.state = None

# def find_best()

def inside_matrix(row, num_rows):
    return row < num_rows and row >= 0


def viterbi(model, sequence, alphabet):

    bg_distribution = (1.0 / (len(alphabet) + 1)) # e.g. 1 in (num letters + one gap) random chance

    model.normalize()
    model.print_model()
    matrix = [[Cell(row, col, float(0)) for col in range(len(sequence) + 1)]
              for row in range(len(model) * 3 + 2)]  # matrix



    best = []

    curr_index = 1

    already_incremented = False


    for row in range(len(matrix)):
        for col in range(len(matrix[0])):

            # skip the begin state (which is row 0)
            if row == 0:
                matrix[row][col].state = "mat"
                continue

            # skip the end state (which is last row)
            elif row == len(matrix) - 1 and col < len(sequence):
                continue

            curr = matrix[row][col]
            curr.state = STATE_TABLE[row % NUM_STATES]

            
            # need col to be incremented by 1
            if col == 0:
                matrix[row][col].state = curr.state 

            #################################
            ## BEGIN STATE handler block
            #################################

            elif curr.state == "mat" and row == 1:
                print sequence[col - 1]
                print model['mat'][sequence[col - 1]][curr_index]
                emission_prob = math.log(model['mat'][sequence[col - 1]][curr_index] / (bg_distribution))
                print emission_prob

                curr.back = matrix[row - 1][col - 1]
                curr.value = emission_prob + math.log(model[('mat', 'mat')][curr_index - 1])

            elif curr.state == "ins" and row == 2:
                curr.back = matrix[row - 2][col - 1]
                curr.value = math.log(model[('mat', 'ins')][curr_index - 1])

            elif curr.state == "del" and row == 3:
                curr.back = matrix[row - 3][col - 1]
                curr.value = math.log(model[('mat', 'del')][curr_index - 1])

                # Hack to get around begin state edge case
                if not already_incremented:
                    curr_index += 1
                    already_incremented = True



            #################################
            ## END STATE handler block
            #################################


            elif row == len(matrix) - 1 and col == len(sequence):

                from_mat = matrix[row - 3][col].value
                from_ins = matrix[row - 2][col].value
                from_del = matrix[row - 1][col].value

                # match was the best
                if from_mat >= from_ins and from_mat >= from_del:
                    curr.back = matrix[row - 3][col]
                    curr.value = from_mat

                # insertion was the best
                elif from_ins > from_mat and from_ins > from_del:
                    curr.back = matrix[row - 2][col]
                    curr.value = from_ins

                # delete was the best
                elif from_del > from_mat and from_del >= from_ins:
                    curr.back = matrix[row - 1][col]
                    curr.value = from_del

                else:
                    print "problem comparing!"
                    sys.exit(1)



            #################################
            ## REGULAR STATE handler block
            ## these are the regular cases in between begin and end states! (we're at >= row 4)
            #################################

            # match state
            elif curr.state == "mat":
                emission_prob = math.log(model['mat'][sequence[col - 1]][curr_index] / (bg_distribution))

                from_mat = matrix[row - 3][col - 1].value + math.log(model[('mat', 'mat')][curr_index - 1]) 
                from_ins = matrix[row - 2][col - 1].value + math.log(model[('ins', 'mat')][curr_index - 1])
                from_del = matrix[row - 1][col - 1].value + math.log(model[('del', 'mat')][curr_index - 1])

                # match was the best
                if from_mat >= from_ins and from_mat >= from_del:
                    curr.back = matrix[row - 3][col - 1]
                    curr.value = emission_prob + from_mat

                # insertion was the best
                elif from_ins > from_mat and from_ins > from_del:
                    curr.back = matrix[row - 2][col - 1]
                    curr.value = emission_prob + from_ins

                # delete was the best
                elif from_del > from_mat and from_del >= from_ins:
                    curr.back = matrix[row - 1][col - 1]
                    curr.value = emission_prob + from_del

                else:
                    print "problem comparing!"
                    sys.exit(1)


            # insert state
            elif curr.state == "ins":  
                from_mat = matrix[row - 1][col - 1].value + math.log(model[('mat', 'ins')][curr_index - 1])
                from_ins = matrix[row][col - 1].value + math.log(model[('ins', 'ins')][curr_index - 1])
                from_del = matrix[row + 1][col - 1].value + math.log(model[('del', 'ins')][curr_index - 1])

                # match was the best
                if from_mat >= from_ins:
                    curr.back = matrix[row - 1][col - 1]
                    curr.value = from_mat

                # insertion was the best
                else:
                    curr.back = matrix[row][col - 1]
                    curr.value = from_ins


            elif curr.state == "del":  # delete state
                from_mat = matrix[row - 5][col].value + math.log(model[('mat', 'del')][curr_index - 1])
                from_ins = matrix[row - 4][col].value + math.log(model[('ins', 'del')][curr_index - 1])
                from_del = matrix[row - 3][col].value + math.log(model[('del', 'del')][curr_index - 1])

                # match was the best
                if from_mat >= from_del:
                    curr.back = matrix[row - 5][col]
                    curr.value = from_mat

                # delete was the best
                else:
                    curr.back = matrix[row - 3][col]
                    curr.value = from_del

                if col == len(sequence) - 1:
                    curr_index += 1

    for row in matrix:
        print row

    curr = matrix[len(matrix) - 1][len(sequence)]

    while curr:
        print "backtracking to row: ", curr.row, "col: ", curr.col, "with value: ", curr.value, "and state is: ", curr.state
        curr = curr.back


if __name__ == "__main__":

    # alphabet = 'AGCT'
    # aligned_sequences = ['AG---C',
    #                      'A-AG-C',
    #                      'AG-AA-',
    #                      '--AAAC',
    #                      'AG---C'
    #                      ]

    alphabet = 'GPAVLIMCFYWHKRQNEDST'

    # aligned_sequences = Align("./training_data.txt")

    aligned_sequences = ['SPADKTNVKAAWGKVGA--HAGEYGAEALERMFLS',
                     'TPEEKSAVTALWGKV----NVDEVGGEALGRLLVV',
                     'SEGEWQLVLHVWAKVEA--DVAGHGQDILIRLFKS',
                     'SADQISTVQASFDKVKG------DPVGILYAVFKA',
                     'SAAEKTKIRSAWAPVYS--TYETSGVDILVKFFTS'
                     ]

    model = build_model(aligned_sequences, alphabet, 1)
    # model.normalize()
    # model.print_model()

    # viterbi(model, "GAT", alphabet)
    viterbi(model, "CCCCCC", alphabet)
# viterbi(M,"TAGATTG")
