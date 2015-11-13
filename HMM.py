from HMMModel import *
from State import *
import math
import sys

NUM_STATES = 3

STATE_TABLE = {1: "mat", 2: "ins", 3: "del"}


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
        return "row: ", self.row, "col: ", self.col, "state: ", self.state, "back row: ", self.back.row, "back col: ", self.back.col, "value: ", self.value

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

# aligned_sequences = Align('training_data.txt') # test
aligned_sequences = ['AG---C',
                     'A-AG-C',
                     'AG-AA-',
                     '--AAAC',
                     'AG---C'
                     ]


def count_symbols(string, symbols, psuedocount=float(0)):
    symbol_counts = dict.fromkeys(symbols, psuedocount)
    for symbol in string:
        symbol_counts[symbol] += 1
    return symbol_counts


def build_model(aligned_sequences, psuedo_count=1, nucleotides=True):
    '''Description
    '''
    alphabet = 'ACGT' if nucleotides == True else 'GPAVLIMCFYWHKRQNEDST'

    states = []
    match_state_count = 3  # change to 0 after check
    for pos in range(len(aligned_sequences[0])):
        emissions = ''
        for seq in aligned_sequences:
            emissions += seq[pos]
        state = State(emissions)
        states.append(state)
        # if state.is_match_state():
        # match_state_count += 1
    # #print 'STATES', ':', states, '\n'

    model = HMMModel(alphabet, match_state_count, psuedo_count)

    prev_mat_states = []
    prev_ins_states = ''
    for i in range(len(states)):
        # if states[i].is_match_state()
        if i == 0 or i == 1 or i == 5:  # match state
            for emission in states[i]:  # update mat emissions
                if emission != '-':
                    model['mat'][emission][len(prev_mat_states) + 1] += 1

            for emission in prev_ins_states:  # update ins emissions
                if emission != '-':
                    model['ins'][emission][len(prev_mat_states)] += 1

            if len(prev_mat_states) > 0:  # update state transitions
                for j in range(len(states[i])):
                    e1 = '%s' % states[prev_mat_states[-1]][j]
                    e2 = '%s' % states[i][j]
                    # #print e1, '->', e2

                    if e1 != '-' and e2 != '-':
                        model[('mat', 'mat')][prev_mat_states[-1] + 1] += 1
                    elif e1 != '-' and e2 == '-':
                        model[('mat', 'del')][prev_mat_states[-1]] += 1
                    elif False:
                        model[('del', 'mat')][prev_mat_states[-1]] += 1
                    elif False:
                        model[('del', 'del')][prev_mat_states[-1]] += 1

            # cleanup
            prev_mat_states.append(i)
            prev_ins_states = ''

        else:  # insert state
            # keep track until next match state
            prev_ins_states += states[i].get_emissions()

    return model


# model = build_model(aligned_sequences,0) # test
# model.#print_model()


M = HMMModel('ACGT', 3, 1)

M['mat']['A'][1] += 4
M['mat']['C'][3] += 4
M['mat']['G'][2] += 3

M['ins']['A'][2] += 6
M['ins']['G'][2] += 1

M[('mat', 'mat')][0] += 4
M[('mat', 'mat')][1] += 3
M[('mat', 'mat')][2] += 2
M[('mat', 'mat')][3] += 4

M[('mat', 'del')][0] += 1
M[('mat', 'del')][1] += 1

M[('mat', 'ins')][2] += 1

M[('ins', 'mat')][2] += 2
M[('ins', 'del')][2] += 1
M[('ins', 'ins')][2] += 4

M[('del', 'mat')][3] += 1
M[('del', 'del')][1] += 1
M[('del', 'ins')][2] += 2


class Sequence(object):

    def __init__(self, letter):
        self.letter = letter
        self.state = None


def viterbi(model, sequence):
    model.normalize()
    matrix = [[Cell(row, col, float(0)) for col in range(len(sequence) + 1)]
              for row in range(len(model) * 3 + 2)]  # matrix

    for row in matrix:
        print row

    # print matrix

    best = []

    # for row in range(len(matrix)):
    # print matrix[row]

    curr_index = 1
    for row in range(len(matrix)):
        for col in range(len(matrix[row])):
            if row > 0 and col > 0 and row < len(model) * 3 + 1:

                # print row,col, '\n', '---'
                curr = matrix[row][col]

                curr.state = STATE_TABLE[(row - 1) % NUM_STATES + 1]

                if row < 4 and col < 1:  # from begin state
                    pass

                elif curr.state == "mat":  # match state
                    # need to fill the matrix now
                    emission_prob = math.log(
                        model['mat'][sequence[col - 1]][curr_index])
                    from_mat = matrix[row - 3][col - 1].value + \
                        math.log(model[('mat', 'mat')][curr_index - 1])
                    from_ins = matrix[row - 2][col - 1].value + \
                        math.log(model[('ins', 'mat')][curr_index - 1])
                    # print model[('del','mat')][curr_index]

                    if curr_index == 1:
                        from_del = 0
                    else:
                        from_del = matrix[row - 1][col - 1].value + \
                            math.log(model[('del', 'mat')][curr_index - 1])


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


                elif curr.state == "ins":  # insert state
                    # print model['ins'][sequence[col-1]]
                    # print row

                    emission_prob = math.log(
                        model['ins'][sequence[col - 1]][curr_index])
                    from_mat = matrix[row - 1][col - 1].value + \
                        math.log(model[('mat', 'ins')][curr_index - 1])
                    from_ins = matrix[row][
                        col - 1].value + math.log(model[('ins', 'ins')][curr_index - 1])

                    if curr_index == 1:
                        from_del = 0
                    else:
                        from_del = matrix[
                            row + 1][col - 1].value + math.log(model[('del', 'ins')][curr_index - 1])

                    # print 'hererrerererererer',
                    # model[('del','ins')][curr_index-1]


                    # match was the best
                    if from_mat >= from_ins and from_mat >= from_del:
                        curr.back = matrix[row - 1][col - 1]
                        curr.value = emission_prob + from_mat

                    # insertion was the best
                    elif from_ins > from_mat and from_ins > from_del:
                        curr.back = matrix[row][col - 1]
                        curr.value = emission_prob + from_ins

                    # delete was the best
                    elif from_del > from_mat and from_del >= from_ins:
                        curr.back = matrix[row + 1][col - 1]
                        curr.value = emission_prob + from_del

                    else:
                        print "problem comparing!"
                        sys.exit(1)


                elif curr.state == "del":  # delete state
                    if curr_index > 1:
                        from_mat = matrix[
                            row - 5][col].value + math.log(model[('mat', 'del')][curr_index - 1])
                        from_ins = matrix[
                            row - 4][col].value + math.log(model[('ins', 'del')][curr_index - 1])
                        from_del = matrix[
                            row - 3][col].value + math.log(model[('del', 'del')][curr_index - 1])



                    # match was the best
                    if from_mat >= from_ins and from_mat >= from_del:
                        curr.back = matrix[row - 5][col - 1]
                        curr.value = emission_prob + from_mat

                    # insertion was the best
                    elif from_ins > from_mat and from_ins > from_del:
                        curr.back = matrix[row - 4][col - 1]
                        curr.value = emission_prob + from_ins

                    # delete was the best
                    elif from_del > from_mat and from_del >= from_ins:
                        curr.back = matrix[row - 3][col - 1]
                        curr.value = emission_prob + from_del

                    else:
                        print "problem comparing!"
                        sys.exit(1)

                    if col == len(sequence) - 1:
                        curr_index += 1


                        # print "incrementation", curr_index
            else:  # to end state
                pass

            # print
    for row in matrix:
        print row


viterbi(M, "ACTGAT")
# viterbi(M,"TAGATTG")
