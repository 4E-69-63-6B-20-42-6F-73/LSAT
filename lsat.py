# from io import StringIO
from io import StringIO

import pynecone as pc

from Bio import SeqIO
from Bio.Align import substitution_matrices
from Bio import Align

# Biopython doesn't have a blosum 65 matrix.
# Created manually from: https://github.com/biogo/biogo/blob/master/align/matrix/matrices/BLOSUM65
# NOTE: permissions must be set correctly
matrix_options = ["BLOSUM65"]

class State(pc.State):
    """The app state."""
    fasta = ""
    processing = False
    empty = True
    result = []
    matrix: str = matrix_options[0]
    ids = []
    gap_penalty : int = 0
    gap_extension_penalty : int = 0

    def process_sequence(self):
        self.processing = True
        self.empty = False

    def get_result(self):
        # Check if the input is in fasta
        if(">" not in self.fasta):
            #It is not fasta, so assume every line contains a sequence.
            real_fasta = ""
            sequences = self.fasta.split()
            for i in range(len(sequences)):
                real_fasta += ">seq" + str(i) + "\n"
                real_fasta += sequences[i].upper() + "\n"
            self.fasta = real_fasta

        # # workaround since Biopython can't create seqences from a string
        sequences_handle = StringIO(self.fasta)
        records = list(SeqIO.parse(sequences_handle, "fasta"))
        self.ids = [str(x.id) for x in records]

        aligments = do_aligment(records, self)
        
        #Create the values so that it can be displayed in a datatable
        #Add the first row
        self.result = [[""]]
        self.result[0].extend(self.ids)

        for i in range(len(records)):
            temp = [self.ids[i]]
            temp.extend(aligments[i])
            self.result.append(temp)
        
        self.processing = False

def do_aligment(fasta: list, state: State):
    result = [ [0]*len(fasta) for i in range(len(fasta))]

    open_gap_penalty = int(state.gap_penalty)
    extend_gap_penalty = int(state.gap_extension_penalty)
    matrix = substitution_matrices.load(state.matrix)

    for i in range(len(fasta)):
        for j in range(len(fasta)):
            sequence1 = fasta[i].seq
            sequence2 = fasta[j].seq

            aligner = Align.PairwiseAligner()
            aligner.mode = "local"
            aligner.substitution_matrix = matrix
            aligner.open_gap_score = open_gap_penalty
            aligner.extend_gap_score = extend_gap_penalty

            result[i][j] = aligner.score(sequence1, sequence2)

            # alignments = pairwise2.align.globalds(sequence1.seq, sequence2.seq ,matrix , int(open_gap_penalty), int(extend_gap_penalty))
            # result[i][j] = round(alignments[0].score, 2)
        
    return result

def index():
    return pc.center(
        pc.vstack(
            pc.vstack(
                pc.heading("L   ocal Sequence Aligment Table", font_size="1.5em"),
                pc.text_area(
                    placeholder="Enter one or more sequences line by line or in fasta format... ", 
                    on_blur=State.set_fasta,
                    height="25em" ),
                pc.form_control(
                    pc.heading("Options", font_size="1.5em"),
                    pc.form_label("Matrix"),
                    pc.select(
                        matrix_options,
                        on_change=State.set_matrix,
                        margin_bottom="1.25em"
                        ),
                    pc.form_label("Gap penalty"),
                    pc.number_input(
                        default_value = -10,
                        on_change=State.set_gap_penalty,
                        margin_bottom="1.25em"
                    ),
                    pc.form_label("Gapextension penalty"),
                    pc.number_input(
                        default_value = -1,
                        on_change=State.set_gap_penalty,
                        margin_bottom="1.25em"
                    ),
                ),
                pc.divider(),
                pc.button(
                    "Align sequences",
                    on_click=[State.process_sequence, State.get_result],
                    width="100%",
                ),
                bg="white",
                padding="2em",
                shadow="lg",
                border_radius="lg",
                width="50em"
            ),

            pc.vstack(
                pc.vstack(
                    pc.heading("Result", font_size="1.5em"),
                    pc.cond(
                        State.empty,
                        pc.text("No aligment yet"),
                        pc.cond(
                            State.processing,
                            pc.circular_progress(is_indeterminate=True),
                            pc.data_table(
                                data=State.result,
                                pagination=False,
                                search=False,
                                sort=False,
                            )
                        ),      
                    )  
                ),
                bg="white",
                padding="2em",
                shadow="lg",
                border_radius="lg",
                min_width="50em"
            ),
            padding="5em",
            width="100%",
            bg="linear-gradient(90deg, rgba(9,92,121,1) 0%, rgba(16,0,255,0.40379901960784315) 100%);",
        )
    )

# Add state and page to the app.
app = pc.App(state=State)
app.add_page(index, title="LSAT")
app.compile()