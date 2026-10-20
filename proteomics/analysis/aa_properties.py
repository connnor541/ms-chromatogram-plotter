from Bio.SeqUtils.ProtParam import ProteinAnalysis

def compute_peptide_properties(sequence: str):

    clean_seq = sequence.strip().upper()
    analysed_seq = ProteinAnalysis(clean_seq)

    mw_value = analysed_seq.molecular_weight()
    pi_value = analysed_seq.isoelectric_point()
    gravy_value = analysed_seq.gravy()

    return mw_value, pi_value, gravy_value


