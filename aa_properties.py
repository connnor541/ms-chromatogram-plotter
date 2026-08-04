from Bio.SeqUtils.ProtParam import ProteinAnalysis

def compute_peptide_properties(sequence: str):

    clean_seq = sequence.strip().upper()
    analysed_seq = ProteinAnalysis(clean_seq)

    mw_value = analysed_seq.molecular_weight()
    pi_value = analysed_seq.isoelectric_point()
    gravy_value = analysed_seq.gravy()

    return mw_value, pi_value, gravy_value


def analyze_peptide(sequence: str):
    # Clean up the sequence (remove whitespaces and convert to uppercase)
    clean_seq = sequence.strip().upper()
    
    # Initialize the ProteinAnalysis class
    analysed_seq = ProteinAnalysis(clean_seq)
    
    # 1. Calculate Theoretical pI
    pi_value = analysed_seq.isoelectric_point()
    
    # 2. Calculate GRAVY (Grand Average of Hydropathy)
    gravy_value = analysed_seq.gravy()
    
    # 3. Calculate Molecular Weight
    mw_value = analysed_seq.molecular_weight()
    
    # Print results
    print(f"Sequence: {clean_seq}")
    print(f"Length:   {len(clean_seq)}")
    print(f"MW:       {mw_value:.2f} Da")
    print(f"pI:       {pi_value:.2f}")
    print(f"GRAVY:    {gravy_value:.2f}")
    print("-" * 30)

# Example usage with a sample peptide sequence (e.g., a custom peptide)
if __name__ == "__main__":
    # Example with multiple sequences
    sequences = ["ACDEFGHIKLMNPQRSTVWY", "GLVTR"]
    for seq in sequences:
        analyze_peptide(seq)