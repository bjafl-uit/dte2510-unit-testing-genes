"""O3.2: Unit testing of find genes.

+----------------------------+
| DTE-2510 FA24 (UiT)        |
|                            |
| Submission:   O3           |
| Title:        Find Genes   |
| Student:      Bjarte Lode  |
|                            |
| Revised version of code    |
| Originally submitted in O2 |
+----------------------------+


Description:
------------
Simple script to find genes in a DNA sequence. The script will read a DNA
sequence from user input, and search for genes in the sequence. The genes are
defined as a sequence of triplets, starting with ATG and ending with one of the
triplets TAA, TAG, or TGA. Valid sequences may not contain ATG. The script will
print the genes found in the DNA sequence, separated by commas.
"""


def find_genes(dna: str, strict_align: bool = False) -> list[str]:
    """Find genes in a DNA sequence.

    The function will search for valid genes in the DNA sequence, defined as a
    sequence that start with ATG and end with one of the sequences TAA, TAG, or
    TGA. The genes returned will not include the start and stop sequences.

    If strict_align parameter is enabled, the function will only
    interpret codons aligned in multiples of three, once the first start codon
    has been found. If not enabled, start codons that are not aligned in a
    multiple of three with the first start codon found will be recognized as a
    valid start for a new gene sequence.

    Args:
        dna (str): DNA sequence. May only contain the letters A, C, T and G.
        strict_align (bool): Enable strict triplet alignment.

    Returns:
        list[str]: List of genes found in the sequence. If no genes are found,
            the list will contain the message "No genes found".

    Raises:
        TypeError: If input is not a string
        ValueError: If input contains invalid characters
    """
    VALID_CODES = "ATCG"
    START_SEQ = "ATG"
    STOP_SEQS = ["TAA", "TAG", "TGA"]
    NO_GENES_MSG = "No genes found"

    if isinstance(dna, str) is False:
        raise TypeError("Input must be a string")
    dna = dna.strip().upper()
    if len(dna.strip(VALID_CODES)) > 0:
        raise ValueError("Invalid DNA sequence. The string may only contain " +
                         "the letters A, C, T and G")

    valid_genes = []
    first_seq_start_pos = seq_start_pos = dna.find(START_SEQ)
    while seq_start_pos != -1 and seq_start_pos < len(dna) - 3:
        cur_seq = []
        # Start from the next triplet (+3), excluding start codon from output.
        # Iterate over the next triplets (i+=3 each iteration)
        for cur_pos in range(seq_start_pos + 3, len(dna), 3):
            triplet = dna[cur_pos:cur_pos+3]
            if triplet in STOP_SEQS:
                if cur_seq:
                    # Store valid sequence as string with no delimiter
                    # Exclude stop codon from output
                    valid_genes.append("".join(cur_seq))
                break
            if triplet == START_SEQ:
                # Start over if another start codon is found
                # Current sequence is invalid
                break
            cur_seq.append(triplet)
        # Find next start codon from (and including) current position
        seq_start_pos = dna.find(START_SEQ, cur_pos)
        # If strict alignment of triplets enabled, skip to next start codon
        # while current sequence is not in alignment with first start codon
        while (strict_align
               and (seq_start_pos - first_seq_start_pos) % 3 != 0
               and seq_start_pos != -1):
            seq_start_pos = dna.find(START_SEQ, seq_start_pos + 1)

    if len(valid_genes) == 0:
        return [NO_GENES_MSG]
    return valid_genes
