"""O3.2: Unit testing of find genes.

+----------------------------+
| DTE-2510 FA24 (UiT)        |
|                            |
| Submission:   O3.2         |
| Title:        Unit testing |
| Student:      Bjarte Lode  |
+----------------------------+

Description:
------------
Unit tests for the find_genes function. The test cases are based on the
requirements from the original assignment in O2. I ended up rewriting parts of
my original solution, because I found a couple of bugs and wasn't happy with
my implementation.

The assignment brief for O3, asked what changes needed to be made to the
original code for it to pass the test cases. First, I had to change the return
value of the function to return a list of strings instead of a list of lists of
strings. Second, in my original solution I thought I made a smart decision by
splitting the list first by the start codon ATG. This wasn't such a great idea
after all! Consider the sequence 'ATG CCC TA[A TG]T'. This should be
interpreted like my spaces indicate, but splitting first made the ATG marked
with square bracets divide the list.

Basically this called fore a complete rewrite of the function. I also made it
able to ignore noisy characters before the start codon, and added a check for
invalid input (non-ACGT characters), to raise a ValueError. In the return
sequences I joined each triplet with ',' in my original solution, but I changed
this to join the triplets into a single string. If multiple genes are found,
they are returned as a list of strings, as indicated in the example test
cases from the assignment brief.

"""

import unittest
from find_genes import find_genes


class Test_genome(unittest.TestCase):
    """Unittest for find_genes.py

    Test cases:
    - Empty string
    - Original string from assignment
    - Only start and stop codons
    - Only start codon
    - Start codon with no stop codon
    - Stop codon with no start codon
    - ATG-TAG sequence
    - ATG-TAA sequence
    - ATG-TGA sequence
    - Adjacent start and stop codons
    - Multiple stop codons
    - Multiple start codons
    - Leading noise, letters before start codon
        not multiple of 3 (like CAATG)
    - Lowercase input
    - Leading and trailing whitespace input
    - Invalid input (non-ACTG characters)
    - Long sequence of genes

    """
    # Tests for invalid input

    def test_invalid_type_raises_error(self):
        """Test invalid input type raises TypeError."""
        self.assertRaises(TypeError, find_genes, None)

    def test_invalid_input_raises_error(self):
        """Test invalid input (any char not in ATCG) raises ValueError."""
        self.assertRaises(ValueError, find_genes,
                          "ATGCYBREQTAA")

    # Tests for handling input with non standard format

    def test_lcase_input_is_handled(self):
        """Test lowercase input is handled like upper case input."""
        result = find_genes("atgatttag")
        self.assertEqual(result[0], "ATT")

    def test_lead_trail_whitespace_input_is_handled(self):
        """Test leading and trailing whitespace is handled (stripped)."""
        result = find_genes("  ATGCTCTAA ")
        self.assertEqual(result[0], "CTC")

    # Tests strings with no valid genes

    def test_emty_string(self):
        """Test empty string as input (no valid genes)."""
        result = find_genes("")
        self.assertEqual(result[0], "No genes found")

    def test_only_start_stop_codons_returns_no_genes_found(self):
        """Test only start and stop codons in sequence (no valid genes)."""
        result = find_genes("ATGATGATGTAATAGTGAATGTGATAA")
        self.assertEqual(result[0], "No genes found")

    def test_only_start_codon_returns_no_gene_found(self):
        """Test only start codon in sequence (no valid genes)."""
        result = find_genes("ATGATGATGATG")
        self.assertEqual(result[0], "No genes found")

    def test_start_no_stop_codon_returns_no_gene_found(self):
        """Test start codon with no stop codon (no valid genes)."""
        result = find_genes("TTCATGTTTTTAGGACGGGGCGTTACT")
        self.assertEqual(result[0], "No genes found")

    def test_stop_no_start_codon_returns_no_gene_found(self):
        """Test stop codon with no start codon (no valid genes)."""
        result = find_genes("ATCACAAGGAGGCTTTAGGTAACG")
        self.assertEqual(result[0], "No genes found")

    # Tests for start and stop codons

    def test_atg_tag_seq(self):
        """Test ATG-TAG sequence."""
        result = find_genes("ATGGGGCGCTAGCAA")
        self.assertEqual(result[0], "GGGCGC")

    def test_atg_taa_seq(self):
        """Test ATG-TAA sequence."""
        result = find_genes("ATGATCGTGCAATAAAGC")
        self.assertEqual(result[0], "ATCGTGCAA")

    def test_atg_tga_seq(self):
        """Test ATG-TGA sequence."""
        result = find_genes("ATGACTATTTGATGATGC")
        self.assertEqual(result[0], "ACTATT")

    def test_adjacent_start_stop(self):
        """Test sequence with adjacent start and stop codons."""
        result = find_genes("GTGATGTAAATGGATTAG")
        self.assertEqual(result[0], "GAT")

    def test_multi_stop(self):
        """Test sequence with multiple stop codons (one valid gene)."""
        result = find_genes("ATGACGATGCTGTAGACGTAA")
        self.assertEqual(result[0], "CTG")

    def test_multi_start(self):
        """Test sequence with multiple start codons (one valid gene)."""
        result = find_genes("ATGACATGGATGCTTAGGTAACGC")
        self.assertEqual(result[0], "CTTAGG")

    def test_double_start(self):
        """Test sequence with double start codon (one valid gene)."""
        result = find_genes("ATGATGACATTGTAGCTT")
        self.assertEqual(result[0], "ACATTG")

    def test_start_not_triplet_pos(self):
        """Test sequence with start codon not at triplet position."""
        result = find_genes("CAATGACATGATGCTGC")
        self.assertEqual(result[0], "ACA")

    # Tests for longer sequences and strict / loose alignment

    def test_start_not_triplet_pos_in_middle(self):
        """Test strict and loose alignment."""
        seq = "ATGACATGATGATGCTGTAG"
        self.assertEqual(find_genes(seq, strict_align=False),
                         ["ACA", "CTG"])
        self.assertEqual(find_genes(seq, strict_align=True),
                         ["ACA"])

    def test_orig_string(self):
        """Test original string from assignment."""
        result = find_genes("TTATGTTTTAAGGATGGGGCGTTAGTT")
        self.assertEqual(result, ['TTT', 'GGGCGT'])

    def test_long_seq_strict_loose_alignment(self):
        """Test long sequence of genes using strict and loose alignment."""
        seq = "GTGATGTAAATGGATTAGATGACATGGATGCTTAGGTAACGC"
        seq += "ATGATCGTGCAATAAGCTTCATGTTTTAAATGGGGCGTTAAGTT"
        self.assertEqual(
            find_genes(seq, strict_align=False),
            ["GAT", "CTTAGG", "ATCGTGCAA", "TTT", "GGGCGT"]
        )
        self.assertEqual(
            find_genes(seq, strict_align=True),
            ["GAT", "CTTAGG", "ATCGTGCAA"]
        )


# Run the tests
if __name__ == '__main__':
    unittest.main()
