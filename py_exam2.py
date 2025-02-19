from collections.abc import Sequence
from typing import Tuple

def generate_number_list(last_number: int = 21) -> Sequence:
    """Question 1
        Return a sequence of numbers, starting at 3, up to (and including, if appropriate), the argument last_number
        
        Example use: generate_number_list()
        Example output: 3 6 9 12 15 18 21
        Example use: generate_number_list(25)
        Example output: 3 6 9 12 15 18 21 24
    """
    # Enter the last number you want to generate a list up to
    
    last_number = input("Enter a number: ")
    last_number = int(last_number) + 1
    numbers = []
    for i in range(3, last_number, 3):
        numbers.append(i)
        return numbers




def lex_sort_file(filename: str = "multi_seqs.txt") -> Sequence:
    """Question 2
        The file with the filename given as an argument contains several amino-acid sequences, one per line.
        Write a function that returns lines in a sequence, sorted in lexicographical order.

        Example use: lex_sort_file()
        Example output:
        ['AALKIDSTVSQDSAWYTATAINKAGRDTTRCKVNVEVEFAEPEPERKLIIPRGTYRAK',
         'AEKTAVTKVVVAADKAKEQELKSRTKEVITTKQEQMHVTHEQIRKETEKTFVPKVV',
         'EAVATGAKEVKQDADKSAAVATVVAAVDMARVREPVISAVEQTAQRTTTTAVHIQPAQEQVRKE',
         'EGRKGLQRIEELERMAHEGALTGVTTDQKEKQKPDIVLYPEPVRVLEGETARF',
         ...
    """
    # Complete the function body below to answer question 2
    
    return

def top_lysine_stats(filename: str = "multi_seqs.txt") -> Tuple[float, str]:
    """Question 3
        Work out which sequence from file `multi_seqs.txt` has the highest percentage of lysine (`K`) residues,
        and return out both the percentage and the sequence.

        Example use: top_lysine_stats()
        Example output:  17.86, AEKTAVTKVVVAADKAKEQELKSRTKEVITTKQEQMHVTHEQIRKETEKTFVPKVV

    """
    # Complete the function body below to answer question 3
    
    return 

def avg_lysine_stats(filename="multi_seqs.txt") -> Tuple[float, float]:
    """Question 4
        Compute the mean and median number of lysines in sequences in `multi_seqs.txt`, and return them

        Example use: avg_lysine_stats()
        Example output:  (4.390243902439025, 4.0)
    """
    # Complete the function body below to answer question 4
    
    return

def plot_lysine_stats(filename: str = "multi_seqs.txt") -> None:
    """Question 5
        Wrte a function that plot the distribution of lysine counts, in the sequences from file `multi_seqs.txt`.

        Example use: plot_lysine_stats()
        Example output:  <plot of the lysine count distribution>
    """
    # Complete the function body below to answer question 5
    
    return

def translate_dna(codons_fname: str = '../data/codons.txt', dna_fname: str = '../data/dna.txt') -> Sequence[str]:
    """Question 6
        File `codons.txt` includes a non-standard codon table, including start and stop codons. Use the codons in that file to translate the DNA sequence in `dna.txt` into multiple protein sequences.
        To do this, you need to:
- Scan the DNA sequence (from its start), until you find a *start codon*
- Translate the codons one at a time until you encounter a *stop codon*
- Continue scanning the DNA sequence until you find the next *start codon*, and so on.

*Note 1:* Ignore the possibilty of different reading frames, i.e. assume the first codon consists of the first 3 letters in the sequence, and so on.

*Note 2:* Do not use Biopython to answer this question!

Here is a **simple example**:

DNA sequence: `CGTATGGGTTCGATGTCGGTCTAACCC` 

`CGT` &mdash; not a *start codon*, so skip it<br>
`ATG` &mdash; *start codon*<br>
`GGT` &mdash; translate to G<br>
`TCG` &mdash; translate to S<br>
`ATG` &mdash; translate to M (identical to *start codon*, but we have already started)<br>
`TCG` &mdash; translate to S<br>
`GTC` &mdash; translate to V<br>
`TAA` &mdash; *stop codon*<br>
`CCC` &mdash; not a *start codon*, so skip it<br>

So for this DNA sequence, the peptide sequence `GSMSV` should be returned.

        Example use: translate_dna()
        Example output:  ['YTSRRSPSSVGF', ...]
    """
    # Complete the function body below to answer question 6
    
    return

def longest_translatable_sequence(codons_fname: str = '../data/codons.txt', dna_fname: str = '../data/dna.txt') -> int:
    """Question 7
    Find and return the length of the longest translatable sequence in `dna.txt`:
    - the sequence starts with a start codon
    - the sequence ends with a stop codon
    - both start and stop codons should be in the same reading frame
    - there cannot be a stop codon in the same reading frame within the sequence
    - consider all three reading frames

    N.B.: Attempt this question last!

            Example use: longest_translatable_sequence()
            Example output:  245
    """
    # Complete the function body below to answer question 7

    return


if __name__ == '__main__':
    print(generate_number_list())
    print(lex_sort_file(filename='data/multi_seqs.txt'))
    print(top_lysine_stats(filename='data/multi_seqs.txt'))
    print(avg_lysine_stats(filename='data/multi_seqs.txt'))
    print(plot_lysine_stats(filename='data/multi_seqs.txt'))
    print(translate_dna(codons_fname='data/codons.txt', dna_fname='data/dna.txt'))
    print(longest_translatable_sequence(codons_fname='data/codons.txt', dna_fname='data/dna.txt'))

