from collections.abc import Sequence
from typing import Tuple, Sequence, List
import numpy as np
from statistics import median
import matplotlib.pyplot as py

def generate_number_list(last_number: int = 21) -> Sequence:
    """
        Returns a sequence of numbers, starting at 3, up to (and including, if appropriate), the argument last_number
        
        Example use: generate_number_list(29)
        Example output: [3, 6, 9, 12, 15, 18, 21, 24, 27]
    """
    # Enter the last number you want to generate a list up to
    numbers = list(range(3, last_number+1, 3))
    return numbers


def lex_sort_file(filename: str = "multi_seqs.txt") -> Sequence:
    """
        Sorts the sequences in a file in lexicographical order and return them as a list.

        Example use: lex_sort_file()
        Example output:
        ['AALKIDSTVSQDSAWYTATAINKAGRDTTRCKVNVEVEFAEPEPERKLIIPRGTYRAK',
         'AEKTAVTKVVVAADKAKEQELKSRTKEVITTKQEQMHVTHEQIRKETEKTFVPKVV',
         'EAVATGAKEVKQDADKSAAVATVVAAVDMARVREPVISAVEQTAQRTTTTAVHIQPAQEQVRKE',
         'EGRKGLQRIEELERMAHEGALTGVTTDQKEKQKPDIVLYPEPVRVLEGETARF',
         ...
    """
    # Read the file line by line
    with open(filename) as f:
        lines = f.readlines()
        # strip spaces and sort the lines lexicographically
        lines = [line.strip() for line in lines]
        lines.sort()

    return "\n".join(lines)


def top_lysine_stats(filename: str = "multi_seqs.txt") -> tuple[float, str]:
    """Question 3
        Work out which sequence from file `multi_seqs.txt` has the highest percentage of lysine (`K`) residues,
        and return out both the percentage and the sequence.

        Example use: top_lysine_stats()
        Example output:  17.86, AEKTAVTKVVVAADKAKEQELKSRTKEVITTKQEQMHVTHEQIRKETEKTFVPKVV

    """
    # Open the file and read the lines and strip the spaces
    with open(filename) as f:
        # Return list of each line of the file as a string
        lines = f.readlines()
        lines = [line.strip() for line in lines]

        # Create a list to store the lysine percentages
        lysine_percentages = []
        for line in lines:
            # Count the number of lysine residues (K) in each sequence (line)
            lysine_count = line.count("K")
            # Calculate the percentage of lysine residues in each sequence
            lysine_percent = lysine_count / len(line) * 100
            lysine_percentages.append(lysine_percent)

        max_percent = max(lysine_percentages)
        
        average_lysine_percent = sum(lysine_percentages) / len(lysine_percentages)
        sequence = lines[lysine_percentages.index(max_percent)]

    return max_percent, sequence


def avg_lysine_stats(filename="multi_seqs.txt") -> Tuple[float, float]:
    """
        Calculate and return the mean and median number of lysines in sequences in a file.

        Example use: avg_lysine_stats()
        Example output:  (4.390243902439025, 4.0)

        Needed packages: numpy, statistics
    """
    # Complete the function body below to answer question 4

    with open(filename) as f:
        lines = f.readlines()
        lines = [line.strip() for line in lines]
        # Create a list for lysine counts
        lysine_counts = []
        for line in lines:
            lysine_count = line.count("K")
            # Add the counts to the list
            lysine_counts.append(lysine_count)
            # Calculate the mean lysine count
            if lysine_counts:
                mean_lysine = sum(lysine_counts) / len(lysine_counts)
                # Calculate the median lysine count 
                lysine_counts.sort()
                median_lysine = median(lysine_counts)
                return median_lysine, mean_lysine  
            # If there are no lysine counts, return 0
            else:
                mean_lysine = 0
                median_lysine = 0

        return 0,0
        
  

def plot_lysine_stats(filename: str = "multi_seqs.txt") -> None:
    """Question 5
        Write a function that plot the distribution of lysine counts, in the sequences from file `multi_seqs.txt`.

        Example use: plot_lysine_stats()
        Example output:  <plot of the lysine count distribution>
        
    """
    # Complete the function body below to answer question 5
    with open(filename) as f:
        lines = f.readlines()
        lines = [line.strip() for line in lines]
        # Create a list for lysine counts
        lysine_counts = []
        for line in lines:
            lysine_count = line.count("K")
            # Add the counts to the list
            lysine_counts.append(lysine_count)
       
        # Plot the distribution of lysine counts
        py.hist(lysine_counts, bins=20)
        py.xlabel("Lysine Count")
        py.ylabel("Frequency")
        py.title("Distribution of Lysine Counts")
    
    return py.show()

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

    # Read the codons file and store the codons in a dictionary
    with open(codons_fname) as f:
        lines = f.readlines()
        lines = [line.strip() for line in lines]
        codons = {}
        for line in lines:
            codon, aa = line.split()
            codons[codon] = aa
    # Define start and stop codons
    start_codon = 'ATG'
    stop_codons = ['TAA', 'TAG', 'TGA']
    # Read the DNA file and store the DNA sequence
    with open(dna_fname) as f:
        dna = f.read().strip()
    # Create a list to store the protein sequences
    protein_sequences = []
    # Iterate through the DNA sequence, finding the start codon

    


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

