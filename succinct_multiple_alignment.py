"""

"""

# Modules import
import os
from succinct_column import SuccinctColumn
from Bio import SeqIO
import pysdsl

# Class definition
class SuccinctMultipleAlignment:

    def __init__(self, fasta_file, nb_columns = 1000, vector="SDVector"):
        """
        Build the succinct multiple alignment as a list of objects SuccinctColumn.

        Parameters:
        -----------
        fasta_file : str
            A FASTA file containing multiple sequences aligned.
        vector : str
            Selection of the class representing the bit vector.

        Return:
        -------
        None
        """
        self.__multialign = []
        self.__size, self.__length = self.fetch_alignment_size(fasta_file)
        for position in range(0, self.__length, nb_columns):
            self.__multialign += self.fetch_column_V2(fasta_file, position, nb_columns, vector)
    def __len__(self):
        return len(self.__multialign)

    @staticmethod
    def fetch_alignment_size(fasta_file):
        """
        Read the FASTA file and store the size and the number of sequences.

        Parameters:
        -----------
        fasta_file : str
            A FASTA file containing multiple sequences aligned.
        
        Return:
        -------
        (int, int):
            - The number of sequences.
            - The length of the sequences (which is supposed to be the same for every sequence).
        """
        with open(fasta_file, "r") as fileIn:
            seq_count = 0  # sequence counter
            for line in fileIn.readlines():
                if line.startswith(">"):
                    if seq_count == 1:
                        align_length = len_seq  # The length of the 1rst sequence is the length of the alignment.
                    seq_count += 1
                    len_seq = 0  # Length of the sequence read.
                else:
                    len_seq += len(line.strip())
        return (seq_count, align_length)

    @staticmethod
    def fetch_column(fasta_file, position):
        """
        Read the FASTA file and store the nucleotide at a specified position in each sequence.

        Parameters:
        -----------
        fasta_file : str
            A FASTA file containing multiple sequences aligned.
        position : int
            The position of the nucleotide to store in each sequence.
        
        Return:
        -------
        str :
            The nucleotides of every sequences at a specified position.
        """
        with open (fasta_file, "r") as fileIn:
            column_seq = ""
            for line in fileIn.readlines():
                if line.startswith(">"):  # New sequence, reset variables
                    len_still_to_read = position
                elif len_still_to_read != -1:  # If the nucleotide of the sequence have already been stored.
                    len_seq = len(line.strip())
                    if len_still_to_read >= len_seq:  # Instead of rebuilding the sequence, we only store the number of nt still to read.
                        len_still_to_read -= len_seq
                    else:
                        column_seq += line[len_still_to_read].upper()
                        len_still_to_read = -1
        return column_seq

    def size_in_bytes(self):
        """
        Return the size in bytes of the entire succinct multiple alignment (sum of the size in bytes of all the SuccinctColumn objects).

        Parameters:
        -----------
        None

        Return:
        -------
        int :
            The size in bytes of the entire succinct multiple alignment.
        """
        return sum([succinct_column.size_in_bytes() for succinct_column in self.__multialign])

    def get_nt(self, seq_index, position):
        """
        Return the nucleotide in the position specified in the sequence of index "seq_index".

        Parameters:
        -----------
        seq_index : int
            The index of the sequence to search in.
        position : int
            The position to look at in the sequence.

        Return:
        -------
        str : The nucleotide in the position specified in the sequence of index "seq_index".
        """
        return self.__multialign[position].get_nt(seq_index)

    def get_sequence(self, seq_index):
        """
        Return the sequence of index "seq_index".

        Parameters:
        ----------
        seq_index : int
            The index of the sequence to search in.
        
        Return:
        -------
        str :
            The sequence of index "seq_index".
        """
        return "".join([self.get_nt(seq_index, position) for position in range(self.__length)])

    def get_vector(self, index):
        return self.__multialign[index].get_vector()

    def get_kept_nucleotide(self, index):
        return self.__multialign[index].get_kept_nucleotide()

    def get_info(self):
        return self.__length, self.__size

    def fetch_column_V2(self, fasta_file, position, nb_column, vector="SDVector"):
        """
        Read the FASTA file and store the nucleotide at a specified position in each sequence.

        Parameters:
        -----------
        fasta_file : str
            A FASTA file containing multiple sequences aligned.
        position : int
            The position of the nucleotide to store in each sequence.

        Return:
        -------

        """
        seq_count = 0
        nt_kept, previous_nt = ['']*nb_column, ['']*nb_column
        bit_vectors =  []
        with open(fasta_file, "r") as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                i = 0
                #for i in range(0, nb_column):
                while i < nb_column and position+i < self.__length:
                    if seq_count == 0:
                        bit_vectors.append(pysdsl.BitVector(self.__size))
                        bit_vectors[i][seq_count] = 1
                        nt_kept[i] += record.seq[position + i].upper()
                        previous_nt[i] = record.seq[position + i]
                    elif seq_count != 0 and previous_nt[i] != record.seq[position + i]:
                        bit_vectors[i][seq_count] = 1
                        nt_kept[i] += record.seq[position + i].upper()
                        previous_nt[i] = record.seq[position + i]
                    i += 1
                seq_count += 1
        sd_vector = []
        for i in range(len(bit_vectors)):
            sd_vector.append(SuccinctColumn((bit_vectors[i], nt_kept[i]), vector=vector))

        return sd_vector
