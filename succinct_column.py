# Modules import
#-*- coding: utf-8 -*-
import pysdsl
from collections import defaultdict 

# Class definition
class SuccinctColumn:

    def __init__(self, bitvector, nt_kept, vector="SDVector"):
        """
        Build a SDVector or a BitVector and a sequence of nucleotides (corresponding to the "1" in the bit sequence) from all the 
        nucleotides in a column.

        Parameters:
        -----------
        bitvector : pysdsl.BitVector
            A bit vector corresponding to a simplified version of multiple alignment.
        nt_kept : str
            Nucleotides corresponding to the '1' in the bit vector.
        vector : str
            Selection of the class representing the bit vector.

        Return:
        -------
        None
        """
        self.__nucleotides = nt_kept
        if vector == "BitVector":
            self.__vector = bitvector
        elif vector == "SDVector":
            self.__vector = pysdsl.SDVector(bitvector)
        # self.__nucleotides = nt_kept
        # if column_seq == "BitVector":
        #     self.__vector.append(column_seq)
        # elif vector == "SDVector":
        #     self.__vector = [pysdsl.SDVector(bit_vector) for bit_vector in column_seq]
        # print(type(self.__vector))
        # print(self.__vector)

    def __len__(self):
        return len(self.__vector)

    @staticmethod
    def seq_to_bitvector_nts(column_seq):
        """
        Builds a sequence of 0 and 1 and a string of nucleotides corresponding to the "1" from the nucleotides in a column. To clarify,
        if in the column a nucleotide if the same as the previous one, the value at the considered position is 0, else it is 1.

        Parameters:
        -----------
        column_seq : str
            All the nucleotides in a column of a multiple alignment.
        
        Return:
        -------
        list(int) :
            List of 0 and 1 representing a simplified version of the column.
        str : 
            The nucleotides corresponding to each "1" in the simplified version of the column.
        """
        bits_list = []
        nt_kept = ""
        previous_nt = ""
        for nt in column_seq:
            if nt != previous_nt:
                bits_list.append(1)
                nt_kept += nt
                previous_nt = nt
            else:
                bits_list.append(0)
        return bits_list, nt_kept

    def size_in_bytes(self):
        """
        Return the size in bytes of the pysdsl vector representing the column of nucleotides.

        Parameters:
        -----------
        None

        Return:
        -------
        int : 
            The size in bytes of the pysdsl vector representing the column of nucleotides
        """
        return self.__vector.size_in_bytes + len(self.__nucleotides)

   def nt_frequency(self, decimals=2):
        """
        Returns the percentage of each nucleotide in the column.

        Parameters:
        -----------
        None

        Return:
        -------
        (float, float, float, float, float) :
            - percentage of A
            - percentage of T
            - percentage of C
            - percentage of G
            - percentage of -
            -percentage of other nucleotide present
        """
        count_one = -1
        nt_count_dict = defaultdict(int)
        
        for byte in self.__vector:
            if byte == 1:
                count_one += 1
            nt_count_dict[self.__nucleotides[count_one]] += 1
            
        length_vector = len(self.__vector)
        print(nt_count_dict)
        return tuple(round(nt_count_dict[char]/float(length_vector), decimals) for char in nt_count_dict)
    
    def get_nt(self, position):
        """ 
        Returns the nucleotide at the position specified in the column (the p-th sequence in the alignment).

        Parameters:
        -----------
        position : int
            The position of the nucleotide in the column.
        
        Return:
        -------
        str :
            The target nucleotide.
        """
        length_vector = len(self.__vector)
        if position == length_vector - 1:
            nt = self.__nucleotides[-1]
        else:
            to_rank = self.__vector.init_rank_1()
            nt = self.__nucleotides[to_rank.rank(position + 1) - 1]
            del to_rank
        return nt

    def get_vector(self):
        """
        Returns the SDVector object corresponding to the compacted representation of the bit vector.

        Parameters:
        -----------
        None

        Return:
        -------
        pysdsl.SDVector :
            The SDVector object corresponding to the compacted representation of the bit vector.
        """
        return self.__vector

    def get_kept_nucleotide(self):
        """
        Returns the nucleotides used to deduce the column's sequence from the bit vector.

        Parameters:
        -----------
        None

        Return:
        -------
        str :
            The nucleotides kept.
        """
        return self.__nucleotides
