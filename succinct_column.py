# Modules import
#-*- coding: utf-8 -*-
import pysdsl
import os
from collections import defaultdict


# Class definition
class SuccinctColumn:

    def __init__(self, bitvector=None, nt_kept=None, vector="SDVector", load=False, dir_path=None, column=0):
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
        if load:
            self.__vector, self.__nucleotides = self.load_from_file(dir_path=dir_path, column_nb=column)
        else:
            self.__nucleotides = nt_kept
            if vector == "BitVector":
                self.__vector = bitvector
            elif vector == "SDVector":
                self.__vector = pysdsl.SDVector(bitvector)

    def __len__(self):
        return len(self.__vector)

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
        
        return [(nt, round(count / float(length_vector), decimals)) for nt, count in nt_count_dict.items()]
    
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

    def store_to_file(self, column_number, output_dir):
        """
        Store the SDVector and the nucleotides in two files.
        Do not use if the bit vector is represented by a pysdsl.BitVector.

        Parameters:
        -----------
        column_number : int
            The index (in python) of the current column in the sequence
        project_name : str
            The name of the directory where the files will be created (Originally designed to store several Succinct_columns in one directory)
        output_dir : str
            The path / the directory where the directory 'project_name' will be created

        Return:
        -------
        None
        """
        self.__vector.store_to_file('{}/{}_column'.format(output_dir, column_number))
        with open('{}/{}.txt'.format(output_dir, column_number), 'w') as fileIn:
            fileIn.write(self.__nucleotides)

    def load_from_file(self, dir_path, column_nb):
        """
        Create a Succinct_column from the files produced by the store_to_file() function.

        Parameters:
        -----------
        dir_path : str
            The path / the save from which to recreate the saved Succinct_column.
        column_nb : int
            The index (in python) of the current column in the sequence

        Return:
        -------
        pysdsl.SDVector :
            The SDVector object corresponding to the compacted representation of the bit vector of the column specified.
        nucleotides : str
            Nucleotides corresponding to the '1' in the bit vector.
        """
        vector = pysdsl.SDVector.load_from_file('{}/{}_column'.format(dir_path, column_nb))
        with open('{}/{}.txt'.format(dir_path, column_nb)) as fileIn:
            nucleotides = fileIn.readline()
        return vector, nucleotides
