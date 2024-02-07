"""

"""
#-*- coding: utf-8 -*-
# Modules import
import os
from succinct_column import SuccinctColumn
from Bio import SeqIO
import pysdsl
import gzip
import csv 

# Class definition
class SuccinctMultipleAlignment:

    def __init__(self, fasta_file, nb_columns=1000, vector="SDVector", compressed=False):
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
        if compressed:
            self.__size, self.__length = self.fetch_alignment_size_compress(fasta_file)
            for position in range(0, self.__length, nb_columns):
                self.__multialign += self.fetch_column_compress(fasta_file, position, nb_columns, vector)
        else:
            self.__size, self.__length = self.fetch_alignment_size(fasta_file)
            for position in range(0, self.__length, nb_columns):
                self.__multialign += self.fetch_column(fasta_file, position, nb_columns, vector)

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
        if not os.path.isfile(fasta_file):
            raise FileNotFoundError(fasta_file)

        with open(fasta_file, "r") as handle:
            seq_count = 0  # sequence counter
            align_length = None
            for record in SeqIO.parse(handle, 'fasta'):
                if align_length is None:
                    align_length = len(record.seq)
                    seq_count += 1
                else:
                    if align_length != len(record.seq):
                        raise ValueError
                    seq_count += 1
        return (seq_count, align_length)

    @staticmethod
    def fetch_alignment_size_compress(fasta_file):
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
        if not os.path.isfile(fasta_file):
            raise FileNotFoundError(fasta_file)

        with gzip.open(fasta_file, 'rt') as handle:
            seq_count = 0  # sequence counter
            align_length = None
            for record in SeqIO.parse(handle, 'fasta'):
                if align_length is None:
                    align_length = len(record.seq)
                    seq_count += 1
                else:
                    if align_length != len(record.seq):
                        raise ValueError
                    seq_count += 1
        return (seq_count, align_length)
    
    def fetch_column(self, fasta_file, position, nb_column, vector="SDVector"):
        """
        Read the FASTA file and store 'nb_column' columns as SuccinctColumn objects in a list.
        To do that, it reads 'nb_column' nucleotides in each sequence, starting at the position 'position' in the file 'fasta_file'.

        Parameters:
        -----------
        fasta_file : str
            A FASTA file containing multiple sequences aligned.
        position : int
            The position in the sequence where the search starts.
        nb_column : int
            The number of columns to build in a single run.

        Return:
        -------
        list(SuccinctColumn) :
            A list of SuccinctColumns objects of size 'nb_column' corresponsing to columns of the multialignment starting at the 
            position 'position'.
        """
        seq_count = 0
        nt_kept, previous_nt = [''] * nb_column, [''] * nb_column
        bit_vectors = []
        with open(fasta_file, "r") as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                i = 0
                while i < nb_column and position + i < self.__length:
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
            sd_vector.append(SuccinctColumn(bit_vectors[i], nt_kept[i], vector=vector))
        del bit_vectors, nt_kept, previous_nt
        return sd_vector

    def fetch_column_compress(self, fasta_file, position, nb_column, vector="SDVector"):
        """
        Read the FASTA file and store 'nb_column' columns as SuccinctColumn objects in a list.
        To do that, it reads 'nb_column' nucleotides in each sequence, starting at the position 'position' in the file 'fasta_file'.

        Parameters:
        -----------
        fasta_file : str
            A FASTA file containing multiple sequences aligned.
        position : int
            The position in the sequence where the search starts.
        nb_column : int
            The number of columns to build in a single run.

        Return:
        -------
        list(SuccinctColumn) :
            A list of SuccinctColumns objects of size 'nb_column' corresponsing to columns of the multialignment starting at the
            position 'position'.
        """
        seq_count = 0
        nt_kept, previous_nt = [''] * nb_column, [''] * nb_column
        bit_vectors = []
        with gzip.open(fasta_file, 'rt') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                i = 0
                while i < nb_column and position + i < self.__length:
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
            sd_vector.append(SuccinctColumn(bit_vectors[i], nt_kept[i], vector=vector))
        del bit_vectors, nt_kept, previous_nt
        return sd_vector

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
        """
        Returns the SDVector object of a SuccinctColumn object corresponding to the 'index'-th column in the multiple alignment.

        Parameters:
        -----------
        index : int
            The column considered.

        Return:
        -------
        pysdsl.SDVector :
            The SDVector object corresponding to the compacted representation of the bit vector of the column specified.
        """
        return self.__multialign[index].get_vector()

    def get_kept_nucleotide(self, index):
        """
        Returns the nucleotides - used to deduce the column's sequence from the bit vector - of the 'index'-th column.

        Parameters:
        -----------
        index : int
            The column consideres.

        Return:
        -------
        str :
            The nucleotides kept in the SuccinctColumn object corresponding to the 'index'-th column in the multiple alignment.
        """
        return self.__multialign[index].get_kept_nucleotide()

    def get_info(self):
        """
        Return general informations such as the alignment length (length of the sequences) and the alignment size (number of sequences).

        Parameters:
        -----------
        None

        Return:
        -------
        int :
            The alignment length (length of the sequences).
        int :
            The alignment size (number of sequences).
        """
        return self.__length, self.__size
    
   def size_to_csv(self, file_name="size.csv", sort_by_size=True):
        """
        Save the size in bytes of each SuccinctColumn object in a CSV file.

        Parameters:
        -----------
        file_name: str, optional
            The name of the CSV file to save the sizes. Default is "size.csv".
        sort_by_size: bool, optional
            If True, the sizes will be sorted in ascending order. Default is True.

        Return:
        -------
        None
        """
        ### la liste de tailles des colonnes 
        sizes = [(i, self.__multialign[i].size_in_bytes()) for i in range(self.__length)]
        ### triee les colonnes par ordre croissant  de taille 
        if sort_by_size:
            sizes.sort(key=lambda x: x[1])
        # ecriture dans le fichier CSV
        with open(file_name, "w") as fileOut:
            # 
            writer = csv.writer(fileOut)
            # ecriture d'en-tetes du csv
            writer.writerow(["Index", "column sorted by size","cumulative column sizes "])
            cumulative_size = 0
            for i, size in sizes:
                ###### pour la partie qui cumule  les tailles des colonnes
                cumulative_size += size
                writer.writerow([i, size, cumulative_size])
