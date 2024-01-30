#-*- coding: utf-8 -*-
# Modules import
import os
from succinct_column import SuccinctColumn
from Bio import SeqIO
import pysdsl
import gzip


class SuccinctMultipleAlignment:

    def validate_fasta_file(self, fasta_file):
        """
        Valide le fichier FASTA en utilisant Biopython.

        Parameters:
        -----------
        fasta_file : str
            Le chemin du fichier FASTA.

        Raises:
        -------
        FileNotFoundError:
            Si le fichier FASTA n'existe pas.
        ValueError:
            Si le fichier FASTA n'est pas au format attendu.
        """
        if not os.path.isfile(fasta_file):
            raise FileNotFoundError(fasta_file)

        # Essaye de lire le fichier avec SeqIO, cela va lever une exception si le format est incorrect.
        try:
            SeqIO.read(fasta_file, "fasta")
        except Exception as e:
            raise ValueError("Le fichier FASTA n'est pas au format attendu.")

    def __init__(self, fasta_file, nb_columns = 100, vector="SDVector"):
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
        #self.validate_fasta_file(fasta_file)  # Valide le fichier FASTA.
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
       
            
        with gzip.open(fasta_file, 'rt') as fileIn:
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
        with gzip.open(fasta_file, 'rt') as fileIn:
            column_seq = ""
            for line in fileIn.readlines():
                if line.startswith(">"):  # New sequence, reset variables
                    len_still_to_read = position
                elif len_still_to_read != -1:  # If the nucleotide of the sequence have already been stored.
                    len_seq = len(line.strip())
                    if len_still_to_read >= len_seq:  # Instead of rebuilding the sequence, we only store the number of nt still to read.
                        len_still_to_read -= len_seq
                    else:
                        nt= line[len_still_to_read].upper()
                        if nt not in {'A','T','C','G','-'}:
                            raise ValueError("There is an invalid nucleotide in the sequence.")
                        column_seq += nt
                        len_still_to_read = -1
        return column_seq
    
    def fetch_column_V2(self, fasta_file, position, nb_column, vector="SDVector"):
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
        n_kept, previous_nt = ['']*nb_column, ['']*nb_column
        bit_vectors = []
        with gzip.open(fasta_file, 'rt') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                i = 0
                sequence=record.seq
                sequence_length=len(sequence)
                # Check if the sequence is long enough for the specified position and number of columns
               
                if position + nb_column > sequence_length:
                    #print(" too short for the specified position and number of columns")
                    continue 
                while i < nb_column:
                    # record.seq[position + i] is the current nucleotide
                    if seq_count == 0:
                        bit_vectors.append(pysdsl.BitVector(self.__size))
                        bit_vectors[i][seq_count] = 1
                        n_kept[i] += record.seq[position + i].upper()
                        previous_nt[i] = record.seq[position + i]
                        if len(bit_vectors) == 0:
                        # Gérer le cas où bit_vectors est vide (peut-être lever une exception ou ajuster le comportement selon votre logique)
                            raise ValueError("La liste bit_vectors est vide.")
                    elif seq_count != 0 and previous_nt[i] != record.seq[position + i]:
                        bit_vectors[i][seq_count] = 1
                        n_kept[i] += record.seq[position + i].upper()
                        previous_nt[i] = record.seq[position + i]
                    i += 1
                seq_count += 1
         
        #print("Longueur de bit_vectors :", len(bit_vectors))
        #print("Index utilisé :", i)
    
        return [SuccinctColumn(bit_vectors[i], n_kept[i], vector=vector) for i in range(min(nb_column, len(bit_vectors)))]

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
