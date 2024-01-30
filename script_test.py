from essai import SuccinctMultipleAlignment
from succinct_column import SuccinctColumn

align = SuccinctMultipleAlignment('1000_seq.fasta.gz', nb_columns=1000)

print("SuccinctMultipleAlignment :")
print('len() :')
print(len(align))
print('size_in_bytes() :')
print(align.size_in_bytes())

# print('fetch_column() :')
fetch_column = align.fetch_column_V2('1000_seq.fasta.gz', 0, 1000)
print(fetch_column[0].get_vector())

print('get_nt() :')
print(align.get_nt(3, 5))

# print('get_sequence() :')
# print(align.get_sequence(3))

# print('get_vector() :')
# print(align.get_vector(0))

print('get_kept_nucleotide() :')
print(align.get_kept_nucleotide(6))

print('get_info() :')
print(align.get_info())


print('\nSuccinctColumn :')
#
column = fetch_column[0]
#
print('size_in_bytes() :')
print(column.size_in_bytes())
# 
# print('len() :')
# print(len(column))
#
# print('get_nt() :')
# print(column.get_nt(3))

# print('get_vector() :')
# print(column.get_vector())

print('nt_frequency() :')
print(column.nt_frequency())
