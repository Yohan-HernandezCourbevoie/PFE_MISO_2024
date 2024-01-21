from succinct_multiple_alignment import SuccinctMultipleAlignment
from succinct_column import SuccinctColumn

align = SuccinctMultipleAlignment('6_seq.fasta')
print("SuccinctMultipleAlignment :")
print('size_in_bytes() :')
print(align.size_in_bytes())

print('fetch_column() :')
fetch_column = align.fetch_column('6_seq.fasta', 3)
print(fetch_column)

print('get_nt() :')
print(align.get_nt(3, 4858))

print('get_sequence() :')
print(align.get_sequence(3))

print('get_vector() :')
print(align.get_vector(1))

print('get_info() :')
print(align.get_info())

print('\nSuccinctColumn :')

column = SuccinctColumn(fetch_column)

print('SuccinctColumn :')
print('seq_to_bytes_nts :')
print(column.seq_to_bytes_nts(fetch_column))

print('size_in_bytes() :')
print(column.size_in_bytes())

print('nt_frequency() :')
print(column.nt_frequency())

print('get_nt() :')
print(column.get_nt(3))

print('get_vector() :')
print(column.get_vector())