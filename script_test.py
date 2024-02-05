from succinct_multiple_alignment import SuccinctMultipleAlignment
from succinct_column import SuccinctColumn
import time
import matplotlib as plt
import os
import psutil

align = SuccinctMultipleAlignment('1000_seq.fasta', nb_columns=1000)

print("SuccinctMultipleAlignment :")
print('len() :')
print(len(align))
print('size_in_bytes() :')
print(align.size_in_bytes())

# print('fetch_column() :')
# fetch_column = align.fetch_column_V2('1000_seq.fasta', 0, 1000)
# print(fetch_column[0].get_vector())

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


# print('\nSuccinctColumn :')
#
# column = fetch_column[0]
#
# print('size_in_bytes() :')
# print(column.size_in_bytes())
# 
# print('len() :')
# print(len(column))
#
# print('get_nt() :')
# print(column.get_nt(3))

# print('get_vector() :')
# print(column.get_vector())

# print('nt_frequency() :')
# print(column.nt_frequency())

def time_profile(func, *args, **kwargs):
    #Enregistre le temps au debut avant l'execution de la fonction
    start_time = time.time()
    ### appel de la fonction
    result = func(*args, **kwargs)
    ### temps apres l'execution
    end_time = time.time()
    # calcul de la duree d'execution
    execution_time = end_time - start_time
    
    return result, execution_time

def memory_profile_psutil(func, *args, **kwargs):
    """ 
    cette fonction mesure le temps, le temps d'execution d'une fonction specifiee en utilisant 
    le module 
    """
    #processus
    process = psutil.Process(os.getpid())
    # print(process)
    # mesure l'utilisation de la memoire avant en convertissant les octets en mega-octets
    memory_avant = process.memory_info().rss / (1024.0 ** 2)  # En mega-octets car 1 mega-octet == 1024 kiloctet
    result = func(*args, **kwargs)
    ##1 mega-octet == 1024 kiloctet
    memory_apres = process.memory_info().rss / (1024.0 ** 2)
    ## mesure de la memoire apres 
    memory_usage_result = memory_avant - memory_apres
    return result, memory_usage_result

def plot_time_memory_psutil(func, *args, **kwargs):
    # Temps
    _, execution_time = time_profile(func, *args, **kwargs)
    
    # Memoire
    memory_usage_result = memory_profile_psutil(func, *args, **kwargs)
    
    # Generation des courbes
    plt.figure(figsize=(10, 5))
    #courbe du temps d'execution 
    plt.subplot(1, 2, 1)
    plt.bar(['Execution Time'], [execution_time], color='blue')
    plt.ylabel('Time (s)')
    plt.title('Execution Time')
    # courbe de la memoire
    plt.subplot(1, 2, 2)
    plt.bar(['Memory Usage'], [memory_usage_result[1]], color='green')
    plt.ylabel('Memory Usage (MiB)')
    plt.title('Memory Usage')

    plt.tight_layout()
    plt.show()

# Exemple d'utilisation :
plot_time_memory_psutil(align.fetch_column, '1000_seq.fasta', 0, 1000)
