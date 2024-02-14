# Modules import
from succinct_multiple_alignment import SuccinctMultipleAlignment
from succinct_column import SuccinctColumn
import argparse
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import psutil


# Functions definition
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
    plt.savefig("img.png")


# Main program
if __name__ == "__main__":
    # Arguments definition.
    parser = argparse.ArgumentParser(description="Compare a file of sequencing reads to a reference genome.")
    parser.add_argument("--file", "-f", help="Le fichier multifasta contenant l'alignement multiple."
                                             "ou bien le fichier de sauvegarde de l'alignement multiple succinct",
                        type=str, required=True, action='store')
    parser.add_argument("--ncols", "-n", help="Le nombre de colonnes lues a chaque lecture du fichier.",
                        type=int, default=1000, action='store')
    parser.add_argument("--compressed", "-c", help="Le fichier est compresse.", action='store_true')
    parser.add_argument("--infos", "-i", help="Ecrit des informations generales sur l'alignement multiple dans le terminal.",
                        action='store_true')
    parser.add_argument("--performance", "-p", help="Affiche les performances du programme en terme de temps et de memoire.",
                        action='store_true')
    parser.add_argument("--save", "-s", help="Stocke l'alignement multiple succinct dans un dossier compresse.",
                        action='store_true')
    parser.add_argument('--save_dir', '-sd', help="Chemin vers le repertoire ou sera strocker l'alignement multiple succinct",
                        type=str, action='store')
    parser.add_argument("--load", "-l", help="Charge l'alignement multiple succinct depuis un fichier de sauvegarde",
                        action='store_true')
    args = parser.parse_args()

    # Assert that the file exists.
    if not os.path.isfile(args.file):
        raise FileNotFoundError(args.file)

    align = SuccinctMultipleAlignment(args.file, nb_columns=args.ncols, compressed=args.compressed, load_file=args.load)
    
    if args.infos:
        print("SuccinctMultipleAlignment :")
        print('infos :')
        print(align.get_info())
        print('size_in_bytes():')
        print(align.size_in_bytes())
    
    if args.performance:
        plot_time_memory_psutil(align.fetch_column, args.file, 0, args.ncols)

    if args.save:
        align.store_to_file(args.save_dir)