
# README  PROJET DE FIN D'ETUDE MISO_2024
**Authors** : Arthur , Yohan et Awa

**Miso Fevrier 2024** 

## 1 . Description générale du projet 
Le projet " **Représentation succincte d'alignements multiples de SARS-CoV-2** "vise à développer une méthode efficace pour représenter de manière concise les alignements multiples de génomes du virus SARS-CoV-2.Cette initiative découle de la nécessité de gérer et d'analyser efficacement une grande quantité de données génomiques du SARS-CoV-2, collectées depuis le début de la pandémie.

L'objectif principal du projet est de concevoir une représentation compacte des alignements multiples qui permette de récupérer rapidement n'importe quel génome ou colonne spécifique tout en optimisant l'utilisation de la mémoire et en réduisant les temps d'accès. Pour atteindre cet objectif, l'idée est de représenter chaque colonne de l'alignement sous forme de vecteurs de bits, où un bit est défini à 1 si le nucléotide dans la colonne diffère du précédent, sinon il est à 0. Les données associées aux positions où les bits sont à 1 sont stockées dans un tableau séparé.

## 2. Configuration compatibles 
Les configurations compatibles avec le programme sont les suivantes 
Versions de Python compatibles :

Python 2.7  ou version ultérieure

Systémes d'exploitation:
Toutes les plateformes prenant en charge Python 2.7 (Windown, linux, macOS)


## 3. Installation des packages

#### Installation de Python27

installer miniforg: <https://github.com/conda-forge/miniforge#install>
Ensuite il faut créez un environnement pour Python 2.7:
```js
mamba create  -y -n py27 python=2.7
or 
conda  create  -y -n py27 python=2.7
```
Activation de l'environnement: 
```js 
 conda activate py27
 or 
 mamba activate py27
 ```
#### Installation de pysdsl

```js 
pip install https://github.com/QratorLabs/pysdsl/releases/download/1.0.0a0/pysdsl-1.0.0a0-cp27-cp27mu-linux_x86_64.whl
```
#### Installation de Biopython pour py27
```js
mamba install biopython=1.76
```

## 4. Fonctionnalités

### Classe 'SuccinctMultipleAlignement`:


#### 1. conctruction de l'alignement multiple succinct: 
Cette classe permet de construire un alignement multiple à partir d'un fichier FASTA contenant les Séquences alignées. On peut specifier le nombre de colonnes à considérer et le type de vecteur (SDVector ou BitVector) à utiliser pour utiliser pour la représentation.


#### 2. Extraction de colonnes specifiques

La methode 'fetch_column'  qui lit le fichier FASTA et stocker les colonnes 'nb_column' sous forme d'objets SuccinctColumn dans une liste
#### 3. Calcul de la taille en bytes de l'alignement 
La méthode 'size_in_bytes()'
renvoie la taille en octets de l'ensemble de l'alignement multiple succinct, ce qui donne une indication de l'espace mémoire occupé par l'alignement.

#### 4.Extraction de nucléotides specifique dans une colonne 
La méthode 'get_nt()'
permet d'extraire un nucléotide spécifique d'une colonne donnée de l'alignement multiple en spécifiant sa position.

####  5. Récupération de la sequence d'une colonne 
La méthode 'get_sequence()' retourner la séquence d'index "seq_index"
#### 6. Récupération du vecteur SDVector 
 La méthode 'get_vector()' permet d'obtenir l'objet SDVector correspondant à la représentation compacte du vecteur de bits d'une colonne donnée.

#### 7.Récupération des nucléotides conservés
 La méthode 'get_kept_nucleotide()' permet de récupérer les nucléotides conservés utilisés pour déduire la séquence de la colonne à partir du vecteur de bits.

#### 8. Sauvegarde en fichier CSV

La méthode 'size_to_csv()' retourne un fichier CSV contient trois colonnes : l'index de la colonne, la taille de la colonne triée par taille (si triée), et la taille cumulative des colonnes 

#### 9. La méthode 'store_to_file(self, output_dir)': 
Cette méthode stocke tous les objets SuccinctColumn de l'SuccinctMultipleAlignment dans un répertoire compressé. Elle prend en paramètre le répertoire de sortie où le fichier compressé sera créé. Les colonnes sont stockées individuellement dans le répertoire, puis le répertoire est compressé en un fichier tar.gz. Un fichier info.txt est également créé pour stocker des informations sur la taille et la longueur des séquences.

#### 10: La méthode 'load_from_file(self, compressed_save)': 
Cette méthode recrée un SuccinctMultipleAlignment à partir des fichiers produits par la méthode store_to_file(). Elle prend en paramètre le chemin du fichier compressé. Les fichiers sont extraits du fichier compressé, les informations sur la taille et la longueur des séquences sont récupérées à partir du fichier info.txt, puis les SuccinctColumn sont chargés à partir des fichiers individuels pour reconstruire l'alignement multiple.

#### 11. La méthode 'find_columns_with_excessive_space()': 
Il identifie les colonnes qui occupent significativement plus d'espace que la taille moyenne des colonnes.


### Classe 'SuccinctColumn' :
#### 1.Construction de l'objet : 
La classe permet de construire un objet SuccinctColumn à partir d'un vecteur de bits (SDVector ou BiTvector ) et des nucléotides correspondant aux "1" dans le vecteur.
#### 2.Calcul de la taille en bytes de l'objet:
La méthode 'size_in_bytes() 'fournit une fonction pour calculer la taille en bytes de l'objet SuccinctColumn, ce qui donne une indication de l'espace mémoire occupé par la colonne de l'alignement.

#### 3. Calcul de la fréquence des nucléotides dans la colonne :
La méthode 'nt_frequency()' permet de calculer la fréquence des différents nucléotides présents dans la colonne, ce qui permet d'évaluer la diversité nucléotidique.
#### 4. Extraction de nucléotides spécifiques dans la colonne :
 La méthode 'get_nt()' permet d'extraire un nucléotide spécifique d'une colonne donnée de l'alignement multiple en spécifiant sa position.

#### 5. Récupération du vecteur SDVector :
La méthode 'get_vector()' permet d'obtenir l'objet SDVector correspondant à la représentation compacte du vecteur de bits de la colonne.

#### 6. Récupération des nucléotides conservés :
La méthode 'get_kept_nucleotide()'renvoie les nucléotides utilisés pour déduire la séquence de la colonne à partir du vecteur de bits.elle initialise un objet init_rank pour accéder aux rangs des 1 dans le vecteur de bits.

#### 7. La méthode 'store_to_file()'
Elle  stocke le vecteur SDVector et les nucléotides dans deux fichiers distinct.

#### 8. La méthode 'load_from_file(self, dir_path, column_nb)':
Elle charge un objet pysdsl.SDVector et les nucléotides à partir des fichiers créés par la méthode store_to_file().


## Exécution  du Programme :

Pour exécuter le programme, utilisez la commande suivante dans votre terminal :
```
python script_test.py --file "exemple.fasta" --ncols 1000 --compressed --infos --performance 
```
