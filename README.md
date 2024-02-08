
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
"""js
mamba install biopython=1.76
"""

## 4. Fonctionnalités

### Classe 'SuccinctMultipleAlignement`:


#### 1. conctruction de l'alignement multiple succinct: 
Cette classe permet de construire un alignement multiple à partir d'un fichier FASTA contenant les Séquences alignées. On peut specifier le nombre de colonnes à considérer et le type de vecteur (SDVector ou BitVector) à utiliser pour utiliser pour la représentation.


#### 2. Extraction de colonnes specifiques

#### 3. Calcul de la taille en bytes de l'alignement 


#### 4.Extraction de nucléotides specifique dans une colonne 


####  5. Récupération de la sequence d'une colonne 

#### 6. Récupération du vecteur SDVector 

#### 7.Récupération des nucléotides conservés

#### 8.Calcul de la fréquence des nucléotides dans une colonne
#### 9. Sauvegarde de la taille des colonnes dans un fichier CSVqaz


####
