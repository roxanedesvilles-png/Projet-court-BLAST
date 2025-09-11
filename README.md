# Mini BLAST Demo

## Description

Ce projet propose une implémentation simplifiée de l’algorithme BLAST (Basic Local Alignment Search Tool), codée en Python.
Le programme reprend les grandes étapes décrites dans la publication originale (Altschul & Lipman, 1990) :

- génération d’une séquence query aléatoire,

- création d’une banque de séquences aléatoires,

- recherche de hits initiaux (k-mers),

- extension en HSP (High-Scoring Segment Pairs),

- identification du MSP (Maximal Scoring Pair).

⚠️ Cette version est volontairement simplifiée : elle ne calcule pas la probabilité p, et les scores sont basés sur un système match/mismatch basique.

### Contenu

Sur la branche principale : 

- mini_blast_demo_main_script.py → script principal (implémentation de l'algorithme)

- README.md → ce fichier   

- Mini-BLAST.ipynb → ébauche du script principal (notebook)

- article_blast.pdf → la publication originale


Sur la branche additionnelle :

- results_1_exec.txt → résultats obtenus après une execution du main script

- mini_blast_prot_test.py → ébauche de script d'implémentation de BLAST pour les séquences protéiques


## Utilisation 

### Exécution de base

python3 ./mini_blast_demo_main_script.py

### Paramètres (facultatifs)

Dans la fonction mini_blast_demo, les arguments permettent d’ajuster :

- query_len : longueur de la query (défaut : 120)

- seq_count : nombre de séquences dans la banque (défaut : 5)

- seq_min, seq_max : longueurs min/max des séquences de la banque

- k_range : plage de longueurs de k-mers testées (ex: (6,12))

- V_range : plage de valeurs de seuil (ex: (20,45))

- X : paramètre d’arrêt de l’extension

- match, mismatch : scores pour les alignements


### Exemple de sortie

Query length = 120, match=5, mismatch=-4, X=12
Query : ATG...

--- Test k=6, V=20 ---
Seq1 : #HSPs=3, MSP score=35 | k-mer=ATGCCA | Query[15:28] vs Seq[67:80]
  HSPs détaillés :
    ATGCCA | score=30 | Q[15:25] vs Seq[67:77]
    ...


### Prérequis

Python ≥ 3.8

Le script ne nécessite aucune dépendance externe.

Optionnel : pour travailler dans un environnement isolé, on peut utiliser uv : 

uv init

uv run python3 mini_blast_demo_main_script.py


## Références

Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ. Basic local alignment search tool. J Mol Biol. 1990 Oct 5;215(3):403-10. doi: 10.1016/S0022-2836(05)80360-2. PMID: 2231712.



