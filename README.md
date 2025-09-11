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

## Contenu

Sur la branche principale : 

- mini_blast_demo_main_script.py → script principal (implémentation de l'algorithme)

- README.md → ce fichier 

- DESVILLES_Roxane_rapport.pdf → le rapport du projet en pdf 

- 
