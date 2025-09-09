import random
from Bio.Align import substitution_matrices

# --- Classe principale BLAST simplifié pour protéines ---
class MiniBlastProtein:
    def __init__(self, query, mini_banque, blosum="BLOSUM62", X=5, V=12):
        self.query = query
        self.mini_banque = mini_banque
        self.X = X
        self.V = V
        self.blosum = substitution_matrices.load(blosum)
    
    def find_hits(self, seq, w, T):
        hits = []
        for i in range(len(self.query) - w + 1):
            word = self.query[i:i+w]
            for j in range(len(seq) - w + 1):
                target = seq[j:j+w]
                score = sum(self.blosum[a,b] for a,b in zip(word, target))
                if score >= T:
                    hits.append((i, j, score, word, target))
        return hits
    
    def extend_hit(self, seq, q_start, s_start, length):
        best_score = score = sum(self.blosum[a,b] for a,b in zip(
            self.query[q_start:q_start+length],
            seq[s_start:s_start+length]
        ))
        best_segment = (q_start, q_start+length, s_start, s_start+length)
        
        # Extension gauche
        i, j = q_start-1, s_start-1
        while i >=0 and j >=0:
            score += self.blosum[self.query[i], seq[j]]
            if score > best_score:
                best_score = score
                best_segment = (i, best_segment[1], j, best_segment[3])
            if best_score - score >= self.X:
                break
            i, j = i-1, j-1
        
        # Extension droite
        i, j = best_segment[1], best_segment[3]
        while i < len(self.query) and j < len(seq):
            score += self.blosum[self.query[i], seq[j]]
            if score > best_score:
                best_score = score
                best_segment = (best_segment[0], i+1, best_segment[2], j+1)
            if best_score - score >= self.X:
                break
            i, j = i+1, j+1
        
        if best_score >= self.V:
            return best_segment, best_score
        else:
            return None, best_score
    
    def run(self, w_values, T_values):
        for w in w_values:
            for T in T_values:
                print(f"\n--- Paramètres w={w}, T={T} ---")
                for name, seq in self.mini_banque.items():
                    hits = self.find_hits(seq, w, T)
                    HSPs = []
                    for q_start, s_start, score, word, target in hits:
                        hsp, hsp_score = self.extend_hit(seq, q_start, s_start, w)
                        if hsp:
                            HSPs.append((hsp, hsp_score, word, target))
                    if HSPs:
                        MSP = max(HSPs, key=lambda x: x[1])
                        print(f">{name} | #HSPs={len(HSPs)} | MSP score={MSP[1]}")
                        print(f"Query[{MSP[0][0]}:{MSP[0][1]}] vs {name}[{MSP[0][2]}:{MSP[0][3]}]")
                        print(f"Query: {self.query[MSP[0][0]:MSP[0][1]]}")
                        print(f"Subj : {seq[MSP[0][2]:MSP[0][3]]}")

# --- Génération d'une query aléatoire ---
amino_acids = "ACDEFGHIKLMNPQRSTVWY"
def random_protein(length=60):
    return ''.join(random.choice(amino_acids) for _ in range(length))

query = random_protein(60)

# --- Mini banque de séquences ---
mini_banque = {
    "Prot1": random_protein(100),
    "Prot2": random_protein(120),
    "Prot3": random_protein(80),
    "Prot4": random_protein(110),
    "Prot5": random_protein(90)
}

# --- Paramètres ---
w_values = [3,4,5]   # longueur des mots
T_values = [10,15,20] # score minimal pour considérer un hit
X = 5
V = 12

# --- Exécution ---
blast = MiniBlastProtein(query, mini_banque, X=X, V=V)
blast.run(w_values, T_values)

