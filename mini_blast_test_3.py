import random

# ------------------------
# Génération aléatoire ADN
# ------------------------
def random_dna(length):
    return ''.join(random.choice("ATGC") for _ in range(length))

# ------------------------
# Extension sans gap
# ------------------------
def extend_hit(query, seq, q_start, s_start, length, X=3, V=5, match=1, mismatch=-1):
    best_score = score = sum(match if query[q_start+i] == seq[s_start+i] else mismatch for i in range(length))
    best_segment = (q_start, q_start+length, s_start, s_start+length)

    # --- Extension à gauche ---
    i, j = q_start-1, s_start-1
    while i >=0 and j >=0:
        score += match if query[i]==seq[j] else mismatch
        if score > best_score:
            best_score = score
            best_segment = (i, best_segment[1], j, best_segment[3])
        if best_score - score >= X:
            break
        i, j = i-1, j-1

    # --- Extension à droite ---
    i, j = best_segment[1], best_segment[3]
    while i < len(query) and j < len(seq):
        score += match if query[i]==seq[j] else mismatch
        if score > best_score:
            best_score = score
            best_segment = (best_segment[0], i+1, best_segment[2], j+1)
        if best_score - score >= X:
            break
        i, j = i+1, j+1

    if best_score >= V:
        return best_segment, best_score
    else:
        return None, best_score

# ------------------------
# Mini-BLAST paramétrique complet
# ------------------------
def mini_blast_demo(query_len=120, seq_count=5, seq_min=140, seq_max=200,
                    k_range=(6,12), V_range=(20,45), X=12, match=5, mismatch=-4):
    # Génération query
    query = random_dna(query_len)
    # Génération banque
    bank = {f"Seq{i+1}": random_dna(random.randint(seq_min, seq_max)) for i in range(seq_count)}

    print(f"Query length = {query_len}, match={match}, mismatch={mismatch}, X={X}")
    print(f"Query : {query}\n")

    ks = list(range(k_range[0], k_range[1]))         # Ex: 6 à 11
    Vs = list(range(V_range[0], V_range[1], 5))     # Ex: 20,25,30,35,40

    for k in ks:
        for V in Vs:
            print(f"--- Test k={k}, V={V} ---")
            hits = {}
            # Recherche hits initiaux
            for name, seq in bank.items():
                hits_seq = []
                for i in range(len(query)-k+1):
                    km = query[i:i+k]
                    pos = seq.find(km)
                    if pos != -1:
                        hits_seq.append((km, pos))
                if hits_seq:
                    hits[name] = hits_seq

            # Extension
            HSPs = {}
            for name, seq_hits in hits.items():
                seq = bank[name]
                hsp_list = []
                for km, pos in seq_hits:
                    q_pos = query.find(km)
                    hsp, score = extend_hit(query, seq, q_pos, pos, len(km), X=X, V=V, match=match, mismatch=mismatch)
                    if hsp:
                        hsp_list.append((km, hsp, score))
                if hsp_list:
                    # MSP = meilleur HSP
                    MSP = max(hsp_list, key=lambda x:x[2])
                    HSPs[name] = {"HSPs": hsp_list, "MSP": MSP}

            # Affichage résumé détaillé
            if HSPs:
                for name, data in HSPs.items():
                    print(f"{name} : #HSPs={len(data['HSPs'])}, MSP score={data['MSP'][2]} | k-mer={data['MSP'][0]} | Query[{data['MSP'][1][0]}:{data['MSP'][1][1]}] vs Seq[{data['MSP'][1][2]}:{data['MSP'][1][3]}]")
                    print("  HSPs détaillés :")
                    for km, seg, score in data['HSPs']:
                        print(f"    {km} | score={score} | Q[{seg[0]}:{seg[1]}] vs Seq[{seg[2]}:{seg[3]}]")
            else:
                print("Aucun HSP trouvé.")
            print()

# ------------------------
# Lancer le script
# ------------------------
if __name__ == "__main__":
    random.seed()  # reseed pour résultat différent à chaque run
    mini_blast_demo()

