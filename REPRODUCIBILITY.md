# Reproducing the twelve finite standard-core checks

The finite input file `certificates/standard_cores.json` contains twelve integer tuples 
`(r, N, B_1, ..., B_{2r-1}, u)`. For each tuple the script constructs

- the endpoint `S=s_{2r-1}(B)`;
- the linear normal operator `Γ_{B,u}(P) = (s_{2r}(B_1,...,B_{2r-1},P), C_{1,B}(P)u, ..., C_{2r-1,B}(P)u)`, where `C_{i,B}(P)` is `s_{2r-1}(B_1,...,B_{i-1},P,B_{i+1},...,B_{2r-1})`.

For each modulus 101 and 103 the program computes `S_B`, constructs the linear operator `Γ_{B,u}` using signed insertion operators, chooses `N^2` independent rows of its matrix by modular elimination, and computes the determinant of that selected square minor. It checks that both `det(S_B)` and this minor determinant are nonzero, and records the pivot-row indices, determinant residues, rank, and `seed_col = S_B u` modulo the prime. It does not separately compute the full bordered standard-polynomial Jacobian.

## Requirements and command

Python 3.11 or later and NumPy 1.26 or later (below 3). From the directory containing `main.tex`, run:

```bash
python -m pip install -r requirements.txt
python code/verify_standard_cores.py
```

The expected terminal summary is `PASS: 12 core pairs; 24 full-rank normal maps and nonzero minors (p=101,103).` The program validates each prescribed `(r,N)` pair exactly once, then overwrites `certificates/standard_core_checks.json` with the computed result file. The computation is entirely deterministic and uses no network connection.

## Mathematical interpretation

Full rank modulo either prime proves that some corresponding integer determinant is nonzero, and hence proves the finite admissibility of the twelve listed core sizes over characteristic-zero fields. It does not by itself prove the result at arbitrary sizes: the manuscript's additive-transfer lemma supplies those infinitely many cases. The mixed-frame detection lemma is symbolic and uses the classical sharp degree-`2N` matrix-identity theorem; it is not a computational search.

The computed outcomes depend on the supplied integer input tuples rather than the stored ranks. The generated result file is reproducible from this package alone.
