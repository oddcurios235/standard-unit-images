# Images of Multilinear Polynomials with a Standard-Polynomial Unit Specialization

**Y. T. Miah** · Preprint dated 8 November 2026

[Read the manuscript (PDF)](main.pdf) · [LaTeX source](main.tex) · [Finite verification](REPRODUCIBILITY.md)

## About

This preprint studies the sets of individual matrix values of ordinary multilinear polynomials whose specialization at a designated unit input is a standard polynomial. For $r\in\{3,4,5\}$, suppose that

$$
f(I,X_1,\ldots,X_{2r})=\lambda s_{2r}(X_1,\ldots,X_{2r}),\qquad \lambda\ne 0
$$

as a formal polynomial. Over an algebraically closed field $K$ of characteristic zero, for every $n\ge r+1$ the full single-evaluation image of $f$ on $M_n(K)$ is $\mathfrak{sl}_n(K)$ when its trace vanishes identically, and $M_n(K)$ otherwise. The theorem applies in degrees 7, 9 and 11, including arbitrary proper perturbations of $s_7$, $s_9$ and $s_{11}$.

The argument builds on the standard-polynomial image classification in Miah's [*The L’vov–Kaplansky Conjecture for Standard Polynomials*](https://hexagonmath.org/2611.00184v1), together with geometric image methods of Tang and Pan. Twelve exact standard-polynomial cores, additive size transfer and a mixed-frame compatibility argument provide the further ingredients.

## Files

- [`main.pdf`](main.pdf), [`main.tex`](main.tex), [`references.bib`](references.bib) and [`main.bbl`](main.bbl): the manuscript and its bibliography.
- [`code/verify_standard_cores.py`](code/verify_standard_cores.py): the finite core verifier.
- [`certificates/standard_cores.json`](certificates/standard_cores.json): twelve integer core inputs, with [`certificates/standard_core_checks.json`](certificates/standard_core_checks.json) containing their recorded modular checks.
- [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md) and [`requirements.txt`](requirements.txt): computational scope, instructions and Python dependency.
- [`CITATION.cff`](CITATION.cff): machine-readable citation metadata.

## Build the manuscript

With a TeX distribution providing the packages used in `main.tex`, pdfLaTeX, latexmk and BibTeX, run from this directory:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

The supplied `main.bbl` also permits a PDF build using just pdfLaTeX (run twice, or until cross-references settle). BibTeX and `references.bib` allow the bibliography to be regenerated when needed.

## Reproduce the finite checks

Use Python 3.11 or later with NumPy 1.26 or later (below 3), and run from the directory containing `main.tex`:

```sh
python -m pip install -r requirements.txt
python code/verify_standard_cores.py
```

The script recomputes the endpoint determinants and full-rank normal-map minors for twelve $(r,N)$ pairs, over each of the primes 101 and 103, and rewrites `certificates/standard_core_checks.json`. These finite certificates are used alongside the manuscript's separate arguments for arbitrary sizes. See [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md) for the formulas and scope of the checks.

## Citation

Y. T. Miah, *Images of Multilinear Polynomials with a Standard-Polynomial Unit Specialization*, preprint dated 8 November 2026. Citation metadata is in [`CITATION.cff`](CITATION.cff).

## Licensing

Original scholarly material and data are licensed under [CC BY 4.0](LICENSES/CC-BY-4.0.txt); original verification code and related scripts are under [BSD 3-Clause](LICENSES/BSD-3-Clause.txt). See [LICENSE.md](LICENSE.md) for the exact file scope.
