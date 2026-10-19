# Geometry of Integers (`foliated.py`)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23456718.svg)](https://doi.org/10.5281/zenodo.23456718)
[![Portal](https://img.shields.io/badge/Interactive_Portal-new--mathematical--concepts.org-blue)](https://www.new-mathematical-concepts.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/Python-3.8%2B-yellow.svg)](https://www.python.org/)

Reference exact-integer Python library for the foundational manuscript:
> **Hassan Baloui (2026).** *Geometry of Integers*. Zenodo. [https://doi.org/10.5281/zenodo.23456718](https://doi.org/10.5281/zenodo.23456718)

Explore the interactive visualizers and widgets at **[www.new-mathematical-concepts.org](https://www.new-mathematical-concepts.org)**.

---

## Overview

Classical Peano arithmetic embeds the integers into a flat, one-dimensional sequence ($m=1$) where an integer's distance from the origin coincides with its scalar magnitude. **Geometry of Integers** generalizes this architecture by organizing the integers into an infinite hierarchy of $m$-degree discrete foliations $\mathbb{N}_m$ and their symmetric double-lattice completions $\mathbb{Z}_m$ ($m \ge 2$), partitioned into horizontal shells $S_{k,m} = [k^m, (k+1)^m - 1]$ aligned along a central median spine.

Decoupling an integer's longitudinal depth along its directional leaves ($\lambda \in \{0, \tfrac{1}{2}, 1\}$) from its ambient scalar value induces a directional discrete transport operator $\widehat{\nabla}_{\lambda,m}$ that simultaneously generates:

1. **A Global $L^\infty$ Chebyshev Geodesic Metric ($\delta_m$):** Exact 1-zigzag geodesic routing across both the *Diagonal-to-Vertical* ($\Delta k \ge \vert{}\Delta g_v\vert{}$) and *Diagonal-to-Shell* ($\vert{}\Delta g_v\vert{} > \Delta k$) regimes.
2. **Non-Commutative Multivalued Foliated Arithmetic ($\oplus_m, \ominus_m, \otimes_m$):** Every arithmetic operation across $\mathbb{Z}_m$ produces a coupled directional triad $\{\lambda = 0, \tfrac{1}{2}, 1\}$ and reduces identically to classical commutative Peano arithmetic in the flat limit $m = 1$ (*Correspondence Principle*).
3. **A Leaf-Local Fundamental Theorem of Arithmetic (Leaf-FTA):** Every integer $a \in \mathbb{N}_m$ ($m \ge 2$) sits at the intersection of three directional leaves (Left, Vertical, Right) and admits a unique prime factorization along each leaf simultaneously, yielding a coupled triad of prime factorizations linked to Goldbach partitions for $m=2$.
4. **Holographic Prime Conservation:** In the quadratic foliation $\mathbb{N}_2$, the count of local leaf-depth primes inside the $n$-th horizontal shell $S_n = [n^2, (n+1)^2 - 1]$ along the diagonal leaves is governed exactly by $\pi(n) + \pi(n-1)$.

---

## Quickstart (Zero Dependencies)

`foliated.py` is a single, self-contained Python module using exact arbitrary-precision integer arithmetic (`math.isqrt` and integer Newton–Raphson root extraction) with **zero external dependencies**.

Download `foliated.py` and run the built-in verification suite directly:

```bash
python3 foliated.py
