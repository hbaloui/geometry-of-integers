"""
foliated.py — Reference Implementation for "Geometry of Integers"
Author: Dr. Hassan Baloui
DOI: https://doi.org/10.5281/zenodo.23456718
Portal: https://www.new-mathematical-concepts.org

A zero-dependency, exact-integer Python library implementing m-degree integer
foliations (N_m and Z_m), directional discrete transport operators, 1-zigzag
Chebyshev geodesics, non-commutative foliated arithmetic, leaf-local prime
factorizations (Leaf-FTA), and holographic prime conservation.
"""

from __future__ import annotations
import math
from dataclasses import dataclass
from typing import Dict, List, Tuple, Union


# ============================================================================
# 1. EXACT ARBITRARY-PRECISION INTEGER UTILITIES
# ============================================================================

def floor_root(n: int, m: int) -> int:
    """
    Return the exact integer floor of the m-th root of n: floor(n**(1/m)).
    Uses exact integer arithmetic (avoiding floating-point inaccuracy).
    """
    if n < 0:
        raise ValueError("floor_root requires n >= 0.")
    if m < 1:
        raise ValueError("Foliation degree m must be >= 1.")
    if n <= 1 or m == 1:
        return n
    if m == 2:
        return math.isqrt(n)

    # Initial bit-length estimate followed by exact integer Newton-Raphson
    x = 1 << ((n.bit_length() - 1) // m + 1)
    while True:
        x_next = ((m - 1) * x + n // (x ** (m - 1))) // m
        if x_next >= x:
            break
        x = x_next
    while (x + 1) ** m <= n:
        x += 1
    while x ** m > n:
        x -= 1
    return x


def ceil_half(x: int) -> int:
    """Return exact ceil(x / 2) for integer x."""
    return (x + 1) // 2


def is_prime(n: int) -> bool:
    """Deterministic primality test for positive integer n."""
    if n < 2:
        return False
    if n in (2, 3):
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    limit = math.isqrt(n)
    i = 5
    while i <= limit:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True


def prime_factorization(n: int) -> Dict[int, int]:
    """
    Return the canonical prime factorization of n >= 1 as {prime: exponent}.
    For n == 1, returns {1: 1} (the unit element).
    """
    if n < 1:
        raise ValueError("Factorization requires a positive integer >= 1.")
    if n == 1:
        return {1: 1}
    factors: Dict[int, int] = {}
    d = 2
    temp = n
    while d * d <= temp:
        while temp % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1 if d == 2 else 2
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors


def format_factors(factors: Dict[int, int]) -> str:
    """Format a prime factorization dict into a readable string like '2^2 * 3'."""
    if factors == {1: 1}:
        return "1 (unit)"
    parts = [f"{p}^{e}" if e > 1 else f"{p}" for p, e in sorted(factors.items())]
    return " * ".join(parts)


# ============================================================================
# 2. DATA STRUCTURES FOR TRIADS, COORDINATES & GEODESICS
# ============================================================================

@dataclass(frozen=True)
class Triad:
    """Represents a directional triad across gauges lambda in {0, 1/2, 1}."""
    left: int       # lambda = 0
    vertical: int   # lambda = 1/2
    right: int      # lambda = 1

    def as_tuple(self) -> Tuple[int, int, int]:
        return (self.left, self.vertical, self.right)

    def as_dict(self) -> Dict[float, int]:
        return {0.0: self.left, 0.5: self.vertical, 1.0: self.right}

    def __repr__(self) -> str:
        if self.left == self.vertical == self.right:
            return f"{{{self.vertical}}} (m=1 singleton)"
        return f"{{{self.left}, {self.vertical}, {self.right}}}"


@dataclass(frozen=True)
class GeodesicPath:
    """Stores the complete 1-zigzag geodesic result between two integers."""
    source: int
    target: int
    m: int
    distance: int
    delta_k: int
    delta_g: int
    lambda_star: float
    regime: str
    turning_point: int
    turning_shell: int
    virtual_point: int
    virtual_shell: int
    diag_steps: int
    second_leg_steps: int
    description: str


# ============================================================================
# 3. CORE FOLIATION ENGINE (N_m and Z_m)
# ============================================================================

class Foliation:
    """
    Represents the m-degree discrete integer foliation Z_m (and positive half N_m).
    When m = 1, reduces identically to the classical Peano number line.
    """

    def __init__(self, m: int = 2):
        if m < 1:
            raise ValueError("Foliation degree m must be >= 1.")
        self.m = m

    def __call__(self, z: int) -> "FoliatedInt":
        return FoliatedInt(z, self)

    # --- Shell & Coordinate Geometry ---

    def shell_index(self, n: int) -> int:
        """Return shell index k = floor(|n|**(1/m))."""
        return floor_root(abs(n), self.m)

    def shell_width(self, k: int) -> int:
        """Return the number of integers in shell S_{k,m}: (k+1)^m - k^m."""
        if k < 0:
            raise ValueError("Shell index k must be >= 0.")
        return (k + 1) ** self.m - k ** self.m

    def half_width(self, k: int) -> int:
        """
        Return half-width H_{k,m} = ((k+1)^m - k^m - 1) // 2.
        For m=2, H_{k,2} = k. For m=1, H_{k,1} = 0.
        """
        return (self.shell_width(k) - 1) // 2

    def shell_bounds(self, k: int) -> Tuple[int, int]:
        """Return (min_element, max_element) of positive shell S_{k,m}."""
        return (k ** self.m, (k + 1) ** self.m - 1)

    def coordinates(self, n: int) -> Dict[str, int]:
        """
        Return full geometric coordinates for n >= 1 in N_m:
        k (shell), r (residue), g_l (left gen), g_v (vertical gen), g_r (right gen).
        """
        if n < 0:
            raise ValueError("Use coordinates() on magnitude |z| >= 0.")
        if n == 0:
            return {"n": 0, "k": 0, "r": 0, "g_l": 0, "g_v": 0, "g_r": 0}
        k = self.shell_index(n)
        r = n - k ** self.m
        H_k = self.half_width(k)
        g_l = r
        g_v = r - H_k
        g_r = 2 * H_k - r
        return {"n": n, "k": k, "r": r, "g_l": g_l, "g_v": g_v, "g_r": g_r}

    def birth_shells(self, n: int) -> Triad:
        """
        Return the birth shell index k_0(n) along each directional leaf (0, 1/2, 1).
        In m=2: k_0(L) = ceil(r/2), k_0(V) = |g_v|, k_0(R) = ceil(g_r/2).
        """
        if n <= 0:
            return Triad(0, 0, 0)
        if self.m == 1:
            return Triad(1, 1, 1)
        c = self.coordinates(n)
        if self.m == 2:
            k0_l = max(1, ceil_half(c["g_l"]))
            k0_v = max(1, abs(c["g_v"]))
            k0_r = max(1, ceil_half(c["g_r"]))
            return Triad(k0_l, k0_v, k0_r)

        # General m >= 3: smallest k0 >= 1 where shell half-width accommodates the generator
        k0_l = 1
        while 2 * self.half_width(k0_l) < c["g_l"]:
            k0_l += 1
        k0_v = 1
        while self.half_width(k0_v) < abs(c["g_v"]):
            k0_v += 1
        k0_r = 1
        while 2 * self.half_width(k0_r) < c["g_r"]:
            k0_r += 1
        return Triad(k0_l, k0_v, k0_r)

    def leaf_depths(self, n: int) -> Triad:
        """Return 0-indexed local leaf depths d_{lambda,m}(n) = k - k_0."""
        if n <= 0:
            return Triad(0, 0, 0)
        k = self.shell_index(n)
        k0 = self.birth_shells(n)
        return Triad(k - k0.left, k - k0.vertical, k - k0.right)

    def leaf_positions(self, n: int) -> Triad:
        """Return 1-indexed leaf positions N_{lambda,m}(n) = d_{lambda,m}(n) + 1."""
        d = self.leaf_depths(n)
        return Triad(d.left + 1, d.vertical + 1, d.right + 1)

    # --- Discrete Transport Operators ---

    def nabla(self, n: int, lam: float, s: int) -> int:
        """
        Positive half-lattice transport operator nabla_{lambda, m}(n, s).
        Shifts n in S_k to shell K = k + s along directional gauge lam in {0, 0.5, 1}.
        """
        if n < 0:
            raise ValueError("nabla requires n >= 0; use widehat_nabla for Z_m.")
        if self.m == 1:
            return n + s
        k = self.shell_index(n)
        r = n - k ** self.m
        K = k + s
        if K < 0:
            raise ValueError(f"Target shell K = {K} is negative; use widehat_nabla.")
        delta_H = self.half_width(K) - self.half_width(k)
        # 2 * lam is an exact integer in {0, 1, 2}
        two_lam = round(2 * lam)
        return (K ** self.m) + r + two_lam * delta_H

    def widehat_nabla(self, z: int, lam: float, s: int) -> int:
        """
        Extended discrete shift operator on the symmetric double-lattice Z_m.
        Resolves inward shifts that exceed local leaf depth by reflecting across 0.
        """
        if self.m == 1:
            return z + s
        if s == 0:
            return z
        if z == 0:
            sgn_s = 1 if s > 0 else -1
            return sgn_s * self.nabla(1, lam, abs(s) - 1)

        eps_z = 1 if z > 0 else -1
        n = abs(z)
        s_tilde = eps_z * s

        depths = self.leaf_depths(n).as_dict()
        d_lam = depths[float(lam)]

        if s_tilde >= -d_lam:
            return eps_z * self.nabla(n, lam, s_tilde)
        else:
            return -eps_z * self.nabla(n, lam, abs(s_tilde))

    # --- Non-Commutative Foliated Arithmetic on Z_m ---

    def add(self, a: int, b: int) -> Triad:
        """Foliated addition: a (+) b = { widehat_nabla(a, lam, b) }."""
        return Triad(
            self.widehat_nabla(a, 0.0, b),
            self.widehat_nabla(a, 0.5, b),
            self.widehat_nabla(a, 1.0, b),
        )

    def sub(self, a: int, b: int) -> Triad:
        """Foliated subtraction: a (-) b = { widehat_nabla(a, lam, -b) }."""
        return self.add(a, -b)

    def mul(self, a: int, b: int) -> Triad:
        """Foliated multiplication: a (x) b = { widehat_nabla(a, lam, a*(b-1)) }."""
        if a == 0 or b == 0:
            return Triad(0, 0, 0)
        shift = a * (b - 1)
        return Triad(
            self.widehat_nabla(a, 0.0, shift),
            self.widehat_nabla(a, 0.5, shift),
            self.widehat_nabla(a, 1.0, shift),
        )

    # --- 1-Zigzag Chebyshev Geodesic Router ---

    def geodesic(self, x: int, y: int) -> GeodesicPath:
        """
        Compute the exact 1-zigzag Chebyshev geodesic between x and y in N_m,
        distinguishing the Diagonal-to-Vertical and Diagonal-to-Shell regimes.
        """
        if x < 1 or y < 1:
            raise ValueError("Geodesic routing in N_m requires positive integers >= 1.")
        u, v = (x, y) if x <= y else (y, x)
        cu = self.coordinates(u)
        cv = self.coordinates(v)

        delta_k = cv["k"] - cu["k"]
        delta_g = cv["g_v"] - cu["g_v"]
        abs_dg = abs(delta_g)

        if delta_g > 0:
            lam_star = 1.0
        elif delta_g < 0:
            lam_star = 0.0
        else:
            lam_star = 0.5

        # Virtual intersection from Equation (2.11)
        virtual_pt = self.nabla(u, lam_star, abs_dg)
        virtual_k = self.shell_index(virtual_pt)

        # Closed-form metric Equation (2.11)
        distance = (delta_k + abs_dg + abs(cv["k"] - virtual_k)) // 2

        # Physical non-overshooting turning point X_actual = nabla(u, min(delta_k, |delta_g|))
        diag_steps = min(delta_k, abs_dg)
        turning_pt = self.nabla(u, lam_star, diag_steps)
        turning_k = self.shell_index(turning_pt)

        if delta_k >= abs_dg:
            regime = "Diagonal-to-Vertical (delta_k >= |delta_g|)"
            second_steps = delta_k - diag_steps
            desc = (
                f"Diagonal (lambda={lam_star:g}) for {diag_steps} steps from {u} to "
                f"X = {turning_pt} (Shell S_{turning_k}, hitting V_{cv['g_v']}), then "
                f"{second_steps} vertical steps (lambda=0.5) down to {v}."
            )
        else:
            regime = "Diagonal-to-Shell (|delta_g| > delta_k)"
            second_steps = abs_dg - diag_steps
            desc = (
                f"Diagonal (lambda={lam_star:g}) for {diag_steps} steps from {u} to "
                f"X = {turning_pt} (hitting target Shell S_{cv['k']}), then "
                f"{second_steps} horizontal intra-shell steps along S_{cv['k']} to {v}."
            )

        return GeodesicPath(
            source=u,
            target=v,
            m=self.m,
            distance=distance,
            delta_k=delta_k,
            delta_g=delta_g,
            lambda_star=lam_star,
            regime=regime,
            turning_point=turning_pt,
            turning_shell=turning_k,
            virtual_point=virtual_pt,
            virtual_shell=virtual_k,
            diag_steps=diag_steps,
            second_leg_steps=second_steps,
            description=desc,
        )

    # --- Leaf-Local FTA & Holographic Prime Conservation ---

    def factorize_triad(self, n: int) -> Dict[str, str]:
        """
        Return the Leaf-Local Fundamental Theorem of Arithmetic triad of prime
        factorizations for n >= 1 across its Left, Vertical, and Right leaves.
        """
        pos = self.leaf_positions(n)
        return {
            "Left (lambda=0)": f"N_0 = {pos.left} -> {format_factors(prime_factorization(pos.left))}",
            "Vertical (lambda=1/2)": f"N_1/2 = {pos.vertical} -> {format_factors(prime_factorization(pos.vertical))}",
            "Right (lambda=1)": f"N_1 = {pos.right} -> {format_factors(prime_factorization(pos.right))}",
        }

    def holographic_shell_primes(self, k: int) -> Dict[str, int]:
        """
        Verify the Holographic Prime Conservation Law in shell S_{k,2} (for m=2):
        Counts leaf-depth primes along Left/Right/Vertical leaves in shell k and
        compares with classical prime counts pi(k) + pi(k-1).
        """
        if self.m != 2:
            raise ValueError("Exact pi(k) + pi(k-1) identity applies to m = 2.")
        lo, hi = self.shell_bounds(k)
        left_primes = sum(1 for a in range(lo, hi + 1) if is_prime(self.leaf_positions(a).left))
        right_primes = sum(1 for a in range(lo, hi + 1) if is_prime(self.leaf_positions(a).right))
        pi_k = sum(1 for j in range(1, k + 1) if is_prime(j))
        pi_k_minus_1 = sum(1 for j in range(1, k) if is_prime(j))
        return {
            "shell_k": k,
            "left_leaf_primes_in_S_k": left_primes,
            "right_leaf_primes_in_S_k": right_primes,
            "pi(k) + pi(k-1)": pi_k + pi_k_minus_1,
        }


class FoliatedInt:
    """Object-oriented wrapper around an integer z in Z_m."""

    def __init__(self, z: int, foliation: Foliation):
        self.z = int(z)
        self.foliation = foliation

    @property
    def coords(self) -> Dict[str, int]:
        return self.foliation.coordinates(abs(self.z))

    @property
    def leaf_positions(self) -> Triad:
        return self.foliation.leaf_positions(abs(self.z))

    def factorize(self) -> Dict[str, str]:
        return self.foliation.factorize_triad(abs(self.z))

    def __add__(self, other: Union["FoliatedInt", int]) -> Triad:
        b = other.z if isinstance(other, FoliatedInt) else int(other)
        return self.foliation.add(self.z, b)

    def __sub__(self, other: Union["FoliatedInt", int]) -> Triad:
        b = other.z if isinstance(other, FoliatedInt) else int(other)
        return self.foliation.sub(self.z, b)

    def __mul__(self, other: Union["FoliatedInt", int]) -> Triad:
        b = other.z if isinstance(other, FoliatedInt) else int(other)
        return self.foliation.mul(self.z, b)

    def __repr__(self) -> str:
        return f"FoliatedInt({self.z}, m={self.foliation.m})"


# ============================================================================
# 4. BUILT-IN SELF-VERIFICATION SUITE (Runs when executed directly)
# ============================================================================

if __name__ == "__main__":
    print("=" * 72)
    print(" GEOMETRY OF INTEGERS — AUTOMATED VERIFICATION SUITE")
    print("=" * 72)

    Z2 = Foliation(m=2)
    Z1 = Foliation(m=1)

    # 1. Verify Non-Commutative Multiplication on Z_2
    print("\n[1] Non-Commutative Multiplication on Z_2:")
    print(f"   -2 (x)_2  4  = {Z2(-2) * Z2(4)}")
    print(f"    2 (x)_2 -4  = {Z2(2) * Z2(-4)}")
    print(f"    2 (x)_2  4  = {Z2(2) * Z2(4)}")
    print(f"    4 (x)_2  2  = {Z2(4) * Z2(2)}")

    # 2. Verify Correspondence Principle (m = 1 Collapse)
    print("\n[2] Correspondence Principle (m = 1 Peano Limit):")
    print(f"   -2 (x)_1  4  = {Z1(-2) * Z1(4)}")
    print(f"    2 (x)_1 -4  = {Z1(2) * Z1(-4)}")

    # 3. Verify 1-Zigzag Geodesic (64 -> 886)
    print("\n[3] 1-Zigzag Geodesic Routing (64 -> 886 in N_2):")
    geo = Z2.geodesic(64, 886)
    print(f"   Regime          : {geo.regime}")
    print(f"   Metric Distance : delta_2(64, 886) = {geo.distance}")
    print(f"   Physical Turn X : {geo.turning_point} (in Shell S_{geo.turning_shell})")
    print(f"   Virtual Turn    : {geo.virtual_point} (in Shell S_{geo.virtual_shell})")
    print(f"   Path Summary    : {geo.description}")

    # 4. Verify Leaf-Local Fundamental Theorem of Arithmetic (n = 21)
    print("\n[4] Leaf-Local FTA Triple Prime Factorization for a = 21 in N_2:")
    for leaf, fac in Z2(21).factorize().items():
        print(f"   {leaf:22s}: {fac}")

    # 5. Verify Holographic Prime Conservation Law in Shell S_10
    print("\n[5] Holographic Prime Conservation Law (Shell S_10 = [100, 120]):")
    holo = Z2.holographic_shell_primes(10)
    print(f"   Left-leaf primes in S_10  : {holo['left_leaf_primes_in_S_k']}")
    print(f"   Right-leaf primes in S_10 : {holo['right_leaf_primes_in_S_k']}")
    print(f"   Classical pi(10) + pi(9)  : {holo['pi(k) + pi(k-1)']}")
    print("=" * 72)
