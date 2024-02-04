"""
https://oeis.org/A000073

a_0 == -(x^5 + x^4 + x^3)/(x^3 + x^2 + x - 1)
a_1 == -(x^5 + x^4 + x^3)/(x^3 + x^2 + x - 1)
a_2 == -(x^5 + x^4 + x^3)/(x^3 + x^2 + x - 1)
a_3 == -(x^5 + x^4 + x^3)/(x^3 + x^2 + x - 1)
a_4 == -(x^4 + x^3)/(x^3 + x^2 + x - 1)
a_5 == -(x^4 + x^3)/(x^3 + x^2 + x - 1)
a_6 == -x^3/(x^3 + x^2 + x - 1)
a_7 == 0
"""

def solve_example():
    a_000, a_001, a_010, a_011, a_100, a_101, a_110, a_111 = [var(f"a_{i}") for i in range(8)]

    sys = [a_000 == x^3 + x*(a_100 + a_000),
           a_001 == x^3 + x*(a_100 + a_000),
           a_010 == x^3 + x*(a_101 + a_001),
           a_011 == x^3 + x*(a_101 + a_001),
           a_100 == x^3 + x*(a_110 + a_010),
           a_101 == x^3 + x*(a_110 + a_010),
           a_110 == x^3 + x*(a_111 + a_011),
           a_111 == 0]

    return (solve(sys, a_000, a_001, a_010, a_011, a_100, a_101, a_110, a_111))
