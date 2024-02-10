import string, math, itertools


def goulden_jackson(bad_words, alphabet=string.ascii_uppercase):
    gfvs = {w: var(f"G_{w}") for w in bad_words}
    eqns = []
    for end_word in bad_words:
        eq = -x^(len(end_word))
        for i in range(1, len(end_word) + 1):
            sub = end_word[:i]
            for source_word in bad_words:
                if source_word.endswith(sub):
                    eq += -x^(len(end_word) - len(sub))*gfvs[source_word]
        eqns.append(eq == 0)

    soln = solve(eqns, *gfvs.values())
    CB = sum(eq.right() for eq in (soln[0] if len(bad_words) > 1 else soln))
    G = 1 / (1 - len(alphabet)*x - CB)
    return G.numerator() / G.denominator()


def PGF(word, alphabet=["1", "0"]):
    G = goulden_jackson([word], alphabet)
    l = len(word)
    b = len(alphabet)

    P = x*(1 - x)*(1/(1 - x) - G(s=x/b))/(x^l)
    return P.numerator() / P.denominator()


def goulden_jackson2(bad_words, alphabet=string.ascii_uppercase, subs=lambda c: x):
    gfvs = {w: var(f"G_{w}") for w in bad_words}
    eqns = []
    for end_word in bad_words:
        eq = -math.prod(subs(c) for c in end_word)
        for i in range(1, len(end_word) + 1):
            sub = end_word[:i]
            mul = math.prod(subs(c) for c in end_word[i:])
            for source_word in bad_words:
                if source_word.endswith(sub):
                    eq += -mul*gfvs[source_word]
        eqns.append(eq == 0)

    soln = solve(eqns, *gfvs.values())
    CB = sum(eq.right() for eq in (soln[0] if len(bad_words) > 1 else soln))
    G = 1 / (1 - sum(subs(c) for c in alphabet) - CB)
    return G.numerator() / G.denominator()


def penney(k, p=9/10):
    subs=lambda c: {"H": p*x, "T": (1 - p)*x}[c]
    w1, w2 = "H"*k, "T"*k
    G1 = goulden_jackson2([w1, w2], "HT", subs)
    G2 = goulden_jackson2([w2], "HT", subs)

    print(G1)
    print(G2)


    H = G2 - G1
    P = (1 - x + subs("T")*(subs("T")^(k - 1) - subs("T")^k)/(1 - subs("T")^k))*H
    return P.numerator() / P.denominator()

y = var("y")
def penney_pure(k, subs=lambda c: {"H": x, "T": y}[c]):
    w1, w2 = "H"*k, "T"*k
    G1 = goulden_jackson2([w1, w2], "HT", subs)
    G2 = goulden_jackson2([w2], "HT", subs)

    print(G1)
    print(G2)

    H = G2 - G1
    adj = subs("T")^k / (1 + subs("T"))
    # suppose w \in S_n,m, then wH \in S_(n + 1),m BUT wT may not be \in S_n,(m + 1) in the
    # case w ends in (k - 1) consecutive Ts
    # The number of w \in S_n,m s.t. w ends in (k - 1) consecutive T's is given by |S_n,(m - k + 1)|
    # minus the number of S_n,(m - k + 1) ending in T

    # let T_n,m be the no of members of S_n,m ending in T
    # let K_n,m be the no of members of S_n,m ending in 4 Ts

    # Then T_n,m = S_{n, m - 1} - K_{n, m - 1}
    # Then K_n,m = S_{n, m - 4} - T_{n, m - 4}

    # T(x, y) = y * S(x, y) - y * K(x, y)
    # K(x, y) = y^4 * S(x, y) - y^4 * T(x, y)

    # K(x, y) = y^4 * S(x, y) - y^4 * (y*S(x, y) - y*K(x, y))
    # K(x, y)(1 - y^5) = y^4 * S(x, y) - y^5 * S(x, y)
    # K(x, y) = (y^4 - y^5)/(1 - y^5) * S(x, y)

    # Want: S_n,m - S_{n - 1, m} - S_{n, m - 1} + K_{n, m - 1}

    P = (1 - subs("H") - subs("T") + subs("T")*(subs("T")^(k - 1) - subs("T")^k)/(1 - subs("T")^k))*H
    return P.numerator() / P.denominator()


def bf(k=5, n=10):
    arr = []
    good, bad = "H"*k, "T"*k
    for i in range(n + 1):
        sm = 0
        for tup in itertools.product("HT", repeat=i):
            s = "".join(tup)
            if good in s and bad not in s:
                sm += 1
        arr.append(sm)
    return arr


def bf_suffix(k=5, n=10):
    arr = []
    good, bad = "H"*k, "T"*k
    for i in range(n + 1):
        sm = 0
        for tup in itertools.product("HT", repeat=i):
            s = "".join(tup)
            if good in s and bad not in s and good not in s[:-1]:
                sm += 1
        arr.append(sm)
    return arr


def bf_suffix2(k=5, n=10):
    arr = []
    good, bad = "H"*k, "T"*k
    for i in range(n + 1):
        sm = 0
        for tup in itertools.product("HTD", repeat=i):
            s = "".join(tup)
            if good in s and bad not in s and good not in s[:-1]:
                sm += 1
        arr.append(sm)
    return arr


z = var("z")
def penney_pure2(k, subs=lambda c: {"H": x, "T": y, "D": z}[c]):
    w1, w2 = "H"*k, "T"*k
    G1 = goulden_jackson2([w1, w2], "HTD", subs)
    G2 = goulden_jackson2([w2], "HTD", subs)

    print(G1)
    print(G2)

    H = G2 - G1
    adj = subs("T")^k / (1 + subs("T"))

    P = (1 - subs("H") - subs("D") - subs("T") + subs("T")*(subs("T")^(k - 1) - subs("T")^k)/(1 - subs("T")^k))*H
    return P.numerator() / P.denominator()



def ponder(p_win, p_lose, k=13):
    w1, w2 = "H"*k, "T"*k
    p_draw = 1 - p_win - p_lose
    subs = lambda c: {"H": p_win*x, "T": p_lose*x, "D": p_draw*x}[c]
    G1 = goulden_jackson2([w1, w2], "HTD", subs)
    G2 = goulden_jackson2([w2], "HTD", subs)
    H = G2 - G1

    P = (1 - subs("H") - subs("D") - subs("T") + subs("T")*(subs("T")^(k - 1) - subs("T")^k)/(1 - subs("T")^k))*H
    return P


def prob():
    base_sum = [sum(1/n * x^i for i in range(1, n + 1)) for n in [4, 6, 8, 12, 20]]
    die = math.prod(base_sum + [sum(1/10 * x^i for i in range(0, 9 + 1))])

    p_win = p_lose = p_draw = 0
    for p, coeff in die.coefficients():
        coeff = int(coeff)
        if is_prime(coeff): p_win += p
        elif not is_prime(coeff) and coeff % 2 == 0: p_lose += p
        else: p_draw += p

    return p_win, p_lose, p_draw
        

def prob2():
    base_sum = [sum(1/n * x^i for i in range(1, n + 1)) for n in [4, 6, 8, 12, 20]]
    die = math.prod(base_sum + [sum(1/10 * x^i for i in range(0, 9 + 1))])

    p_win = p_lose = p_draw = 0
    for p, coeff in die.coefficients():
        coeff = int(coeff)
        if is_prime(coeff): p_win += p
        elif not is_prime(coeff) and coeff % 2 == 1: p_lose += p
        else: p_draw += p

    return p_win, p_lose, p_draw
        
