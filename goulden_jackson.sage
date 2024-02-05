import string


def goulden_jackson(bad_words, alphabet=string.ascii_uppercase):
    s, gfvs = var("s"), {w: var(f"G_{w}") for w in bad_words}
    eqns = []
    for end_word in bad_words:
        eq = -s^(len(end_word))
        for i in range(1, len(end_word) + 1):
            sub = end_word[:i]
            for source_word in bad_words:
                if source_word.endswith(sub):
                    eq += -s^(len(end_word) - len(sub))*gfvs[source_word]
        eqns.append(eq == 0)

    soln = solve(eqns, *gfvs.values())
    CB = sum(eq.right() for eq in (soln[0] if len(bad_words) > 1 else soln))
    G = 1 / (1 - len(alphabet)*s - CB)
    return G.numerator() / G.denominator()


def PGF(word, alphabet=["1", "0"]):
    G = goulden_jackson([word], alphabet)
    l = len(word)
    b = len(alphabet)

    P = x*(1 - x)*(1/(1 - x) - G(s=x/b))/(x^l)
    return P.numerator() / P.denominator()
