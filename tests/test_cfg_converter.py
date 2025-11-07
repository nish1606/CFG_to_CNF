import pytest
from cfg_converter import parse_grammar_from_lines, CFGtoChomsky

SAMPLE = [
    "S → ASB",
    "A → aAS | a | ε",
    "B → SbS | A | bb",
]

def test_parse_grammar():
    cfg = parse_grammar_from_lines(SAMPLE)
    assert 'S' in cfg.productions
    assert 'A' in cfg.productions
    assert 'B' in cfg.productions
    assert any(p == ['A','S','B'] for p in cfg.productions['S'])

def is_cnf(cfg):
    # Returns True if grammar is in CNF (A->a or A->BC)
    for lhs, prods in cfg.productions.items():
        for rhs in prods:
            if len(rhs) == 1:
                sym = rhs[0]
                if sym == 'ε':
                    return False
                # allow single terminal
                if not (sym.islower() or (len(sym)==1 and not sym.isupper())):
                    return False
            elif len(rhs) == 2:
                if not (rhs[0].isupper() and rhs[1].isupper()):
                    return False
            else:
                return False
    return True


def test_conversion_to_cnf():
    cfg = parse_grammar_from_lines(SAMPLE)
    conv = CFGtoChomsky(cfg)
    cnf = conv.convert(verbose=False)
    assert is_cnf(cnf), "Converted grammar is not CNF"
