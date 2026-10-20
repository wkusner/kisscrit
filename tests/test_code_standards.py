"""Code standards (transfer prompt v1.2, 'readable by a referee'; checklist F9).

Status: untrusted (a test of the code's documentation, not of the mathematics).

Every Python file under kisscrit/ and scripts/ has a module docstring with a line 'Status: trusted' or
'Status: untrusted'.  Files of the trusted base must say trusted and carry an 'Implements:' line naming the
lemmas or spec sections they implement, and every function in them has a docstring (the mathematical claim).
Shell scripts under scripts/full/ start with a comment header.
"""
import ast, glob, os, pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]

TRUSTED = [
    'kisscrit/core/interval.py', 'kisscrit/core/model.py', 'kisscrit/core/checker.py', 'kisscrit/core/certificate.py',
    'scripts/completeness/check_certificates.py', 'scripts/completeness/lemmas.py',
    'kisscrit/completeness/candidates.py', 'kisscrit/completeness/canon.py', 'scripts/r1/exact_lambda.py',
]


def _python_files():
    out = []
    for pat in ('kisscrit/**/*.py', 'scripts/**/*.py'):
        out += [p for p in glob.glob(str(ROOT / pat), recursive=True)]
    return sorted(os.path.relpath(p, ROOT) for p in out)


def test_every_python_file_has_a_status_line():
    bad = []
    for rel in _python_files():
        doc = ast.get_docstring(ast.parse((ROOT / rel).read_text()))
        head = '\n'.join((doc or '').splitlines()[:25])
        if 'Status: trusted' not in head and 'Status: untrusted' not in head:
            bad.append(rel)
    assert not bad, 'files without a Status line in the module docstring: %s' % bad


def test_trusted_files_say_trusted_and_name_what_they_implement():
    for rel in TRUSTED:
        doc = ast.get_docstring(ast.parse((ROOT / rel).read_text())) or ''
        assert 'Status: trusted' in doc, rel
        assert 'Implements:' in doc, rel


def test_trusted_functions_have_docstrings():
    missing = []
    for rel in TRUSTED:
        tree = ast.parse((ROOT / rel).read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.FunctionDef) and not ast.get_docstring(node):
                missing.append('%s:%s' % (rel, node.name))
    assert not missing, missing


def test_untrusted_files_do_not_claim_trust():
    for rel in _python_files():
        if rel in TRUSTED: continue
        doc = ast.get_docstring(ast.parse((ROOT / rel).read_text())) or ''
        assert 'Status: trusted' not in doc, '%s claims to be trusted but is not in the trusted list' % rel


def test_shell_scripts_have_headers():
    for p in glob.glob(str(ROOT / 'scripts' / 'full' / '*.sh')):
        lines = open(p).read().splitlines()[:6]
        assert any(l.startswith('#') and len(l) > 10 for l in lines[1:]), p


TRUSTED_C = ['csrc/kc_arith.c', 'csrc/kc_arith.h', 'csrc/kc_model.c', 'csrc/kc_model.h', 'csrc/kc_replay.c', 'csrc/kc_replay.h', 'csrc/kc_leaf.c', 'csrc/kc_leaf.h']


def test_c_files_have_status_headers():
    """Every C source and header starts with a comment that says Status: trusted or untrusted; the trusted
    ones say trusted and name what they implement."""
    for p in glob.glob(str(ROOT / 'csrc' / '*.[ch]')):
        rel = os.path.relpath(p, ROOT)
        head = '\n'.join(open(p).read().splitlines()[:30])
        assert head.lstrip().startswith('/*'), '%s does not start with a header comment' % rel
        assert 'Status: trusted' in head or 'Status: untrusted' in head, '%s has no Status line' % rel
        if rel in TRUSTED_C:
            assert 'Status: trusted' in head and 'Implements:' in head, rel
        else:
            assert 'Status: trusted' not in head, '%s claims to be trusted but is not in the trusted list' % rel
