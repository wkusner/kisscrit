"""V6: the checker (kisscrit/core/*) imports search code only as an untrusted oracle: the only module it may
import from kisscrit.completeness is `fast`, and only inside Checker.test_linear_bound."""
import ast, pathlib

CORE = pathlib.Path(__file__).resolve().parents[1] / 'kisscrit' / 'core'


def _imports(tree):
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for a in node.names: yield node, a.name
        elif isinstance(node, ast.ImportFrom):
            mod = node.module or ''
            for a in node.names: yield node, (mod + '.' + a.name) if mod else a.name


def test_core_imports_only_fast_from_completeness():
    for path in CORE.glob('*.py'):
        tree = ast.parse(path.read_text())
        for node, name in _imports(tree):
            if 'completeness' in name:
                assert name.endswith('completeness.fast') or name == 'fast', '%s imports %s' % (path.name, name)
                assert path.name == 'checker.py'
                # the import must sit inside test_linear_bound
                enclosing = [f.name for f in ast.walk(tree) if isinstance(f, ast.FunctionDef)
                             and f.lineno <= node.lineno <= max(getattr(f, 'end_lineno', f.lineno), f.lineno)]
                assert 'test_linear_bound' in enclosing, '%s: completeness import outside test_linear_bound' % path.name
        src = path.read_text()
        for i, line in enumerate(src.splitlines(), 1):
            if 'fast.' in line and path.name == 'checker.py':
                enclosing = [f.name for f in ast.walk(tree) if isinstance(f, ast.FunctionDef)
                             and f.lineno <= i <= max(getattr(f, 'end_lineno', f.lineno), f.lineno)]
                assert 'test_linear_bound' in enclosing, 'checker.py line %d uses fast outside test_linear_bound' % i
