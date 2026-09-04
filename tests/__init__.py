"""Test suite for flowfreq.

This file is required: test modules import fixtures by absolute path
(``from tests.fixtures... import ...``), which only resolves when
``tests`` is a package and the repository root is therefore what pytest adds
to ``sys.path``. Without it, ``python -m pytest`` still works -- it puts the
CWD on the path -- but the ``pytest`` console script used in CI does not, and
collection fails with ``ModuleNotFoundError: No module named 'tests'``.
"""
