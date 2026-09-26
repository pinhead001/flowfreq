"""Parse NSS regression-equation strings into :class:`RegressionEquation` form.

NSS (``/nssservices/Scenarios/Estimate``) returns each estimate with the literal
equation it evaluated, for example::

    3.846*DRNAREA^0.745*10^(0.032*PRECPRIS10)/10^(0.0078*CANOPY_PCT)
    27.9*DRNAREA^0.523*(MINBELEV/1000)^(-0.098)
    0.00953*(DRNAREA)^(0.392)*(0.001*(RELIEF))^(3.36)
    5.0*CONTDA^0.674*(EL5000+1)^(-0.026)

This module turns such a string into the log-linear form of
:mod:`flowfreq.regression.equations`::

    log10(Q) = b0 + sum_i b_i * T_i(X_i),   T_i in {identity, log10, log10_plus1}

A term maps as follows:

- ``X^b`` or ``(X)^(b)`` becomes ``log10`` with coefficient ``b``.
- ``(X+1)^b`` becomes ``log10_plus1``.
- ``10^(c*X)`` becomes ``identity`` with ``c``, and ``/10^(c*X)`` gives ``-c``.
- Constant factors, ``10^(k)`` included, fold into ``b0``.
- A scale inside a power term, ``(s*X)^b`` or ``(X/s)^b``, is algebraically
  ``s^b * X^b``. It folds into the intercept as ``b*log10(s)``, so the equation
  stays exact.

``max(0, E)`` is accepted when ``E`` is itself a product of such terms, because
that product is strictly positive wherever its terms are defined.

**Anything else is refused, never approximated.** That covers an offset other
than ``+1`` (``(X+0.01)^b``, ``(X-20)^b``, ``(X/100+1)^b``), a subtracted
constant (``max(0, E-1)``), logistic ``e^(...)`` forms, and a variable that
appears under two different transforms. Such a string raises
:class:`NSSEquationError`. :func:`evaluate_expression` still evaluates it
directly, so a snapshot can record what NSS computes without claiming to
represent it.

Roadmap: ``docs/MASTER_ROADMAP.md`` §3.1, issue #35.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from typing import Dict, List, Mapping, Optional, Sequence, Tuple, Union

from flowfreq.regression.equations import TRANSFORMS, Citation, RegressionEquation, Variable

__all__ = [
    "NSSEquationError",
    "ParsedEquation",
    "evaluate_expression",
    "parse_equation",
    "parse_expression",
    "statistic_code_to_aep",
    "to_regression_equation",
    "variables_in",
]


class NSSEquationError(ValueError):
    """An NSS equation string is malformed, or cannot be represented exactly."""


# ---------------------------------------------------------------------------
# Tokenizer and parser. The AST is plain tuples:
#   ("num", float) | ("var", str) | ("neg", node)
#   ("bin", op, left, right) with op in + - * / ^
#   ("call", name, (args...))
# ---------------------------------------------------------------------------

Node = Tuple
_TOKEN = re.compile(
    r"\s*(?:(?P<num>(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)"
    r"|(?P<name>[A-Za-z_][A-Za-z0-9_]*)"
    r"|(?P<op>[-+*/^(),]))"
)
_VARIABLE = re.compile(r"^[A-Z][A-Z0-9_]*$")
_FUNCTIONS = frozenset({"max"})


def _tokenize(text: str) -> List[Tuple[str, str]]:
    tokens: List[Tuple[str, str]] = []
    pos = 0
    text = text.rstrip()
    while pos < len(text):
        m = _TOKEN.match(text, pos)
        if m is None or m.end() == pos:
            raise NSSEquationError(f"unexpected character {text[pos:pos + 10]!r} in {text!r}")
        kind = m.lastgroup
        assert kind is not None
        tokens.append((kind, m.group(kind)))
        pos = m.end()
    return tokens


class _Parser:
    """Recursive descent. ``^`` binds tightest and is right-associative;
    unary minus binds looser than ``^``, so ``-2^2`` is ``-(2^2)``."""

    def __init__(self, text: str) -> None:
        self.text = text
        self.tokens = _tokenize(text)
        self.i = 0

    def _peek(self) -> Optional[Tuple[str, str]]:
        return self.tokens[self.i] if self.i < len(self.tokens) else None

    def _take(self, value: Optional[str] = None) -> Tuple[str, str]:
        tok = self._peek()
        if tok is None:
            raise NSSEquationError(f"unexpected end of equation {self.text!r}")
        if value is not None and tok[1] != value:
            raise NSSEquationError(f"expected {value!r}, found {tok[1]!r} in {self.text!r}")
        self.i += 1
        return tok

    def parse(self) -> Node:
        if not self.tokens:
            raise NSSEquationError("empty equation")
        node = self._expr()
        extra = self._peek()
        if extra is not None:
            raise NSSEquationError(f"trailing {extra[1]!r} in {self.text!r}")
        return node

    def _expr(self) -> Node:
        node = self._term()
        while (tok := self._peek()) is not None and tok[1] in "+-" and tok[0] == "op":
            self.i += 1
            node = ("bin", tok[1], node, self._term())
        return node

    def _term(self) -> Node:
        node = self._unary()
        while (tok := self._peek()) is not None and tok[1] in "*/" and tok[0] == "op":
            self.i += 1
            node = ("bin", tok[1], node, self._unary())
        return node

    def _unary(self) -> Node:
        tok = self._peek()
        if tok is not None and tok == ("op", "-"):
            self.i += 1
            return ("neg", self._unary())
        if tok is not None and tok == ("op", "+"):
            self.i += 1
            return self._unary()
        return self._power()

    def _power(self) -> Node:
        base = self._primary()
        tok = self._peek()
        if tok is not None and tok == ("op", "^"):
            self.i += 1
            return ("bin", "^", base, self._unary())
        return base

    def _primary(self) -> Node:
        kind, value = self._take()
        if kind == "num":
            return ("num", float(value))
        if kind == "name":
            nxt = self._peek()
            if nxt is not None and nxt == ("op", "("):
                self.i += 1
                args = [self._expr()]
                while self._peek() == ("op", ","):
                    self.i += 1
                    args.append(self._expr())
                self._take(")")
                return ("call", value, tuple(args))
            return ("var", value)
        if value == "(":
            node = self._expr()
            self._take(")")
            return node
        raise NSSEquationError(f"unexpected {value!r} in {self.text!r}")


def parse_expression(text: str) -> Node:
    """Parse an NSS equation string into a tuple AST.

    Raises
    ------
    NSSEquationError
        If the string is not a well-formed arithmetic expression.
    """
    return _Parser(text).parse()


def _variables(node: Node) -> List[str]:
    kind = node[0]
    if kind == "var":
        return [node[1]]
    if kind == "neg":
        return _variables(node[1])
    if kind == "bin":
        return _variables(node[2]) + _variables(node[3])
    if kind == "call":
        return [v for a in node[2] for v in _variables(a)]
    return []


def evaluate_expression(equation: Union[str, Node], values: Mapping[str, float]) -> float:
    """Evaluate an NSS equation string directly, with no log-linear rewriting.

    This is the independent check on :func:`parse_equation`: the two must agree.

    Raises
    ------
    KeyError
        A variable has no value.
    NSSEquationError
        An unknown function, or a power that is not a real number.
    """
    node = parse_expression(equation) if isinstance(equation, str) else equation
    return _eval(node, values)


def _eval(node: Node, values: Mapping[str, float]) -> float:
    kind = node[0]
    if kind == "num":
        return float(node[1])
    if kind == "var":
        if node[1] == "e":
            return math.e
        return float(values[node[1]])
    if kind == "neg":
        return -_eval(node[1], values)
    if kind == "call":
        args = [_eval(a, values) for a in node[2]]
        if node[1] == "max":
            return max(args)
        if node[1] == "exp" and len(args) == 1:
            return math.exp(args[0])
        raise NSSEquationError(f"unknown function {node[1]!r}")
    op, a, b = node[1], _eval(node[2], values), _eval(node[3], values)
    if op == "+":
        return a + b
    if op == "-":
        return a - b
    if op == "*":
        return a * b
    if op == "/":
        return a / b
    result = a**b
    if isinstance(result, complex):
        raise NSSEquationError(f"{a}^{b} is not a real number")
    return float(result)


# ---------------------------------------------------------------------------
# Log-linearization
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ParsedEquation:
    """An NSS equation rewritten exactly as ``log10(Q) = b0 + sum b_i*T_i(X_i)``.

    Attributes
    ----------
    equation : str
        The NSS string it came from.
    intercept : float
        ``b0``, log10 space.
    codes, transforms, coefficients : tuple
        One entry per variable, in order of first appearance in the string.
    notes : tuple of str
        Rewrites a reviewer should know about, such as a scale folded into ``b0``.
    """

    equation: str
    intercept: float
    codes: Tuple[str, ...]
    transforms: Tuple[str, ...]
    coefficients: Tuple[float, ...]
    notes: Tuple[str, ...] = ()

    def log10_value(self, values: Mapping[str, float]) -> float:
        """Evaluate ``log10(Q)`` at the given characteristics."""
        return self.intercept + sum(
            c * TRANSFORMS[t](float(values[k]))
            for k, t, c in zip(self.codes, self.transforms, self.coefficients)
        )

    def to_dict(self) -> Dict[str, object]:
        return {
            "intercept": self.intercept,
            "variables": [{"code": k, "transform": t} for k, t in zip(self.codes, self.transforms)],
            "coefficients": list(self.coefficients),
            "notes": list(self.notes),
        }


class _LogLin:
    """log10 of an expression: ``const + sum coef * T(X)``."""

    def __init__(self) -> None:
        self.const = 0.0
        self.terms: Dict[str, Tuple[str, float]] = {}
        self.order: List[str] = []
        self.notes: List[str] = []

    @classmethod
    def constant(cls, value: float) -> "_LogLin":
        out = cls()
        out.const = value
        return out

    @classmethod
    def term(cls, code: str, transform: str, coef: float = 1.0) -> "_LogLin":
        out = cls()
        out.terms[code] = (transform, coef)
        out.order.append(code)
        return out

    def scaled(self, p: float) -> "_LogLin":
        out = _LogLin()
        out.const = self.const * p
        out.terms = {k: (t, c * p) for k, (t, c) in self.terms.items()}
        out.order = list(self.order)
        out.notes = list(self.notes)
        return out

    def plus(self, other: "_LogLin", text: str) -> "_LogLin":
        out = self.scaled(1.0)
        out.const += other.const
        for k in other.order:
            t, c = other.terms[k]
            if k in out.terms:
                t0, c0 = out.terms[k]
                if t0 != t:
                    raise NSSEquationError(
                        f"{k} appears under two transforms ({t0}, {t}) in {text!r}; "
                        "RegressionEquation allows one term per variable"
                    )
                out.terms[k] = (t, c0 + c)
            else:
                out.terms[k] = (t, c)
                out.order.append(k)
        out.notes += other.notes
        return out


def _const_value(node: Node) -> Optional[float]:
    """The value of a numeric literal, possibly negated; None if not constant."""
    if node[0] == "num":
        return float(node[1])
    if node[0] == "neg":
        inner = _const_value(node[1])
        return None if inner is None else -inner
    return None


def _check_variable(name: str, text: str) -> None:
    if not _VARIABLE.match(name):
        raise NSSEquationError(f"unsupported identifier {name!r} in {text!r}")


def _linear(node: Node, text: str) -> Tuple[float, Dict[str, float], List[str]]:
    """``node`` as ``const + sum c*X`` (the exponent of a ``10^(...)`` term)."""
    kind = node[0]
    if kind == "num":
        return float(node[1]), {}, []
    if kind == "var":
        _check_variable(node[1], text)
        return 0.0, {node[1]: 1.0}, [node[1]]
    if kind == "neg":
        c, d, o = _linear(node[1], text)
        return -c, {k: -v for k, v in d.items()}, o
    if kind == "bin" and node[1] in "+-":
        c1, d1, o1 = _linear(node[2], text)
        c2, d2, o2 = _linear(node[3], text)
        sign = 1.0 if node[1] == "+" else -1.0
        d = dict(d1)
        for k, v in d2.items():
            d[k] = d.get(k, 0.0) + sign * v
        return c1 + sign * c2, d, o1 + [k for k in o2 if k not in o1]
    if kind == "bin" and node[1] in "*/":
        left, right = _const_value(node[2]), _const_value(node[3])
        if node[1] == "*" and left is not None:
            c, d, o = _linear(node[3], text)
            return left * c, {k: left * v for k, v in d.items()}, o
        if right is not None:
            s = right if node[1] == "*" else 1.0 / right
            c, d, o = _linear(node[2], text)
            return s * c, {k: s * v for k, v in d.items()}, o
    raise NSSEquationError(f"exponent of 10 is not linear in the variables in {text!r}")


def _loglin(node: Node, text: str) -> _LogLin:
    """log10 of ``node`` in the representable form, or raise."""
    kind = node[0]
    if kind == "num":
        if node[1] <= 0:
            raise NSSEquationError(f"non-positive factor {node[1]} in {text!r}")
        return _LogLin.constant(math.log10(node[1]))
    if kind == "var":
        _check_variable(node[1], text)
        return _LogLin.term(node[1], "log10")
    if kind == "neg":
        raise NSSEquationError(f"negated factor in {text!r}")
    if kind == "call":
        args = node[2]
        if node[1] == "max" and len(args) == 2 and _const_value(args[0]) == 0.0:
            inner = _loglin(args[1], text)
            inner.notes.append(
                "max(0, E) dropped: E is a product of positive terms, so max is the identity"
            )
            return inner
        raise NSSEquationError(f"unsupported function call {node[1]!r} in {text!r}")

    op, left, right = node[1], node[2], node[3]
    if op == "*":
        return _loglin(left, text).plus(_loglin(right, text), text)
    if op == "/":
        return _loglin(left, text).plus(_loglin(right, text).scaled(-1.0), text)
    if op == "^":
        if _const_value(left) == 10.0:
            c, d, o = _linear(right, text)
            out = _LogLin.constant(c)
            for k in o:
                out = out.plus(_LogLin.term(k, "identity", d[k]), text)
            return out
        p = _const_value(right)
        if p is None:
            raise NSSEquationError(f"non-constant exponent on a base other than 10 in {text!r}")
        return _base(left, text).scaled(p)
    # A sum: only "X + 1" is representable (log10_plus1); "E - 1" and the like are not.
    return _base(node, text)


def _base(node: Node, text: str) -> _LogLin:
    """log10 of the base of a power term. Adds the one offset allowed: ``X + 1``."""
    if node[0] == "bin" and node[1] in "+-":
        left, right = node[2], node[3]
        if node[1] == "+" and left[0] == "var" and _const_value(right) == 1.0:
            _check_variable(left[1], text)
            return _LogLin.term(left[1], "log10_plus1")
        if node[1] == "+" and right[0] == "var" and _const_value(left) == 1.0:
            _check_variable(right[1], text)
            return _LogLin.term(right[1], "log10_plus1")
        raise NSSEquationError(
            f"additive offset cannot be represented exactly (only X+1, as "
            f"log10_plus1) in {text!r}"
        )
    out = _loglin(node, text)
    if node[0] == "bin" and node[1] in "*/" and out.const != 0.0 and out.terms:
        out.notes.append("scale inside a power term folded into the intercept as b*log10(scale)")
    return out


def parse_equation(equation: str) -> ParsedEquation:
    """Rewrite an NSS equation string exactly into log-linear form.

    Parameters
    ----------
    equation : str
        The ``equation`` field of an NSS Estimate result.

    Returns
    -------
    ParsedEquation

    Raises
    ------
    NSSEquationError
        If the string is malformed or has a term the log-linear schema cannot
        represent exactly.
    """
    ll = _loglin(parse_expression(equation), equation)
    if not ll.terms:
        raise NSSEquationError(f"equation has no variables: {equation!r}")
    return ParsedEquation(
        equation=equation,
        intercept=ll.const,
        codes=tuple(ll.order),
        transforms=tuple(ll.terms[k][0] for k in ll.order),
        coefficients=tuple(ll.terms[k][1] for k in ll.order),
        notes=tuple(dict.fromkeys(ll.notes)),
    )


def to_regression_equation(
    parsed: ParsedEquation,
    *,
    region_code: str,
    aep: float,
    citation: Citation,
    region_name: str = "",
    limits: Optional[Mapping[str, Tuple[Optional[float], Optional[float]]]] = None,
    units: Optional[Mapping[str, str]] = None,
) -> RegressionEquation:
    """Build a :class:`RegressionEquation` from a parsed NSS string.

    NSS supplies no report table, covariance, model-error variance or site
    count, so the result carries only what NSS does. ``sep_log`` is left unset
    because NSS's ``sep`` field is not confirmed to be the report's log10 SEP.

    Parameters
    ----------
    limits : mapping of code to (min, max), optional
    units : mapping of code to unit string, optional
    """
    limits = limits or {}
    units = units or {}
    variables = tuple(
        Variable(
            code=k,
            transform=t,
            units=units.get(k, ""),
            minimum=limits.get(k, (None, None))[0],
            maximum=limits.get(k, (None, None))[1],
        )
        for k, t in zip(parsed.codes, parsed.transforms)
    )
    return RegressionEquation(
        region_code=region_code,
        region_name=region_name,
        aep=aep,
        intercept=parsed.intercept,
        variables=variables,
        coefficients=parsed.coefficients,
        citation=citation,
    )


# ---------------------------------------------------------------------------
# Statistic codes
# ---------------------------------------------------------------------------

# PK50AEP, PK66_7AEP, PK0_2AEP; MT prefixes (AC, BW, RS) and the 10-character
# truncation NSS applies to some codes (ACPK66_7AE).
_AEP_CODE = re.compile(r"^(?P<prefix>[A-Z]*?)PK(?P<pct>\d+(?:_\d+)?)(?:AEP|AE|A)$")
# Older recurrence-interval codes: PK2, PK100, PK1_5 (the 1.5-year flood).
_RI_CODE = re.compile(r"^(?P<prefix>[A-Z]*?)PK(?P<years>\d+(?:_\d+)?)$")


def statistic_code_to_aep(code: str) -> float:
    """Annual exceedance probability of an NSS peak-flow statistic code.

    ``PK50AEP`` -> 0.5, ``PK66_7AEP`` -> 0.667, ``PK0_2AEP`` -> 0.002. An
    underscore is the decimal point. A recurrence-interval code (``PK100``)
    maps to ``1/T``.

    Raises
    ------
    ValueError
        If the code is not a peak-flow AEP or recurrence-interval code.
    """
    code = code.strip().upper()
    m = _AEP_CODE.match(code)
    if m:
        aep = float(m.group("pct").replace("_", ".")) / 100.0
    else:
        m = _RI_CODE.match(code)
        if not m:
            raise ValueError(f"{code!r} is not an NSS peak-flow statistic code")
        years = float(m.group("years").replace("_", "."))
        if years <= 1.0:
            raise ValueError(f"{code!r}: recurrence interval must exceed 1 year")
        aep = 1.0 / years
    if not 0.0 < aep < 1.0:
        raise ValueError(f"{code!r} gives AEP {aep}, outside (0, 1)")
    return aep


def variables_in(equation: str) -> Sequence[str]:
    """Variable codes an equation string references, in order, without duplicates."""
    return list(dict.fromkeys(_variables(parse_expression(equation))))
