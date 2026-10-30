"""Parse NSS regression-equation strings into :class:`RegressionEquation` form.

NSS (``/nssservices/Scenarios/Estimate``) returns each estimate with the literal
equation it evaluated, for example::

    3.846*DRNAREA^0.745*10^(0.032*PRECPRIS10)/10^(0.0078*CANOPY_PCT)
    27.9*DRNAREA^0.523*(MINBELEV/1000)^(-0.098)
    0.00953*(DRNAREA)^(0.392)*(0.001*(RELIEF))^(3.36)
    5.0*CONTDA^0.674*(EL5000+1)^(-0.026)

This module turns such a string into the log-linear form of
:mod:`flowfreq.regression.equations`::

    log10(Q) = b0 + sum_i b_i * T_i(scale_i * X_i + offset_i),
    T_i in {identity, log10, log10_plus1, power}

A term maps as follows:

- ``X^b`` or ``(X)^(b)`` becomes ``log10`` with coefficient ``b``.
- ``(X+1)^b`` becomes ``log10_plus1``.
- Any other affine base in one variable, ``(X+c)^b``, ``(X-c)^b``,
  ``(X/k+c)^b`` or ``(s*X+c)^b``, becomes ``log10`` with ``scale`` and
  ``offset`` (:class:`~flowfreq.regression.equations.Variable`), stored as
  written. A negated base, ``(-1*(X-c))^b`` or ``(c-X)^b``, has ``scale=-1``.
- ``10^(c*X)`` becomes ``identity`` with ``c``, and ``/10^(c*X)`` gives ``-c``.
- ``10^(a - b*X^(k))``, a constant power of a variable inside the exponent of
  10 (Arizona region 5), becomes ``power`` with ``exponent=k`` and coefficient
  ``-b``; ``a`` folds into ``b0``.
- Constant factors, ``10^(k)`` included, fold into ``b0``.
- A scale inside a power term, ``(s*X)^b`` or ``(X/s)^b``, is algebraically
  ``s^b * X^b``. It folds into the intercept as ``b*log10(s)``, so the equation
  stays exact.

``max(0, E)`` is accepted when ``E`` is itself a product of such terms, because
that product is strictly positive wherever its terms are defined.

**Anything else is refused, never approximated.** That covers a subtracted
constant outside a power term (``max(0, E-1)``), a base with two variables or a
negative scale, logistic ``e^(...)`` forms, and a variable that appears under
two different transforms (or the same transform with two different rescalings). Such a string raises
:class:`NSSEquationError`. :func:`evaluate_expression` still evaluates it
directly, so a snapshot can record what NSS computes without claiming to
represent it.

One non-log-linear form is recognised separately. NSS writes a printed constant
``Q=0`` (Nevada region 6's 2-year equation, WSP 2433 Table 10) as a product with a
literal zero factor, ``0*(0*(DRNAREA))^(0)*(0*(ELEV))^(0)``, which
:func:`parse_equation` refuses (``non-positive factor``).
:func:`parse_constant_equation` reads it as a :class:`ParsedConstant`, and
:func:`to_regression_equation` turns that into ``form="constant"`` with
``constant_cfs=0``, its variables kept for the range check only.

Roadmap: ``docs/MASTER_ROADMAP.md`` §3.1, issue #35.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from typing import Dict, List, Mapping, Optional, Sequence, Tuple, Union

from flowfreq.regression.equations import Citation, RegressionEquation, Variable

__all__ = [
    "NSSEquationError",
    "ParsedConstant",
    "ParsedEquation",
    "evaluate_expression",
    "parse_constant_equation",
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
# NSS codes may carry lower-case letters (CSLBlue, CSL10_85fm, SOILCorD).
_VARIABLE = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")
_FUNCTIONS = frozenset({"max"})
_LOG_FUNCTIONS = frozenset({"logN", "LogN"})
# Names that are never a variable code: functions, and the constant e and
# function exp that _eval knows.
_RESERVED = _FUNCTIONS | _LOG_FUNCTIONS | frozenset({"e", "exp"})


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
        if node[1] in _LOG_FUNCTIONS and len(args) == 2:
            return math.log(args[0], args[1])
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
    scales, offsets : tuple of float
        Per-variable :attr:`Variable.scale` and :attr:`Variable.offset`. Empty
        means all defaults (1 and 0).
    exponents : tuple of float or None
        Per-variable :attr:`Variable.exponent`, set only for ``power``. Empty
        means none.
    times : tuple of str or None
        Per-variable :attr:`Variable.times`, set only for an interaction term (a
        power whose exponent is affine in other characteristics). Empty means none.
    notes : tuple of str
        Rewrites a reviewer should know about, such as a scale folded into ``b0``.
    """

    equation: str
    intercept: float
    codes: Tuple[str, ...]
    transforms: Tuple[str, ...]
    coefficients: Tuple[float, ...]
    notes: Tuple[str, ...] = ()
    scales: Tuple[float, ...] = ()
    offsets: Tuple[float, ...] = ()
    exponents: Tuple[Optional[float], ...] = ()
    times: Tuple[Optional[str], ...] = ()

    def variables(self) -> Tuple[Variable, ...]:
        """The variables as :class:`Variable` objects, without limits or units."""
        n = len(self.codes)
        scales = self.scales or (1.0,) * n
        offsets = self.offsets or (0.0,) * n
        exponents: Tuple[Optional[float], ...] = self.exponents or (None,) * n
        times: Tuple[Optional[str], ...] = self.times or (None,) * n
        return tuple(
            Variable(code=k, transform=t, scale=sc, offset=off, exponent=ex, times=y)
            for k, t, sc, off, ex, y in zip(
                self.codes, self.transforms, scales, offsets, exponents, times
            )
        )

    def log10_value(self, values: Mapping[str, float]) -> float:
        """Evaluate ``log10(Q)`` at the given characteristics."""
        return self.intercept + sum(
            c * v.value(values) for v, c in zip(self.variables(), self.coefficients)
        )

    def to_dict(self) -> Dict[str, object]:
        return {
            "intercept": self.intercept,
            "variables": [v.to_dict() for v in self.variables()],
            "coefficients": list(self.coefficients),
            "notes": list(self.notes),
        }


# A term's form: (transform, scale, offset, exponent).
_Form = Tuple[str, float, float, Optional[float]]
# A term's key: (code, times), times set only on an interaction term.
_Key = Tuple[str, Optional[str]]


def _key_text(key: _Key) -> str:
    return key[0] if key[1] is None else f"{key[0]} (times {key[1]})"


class _LogLin:
    """log10 of an expression: ``const + sum coef * T(scale*X + offset) [* Y]``."""

    def __init__(self) -> None:
        self.const = 0.0
        self.terms: Dict[_Key, Tuple[_Form, float]] = {}
        self.order: List[_Key] = []
        self.notes: List[str] = []

    @classmethod
    def constant(cls, value: float) -> "_LogLin":
        out = cls()
        out.const = value
        return out

    @classmethod
    def term(
        cls,
        code: str,
        transform: str,
        coef: float = 1.0,
        scale: float = 1.0,
        offset: float = 0.0,
        exponent: Optional[float] = None,
        times: Optional[str] = None,
    ) -> "_LogLin":
        out = cls()
        out.terms[(code, times)] = ((transform, scale, offset, exponent), coef)
        out.order.append((code, times))
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
                        f"{_key_text(k)} appears under two transforms ({_form_text(t0)}, "
                        f"{_form_text(t)}) in {text!r}; RegressionEquation allows one "
                        "term per variable"
                    )
                out.terms[k] = (t, c0 + c)
            else:
                out.terms[k] = (t, c)
                out.order.append(k)
        out.notes += other.notes
        return out


def _form_text(form: _Form) -> str:
    t, scale, offset, exponent = form
    if exponent is not None:
        t = f"{t}^{exponent:g}"
    if scale == 1.0 and offset == 0.0:
        return t
    return f"{t}(scale={scale:g}, offset={offset:g})"


def _const_value(node: Node) -> Optional[float]:
    """The value of a numeric literal, possibly negated; None if not constant."""
    if node[0] == "num":
        return float(node[1])
    if node[0] == "neg":
        inner = _const_value(node[1])
        return None if inner is None else -inner
    return None


def _check_variable(name: str, text: str) -> None:
    if not _VARIABLE.match(name) or name in _RESERVED:
        raise NSSEquationError(f"unsupported identifier {name!r} in {text!r}")


def _const_product(node: Node) -> Optional[float]:
    """The value of a literal or a product/quotient of literals; None otherwise."""
    value = _const_value(node)
    if value is not None or node[0] != "bin" or node[1] not in "*/":
        return value
    left, right = _const_product(node[2]), _const_product(node[3])
    if left is None or right is None or (node[1] == "/" and right == 0.0):
        return None
    return left * right if node[1] == "*" else left / right


_LogForm = Tuple[str, float, float]


@dataclass
class _Linear:
    """``_linear``'s result: ``const + sum coefs[X] * T(X)`` over ``order``.

    ``T(X)`` is ``X`` itself unless ``X`` enters as a power ``X^k`` (``powers[X]``
    is ``k``) or as a logarithm (``logs[X]`` is the ``(transform, scale, offset)``
    of a log10 term, as :class:`Variable` defines them).
    """

    const: float = 0.0
    coefs: Dict[str, float] = field(default_factory=dict)
    order: List[str] = field(default_factory=list)
    powers: Dict[str, float] = field(default_factory=dict)
    logs: Dict[str, _LogForm] = field(default_factory=dict)
    notes: List[str] = field(default_factory=list)

    def is_identity(self, code: str) -> bool:
        return code not in self.powers and code not in self.logs

    def scaled(self, s: float) -> "_Linear":
        return _Linear(
            self.const * s,
            {k: s * v for k, v in self.coefs.items()},
            list(self.order),
            dict(self.powers),
            dict(self.logs),
            list(self.notes),
        )

    def plus(self, other: "_Linear", sign: float, text: str) -> "_Linear":
        out = self.scaled(1.0)
        out.const += sign * other.const
        for k in other.order:
            if k in out.coefs:
                if (out.powers.get(k), out.logs.get(k)) != (
                    other.powers.get(k),
                    other.logs.get(k),
                ):
                    raise NSSEquationError(
                        f"{k} appears in two different forms in an exponent of 10 in {text!r}"
                    )
                out.coefs[k] += sign * other.coefs[k]
            else:
                out.coefs[k] = sign * other.coefs[k]
                out.order.append(k)
        out.powers.update(other.powers)
        out.logs.update(other.logs)
        out.notes += other.notes
        return out


def _log_term(node: Node, text: str) -> _Linear:
    """NSS's ``logN(E, b)``: ``log10(E) / log10(b)``, with ``E`` read as :func:`_base` does."""
    args = node[2]
    b = _const_value(args[1]) if len(args) == 2 else None
    if b is None or not b > 0 or b == 1.0:
        raise NSSEquationError(
            f"{node[1]} needs two arguments, the second a positive constant base other "
            f"than 1, in {text!r}"
        )
    inner = _base(args[0], text)
    if any(k[1] is not None for k in inner.order):
        raise NSSEquationError(f"interaction term inside {node[1]} in {text!r}")
    s = 1.0 / math.log10(b)
    out = _Linear(const=inner.const * s, notes=list(inner.notes))
    if b != 10.0:
        out.notes.append("logN(E, b) rewritten as log10(E)/log10(b)")
    for key in inner.order:
        (transform, scale, offset, exponent), coef = inner.terms[key]
        k = key[0]
        out.coefs[k] = coef * s
        out.order.append(k)
        if transform == "power" and exponent is not None:
            out.powers[k] = exponent
        elif transform != "identity":
            out.logs[k] = (transform, scale, offset)
    return out


def _linear(node: Node, text: str) -> _Linear:
    """``node`` as ``const + sum c*T(X)`` (the exponent of a ``10^(...)`` term).

    ``T(X)`` is ``X``, a power ``X^k``, or ``logN(E, b)`` on an affine base ``E``.
    """
    kind = node[0]
    if kind == "num":
        return _Linear(const=float(node[1]))
    if kind == "var":
        _check_variable(node[1], text)
        return _Linear(coefs={node[1]: 1.0}, order=[node[1]])
    if kind == "neg":
        return _linear(node[1], text).scaled(-1.0)
    if kind == "call" and node[1] in _LOG_FUNCTIONS:
        return _log_term(node, text)
    if kind == "bin" and node[1] in "+-":
        sign = 1.0 if node[1] == "+" else -1.0
        return _linear(node[2], text).plus(_linear(node[3], text), sign, text)
    if kind == "bin" and node[1] in "*/":
        left, right = _const_product(node[2]), _const_product(node[3])
        if node[1] == "*" and left is not None:
            s, factor, operand = left, node[2], node[3]
        elif right is not None and (node[1] == "*" or right != 0.0):
            s = right if node[1] == "*" else 1.0 / right
            factor, operand = node[3], node[2]
        else:
            raise NSSEquationError(f"exponent of 10 is not linear in the variables in {text!r}")
        inner = _linear(operand, text)
        identity = [k for k in inner.order if inner.is_identity(k)]
        if identity and (
            _const_value(factor) is None or any(abs(inner.coefs[k]) != 1.0 for k in identity)
        ):
            inner.notes.append(
                "constant factors folded into the coefficient of " + ", ".join(identity)
            )
        return inner.scaled(s)
    if kind == "bin" and node[1] == "^" and node[2][0] == "var":
        p = _const_value(node[3])
        if p is not None and p != 0.0:
            code = node[2][1]
            _check_variable(code, text)
            powers = {} if p == 1.0 else {code: p}
            return _Linear(coefs={code: 1.0}, order=[code], powers=powers)
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
            lin = _linear(right, text)
            out = _LogLin.constant(lin.const)
            out.notes = list(lin.notes)
            for k in lin.order:
                if k in lin.logs:
                    transform, scale, offset = lin.logs[k]
                    term = _LogLin.term(k, transform, lin.coefs[k], scale, offset)
                elif k in lin.powers:
                    term = _LogLin.term(k, "power", lin.coefs[k], exponent=lin.powers[k])
                else:
                    term = _LogLin.term(k, "identity", lin.coefs[k])
                out = out.plus(term, text)
            return out
        p = _const_value(right)
        if p is None:
            return _affine_power(left, right, text)
        return _base(left, text).scaled(p)
    # A bare sum: representable only as an affine base in one variable.
    return _base(node, text)


def _affine_power(base_node: Node, exponent_node: Node, text: str) -> _LogLin:
    """log10 of ``E^(c0 + sum c_i*Y_i)``: ``c0*log10(E) + sum c_i*Y_i*log10(E)``.

    The exponent must be affine in plain characteristics (no powers or logs of
    them). Each ``Y_i`` term becomes one interaction term per term of ``log10(E)``,
    with :attr:`Variable.times` set to ``Y_i``; a constant factor ``k`` inside ``E``
    becomes an identity term ``c_i*log10(k)*Y_i``.
    """
    refusal = f"non-constant exponent on a base other than 10 in {text!r}"
    try:
        lin = _linear(exponent_node, text)
    except NSSEquationError:
        raise NSSEquationError(refusal) from None
    if lin.powers or lin.logs:
        raise NSSEquationError(refusal + "; the exponent is not affine in the variables")
    base = _base(base_node, text)
    if any(k[1] is not None for k in base.order):
        raise NSSEquationError(refusal + "; the base already carries an interaction")
    out = base.scaled(lin.const)
    for y in lin.order:
        c = lin.coefs[y]
        if base.const != 0.0:
            out = out.plus(_LogLin.term(y, "identity", c * base.const), text)
        for key in base.order:
            (transform, scale, offset, exponent), coef = base.terms[key]
            term = _LogLin.term(key[0], transform, coef * c, scale, offset, exponent, times=y)
            out = out.plus(term, text)
    out.notes.append("exponent affine in characteristics stored as interaction terms (times)")
    return out


def _negated(node: Node) -> bool:
    """Whether ``node`` is ``-E`` or ``E`` times or over a negative constant."""
    if node[0] == "neg":
        return True
    if node[0] != "bin" or node[1] not in "*/":
        return False
    left, right = _const_product(node[2]), _const_product(node[3])
    return (left is not None and left < 0) or (right is not None and right < 0)


def _base(node: Node, text: str) -> _LogLin:
    """log10 of the base of a power term, including an affine base ``s*X + c``.

    ``X + 1`` keeps its historical form, ``log10_plus1``. Any other sum in one
    variable becomes ``log10`` with ``scale=s`` and ``offset=c``, as written. A
    negated base, ``-E`` or ``-k*E``, is affine too: ``(-1*(X-13))`` and
    ``(13-X)`` both give ``scale=-1``, ``offset=13``.
    """
    negated = _negated(node)
    if negated or (node[0] == "bin" and node[1] in "+-"):
        refusal = (
            f"negated base that is not affine in one variable in {text!r}"
            if negated
            else f"additive offset in a base that is not affine in one variable in {text!r}"
        )
        try:
            lin = _linear(node, text)
        except NSSEquationError:
            raise NSSEquationError(refusal) from None
        c, d, order = lin.const, lin.coefs, lin.order
        if lin.powers or lin.logs:
            raise NSSEquationError(refusal)
        if len(order) != 1:
            raise NSSEquationError(
                f"additive offset in a base with {len(order)} variables in {text!r}; "
                "only an affine base in one variable is representable"
            )
        code = order[0]
        scale = d[code]
        if scale == 0.0:
            raise NSSEquationError(f"zero scale on {code} in {text!r}")
        if scale == 1.0 and c == 1.0:
            return _LogLin.term(code, "log10_plus1")
        return _LogLin.term(code, "log10", scale=scale, offset=c)
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
    times = tuple(k[1] for k in ll.order)
    return ParsedEquation(
        equation=equation,
        intercept=ll.const,
        codes=tuple(k[0] for k in ll.order),
        transforms=tuple(ll.terms[k][0][0] for k in ll.order),
        coefficients=tuple(ll.terms[k][1] for k in ll.order),
        notes=tuple(dict.fromkeys(ll.notes)),
        scales=tuple(ll.terms[k][0][1] for k in ll.order),
        offsets=tuple(ll.terms[k][0][2] for k in ll.order),
        exponents=tuple(ll.terms[k][0][3] for k in ll.order),
        times=times if any(y is not None for y in times) else (),
    )


@dataclass(frozen=True)
class ParsedConstant:
    """An NSS equation that is a constant flow, such as NSS's form of a printed ``Q=0``.

    Attributes
    ----------
    equation : str
        The NSS string it came from.
    value : float
        The constant flow, cfs.
    codes : tuple of str
        The variables the string names, in order of first appearance. They carry the
        region's limits for the range check and no coefficient.
    """

    equation: str
    value: float
    codes: Tuple[str, ...]

    def variables(self) -> Tuple[Variable, ...]:
        """The variables as :class:`Variable` objects, without limits or units."""
        return tuple(Variable(code=k, transform="log10") for k in self.codes)

    def to_dict(self) -> Dict[str, object]:
        return {
            "form": "constant",
            "constant_cfs": self.value,
            "variables": [v.to_dict() for v in self.variables()],
        }


def _factors(node: Node) -> List[Node]:
    """The factors of a product ``a*b*...``; any other node is one factor."""
    if node[0] == "bin" and node[1] == "*":
        return _factors(node[2]) + _factors(node[3])
    return [node]


def parse_constant_equation(equation: str) -> ParsedConstant:
    """Read an NSS string that is a product with a literal zero factor as ``Q=0``.

    Only that form is accepted: the top-level expression must be a product (``*``
    only, no division) one of whose factors is the literal ``0``. Its other factors are
    ignored, as in ``0*(0*(DRNAREA))^(0)*(0*(ELEV))^(0)``, so the result is exactly 0
    wherever NSS's own evaluation is defined.

    Raises
    ------
    NSSEquationError
        If the string is malformed, is not such a product, or names an unsupported
        identifier.
    """
    node = parse_expression(equation)
    factors = _factors(node)
    if len(factors) < 2 or not any(_const_value(f) == 0.0 for f in factors):
        raise NSSEquationError(f"not a product with a literal zero factor: {equation!r}")
    codes = list(dict.fromkeys(_variables(node)))
    for code in codes:
        _check_variable(code, equation)
    return ParsedConstant(equation=equation, value=0.0, codes=tuple(codes))


def to_regression_equation(
    parsed: Union[ParsedEquation, ParsedConstant],
    *,
    region_code: str,
    aep: float,
    citation: Citation,
    region_name: str = "",
    limits: Optional[Mapping[str, Tuple[Optional[float], Optional[float]]]] = None,
    units: Optional[Mapping[str, str]] = None,
) -> RegressionEquation:
    """Build a :class:`RegressionEquation` from a parsed NSS string.

    A :class:`ParsedConstant` gives ``form="constant"`` with zero coefficients.
    NSS supplies no report table, covariance, model-error variance or site
    count, so the result carries only what NSS does. The ``times`` characteristic of
    an interaction term that is not also a variable gets its limits as an
    ``applicability`` limit. ``sep_log`` is left unset
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
            code=v.code,
            transform=v.transform,
            units=units.get(v.code, ""),
            minimum=limits.get(v.code, (None, None))[0],
            maximum=limits.get(v.code, (None, None))[1],
            scale=v.scale,
            offset=v.offset,
            exponent=v.exponent,
            times=v.times,
        )
        for v in parsed.variables()
    )
    if isinstance(parsed, ParsedConstant):
        return RegressionEquation(
            region_code=region_code,
            region_name=region_name,
            aep=aep,
            intercept=None,
            variables=variables,
            coefficients=(0.0,) * len(variables),
            citation=citation,
            form="constant",
            constant_cfs=parsed.value,
        )
    # An interaction's ``times`` characteristic that is not a variable in its own right
    # keeps its published range as an applicability limit, so it is still checked.
    codes = {v.code for v in variables}
    times_codes = dict.fromkeys(v.times for v in variables if v.times is not None)
    applicability = tuple(
        Variable(
            code=y,
            transform="identity",
            units=units.get(y, ""),
            minimum=limits[y][0],
            maximum=limits[y][1],
        )
        for y in times_codes
        if y not in codes and y in limits and limits[y] != (None, None)
    )
    return RegressionEquation(
        region_code=region_code,
        region_name=region_name,
        aep=aep,
        intercept=parsed.intercept,
        variables=variables,
        coefficients=parsed.coefficients,
        citation=citation,
        applicability=applicability,
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
