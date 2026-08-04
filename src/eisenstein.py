"""Exact arithmetic in Z[omega], where omega^2 + omega + 1 = 0."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, order=True)
class Eisenstein:
    """Represent ``a + b*omega`` by the integer pair ``(a,b)``."""

    a: int = 0
    b: int = 0

    def __add__(self, other: object) -> "Eisenstein":
        if isinstance(other, int):
            return Eisenstein(self.a + other, self.b)
        if isinstance(other, Eisenstein):
            return Eisenstein(self.a + other.a, self.b + other.b)
        return NotImplemented

    __radd__ = __add__

    def __sub__(self, other: object) -> "Eisenstein":
        if isinstance(other, int):
            return Eisenstein(self.a - other, self.b)
        if isinstance(other, Eisenstein):
            return Eisenstein(self.a - other.a, self.b - other.b)
        return NotImplemented

    def __rsub__(self, other: object) -> "Eisenstein":
        if isinstance(other, int):
            return Eisenstein(other - self.a, -self.b)
        if isinstance(other, Eisenstein):
            return other - self
        return NotImplemented

    def __neg__(self) -> "Eisenstein":
        return Eisenstein(-self.a, -self.b)

    def __mul__(self, other: object) -> "Eisenstein":
        if isinstance(other, int):
            return Eisenstein(self.a * other, self.b * other)
        if isinstance(other, Eisenstein):
            # (a+bw)(c+dw)=(ac-bd)+(ad+bc-bd)w because w^2=-1-w.
            return Eisenstein(
                self.a * other.a - self.b * other.b,
                self.a * other.b
                + self.b * other.a
                - self.b * other.b,
            )
        return NotImplemented

    __rmul__ = __mul__

    def conjugate(self) -> "Eisenstein":
        # conjugate(w)=w^2=-1-w.
        return Eisenstein(self.a - self.b, -self.b)

    def norm(self) -> int:
        return self.a * self.a - self.a * self.b + self.b * self.b

    def as_pair(self) -> tuple[int, int]:
        return self.a, self.b


ZERO = Eisenstein(0, 0)
ONE = Eisenstein(1, 0)
OMEGA = Eisenstein(0, 1)
OMEGA2 = Eisenstein(-1, -1)
ZETA6 = Eisenstein(1, 1)
ONE_MINUS_OMEGA = Eisenstein(1, -1)

# zeta_6^k for k=0,...,5.
SIXTH_ROOTS: tuple[Eisenstein, ...] = (
    Eisenstein(1, 0),
    Eisenstein(1, 1),
    Eisenstein(0, 1),
    Eisenstein(-1, 0),
    Eisenstein(-1, -1),
    Eisenstein(0, -1),
)
