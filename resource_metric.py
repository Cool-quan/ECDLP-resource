from dataclasses import dataclass

import numpy as np


METRIC_NAMES = ("Toffoli depth", "Toffoli count", "CNOT depth", "CNOT count")

ADDITION = "addition"
COMPARISON = "comparison"
TABLE = "table"
MISC = "misc"
DELTA = "delta"

CATEGORIES = (ADDITION, COMPARISON, TABLE, MISC, DELTA)


def _as_values(values):
    arr = np.array(values, dtype=object)
    if arr.shape != (4,):
        raise ValueError("Resource metrics must have four entries.")
    return arr


def _zero_parts():
    return {category: np.array([0, 0, 0, 0], dtype=object) for category in CATEGORIES}


@dataclass(frozen=True)
class ResourceMetric:
    """
    Four circuit metrics plus per-category contributions.

    The value order follows the original scripts:
    Toffoli depth, Toffoli count, CNOT depth, CNOT count.
    """

    values: np.ndarray
    parts: dict

    __array_priority__ = 1000

    def __init__(self, values, parts=None):
        object.__setattr__(self, "values", _as_values(values))
        normalized = _zero_parts()
        if parts is not None:
            for category, contribution in parts.items():
                if category not in CATEGORIES:
                    raise ValueError(f"Unknown metric category: {category}")
                normalized[category] = _as_values(contribution)
        object.__setattr__(self, "parts", normalized)

    @classmethod
    def zero(cls):
        return cls([0, 0, 0, 0])

    @classmethod
    def basic(cls, values, category):
        parts = _zero_parts()
        parts[category] = _as_values(values)
        return cls(values, parts)

    @classmethod
    def composite(cls, terms):
        total = cls.zero()
        for term in terms:
            total = total + term
        return total

    @classmethod
    def delta(cls, values):
        return cls.basic(values, DELTA)

    @classmethod
    def from_axes(cls, axis_sources):
        if len(axis_sources) != 4:
            raise ValueError("Exactly four axis sources are required.")

        values = np.array([0, 0, 0, 0], dtype=object)
        parts = _zero_parts()
        for axis, source in enumerate(axis_sources):
            metric = ensure_metric(source)
            values[axis] = metric.values[axis]
            for category in CATEGORIES:
                parts[category][axis] = metric.parts[category][axis]
        return cls(values, parts)

    def with_cnot(self, cnot):
        """Keep the Toffoli metrics and take the CNOT metrics from cnot."""
        return self.from_axes([self, self, cnot, cnot])

    def reclassify(self, category):
        return ResourceMetric.basic(self.values, category)

    def only(self, axis):
        values = np.array([0, 0, 0, 0], dtype=object)
        values[axis] = self.values[axis]
        parts = _zero_parts()
        for category in CATEGORIES:
            parts[category][axis] = self.parts[category][axis]
        return ResourceMetric(values, parts)
    
    def only_value(self, axis):
        return self.values[axis]

    def to_array(self):
        return self.values.copy()

    def tolist(self):
        return self.values.tolist()

    def proportions(self):
        result = {}
        for axis, name in enumerate(METRIC_NAMES):
            total = self.values[axis]
            if total == 0:
                result[name] = {category: 0.0 for category in CATEGORIES}
            else:
                result[name] = {
                    category: float(self.parts[category][axis] / total)
                    for category in CATEGORIES
                }
        return result

    def ratio_row(self):
        proportions = self.proportions()
        return [
            proportions[metric_name][category]
            for metric_name in METRIC_NAMES
            for category in CATEGORIES
        ]

    @staticmethod
    def ratio_headers():
        return [
            f"{metric_name} {category} ratio"
            for metric_name in METRIC_NAMES
            for category in CATEGORIES
        ]

    def csv_row(self):
        return self.tolist() + self.ratio_row()

    def format_summary(self, digits=2):
        value_text = ", ".join(
            f"{name}={self.values[i]}" for i, name in enumerate(METRIC_NAMES)
        )
        pieces = []
        proportions = self.proportions()
        for metric_name in METRIC_NAMES:
            ratio_text = ", ".join(
                f"{category}={proportions[metric_name][category] * 100:.{digits}f}%"
                for category in CATEGORIES
            )
            pieces.append(f"{metric_name}: {ratio_text}")
        return f"{value_text}\n  " + "\n  ".join(pieces)

    def _binary(self, other, value_op, part_op):
        other_metric = ensure_metric(other)
        values = value_op(self.values, other_metric.values)
        parts = {
            category: part_op(self.parts[category], other_metric.parts[category])
            for category in CATEGORIES
        }
        return ResourceMetric(values, parts)

    def __add__(self, other):
        return self._binary(other, lambda a, b: a + b, lambda a, b: a + b)

    def __radd__(self, other):
        return ensure_metric(other).__add__(self)

    def __sub__(self, other):
        return self._binary(other, lambda a, b: a - b, lambda a, b: a - b)

    def __rsub__(self, other):
        return ensure_metric(other).__sub__(self)

    def __neg__(self):
        return ResourceMetric(
            -self.values,
            {category: -self.parts[category] for category in CATEGORIES},
        )

    def __mul__(self, other):
        if not isinstance(other, (int, float, np.integer, np.floating)):
            return NotImplemented
        return ResourceMetric(
            self.values * other,
            {category: self.parts[category] * other for category in CATEGORIES},
        )

    def __rmul__(self, other):
        return self.__mul__(other)

    def __floordiv__(self, other):
        if not isinstance(other, (int, np.integer)):
            return NotImplemented

        values = self.values // other
        parts = {
            category: self.parts[category] // other
            for category in CATEGORIES
        }
        remainder = values - sum(parts.values(), np.array([0, 0, 0, 0], dtype=object))
        parts[DELTA] = parts[DELTA] + remainder
        return ResourceMetric(values, parts)

    def __getitem__(self, key):
        return self.values[key]

    def __iter__(self):
        return iter(self.values)

    def __len__(self):
        return len(self.values)

    def __array__(self, dtype=None):
        if dtype is None:
            return np.asarray(self.values)
        return np.asarray(self.values, dtype=dtype)

    def __repr__(self):
        return f"ResourceMetric({self.values.tolist()})"


def ensure_metric(value):
    if isinstance(value, ResourceMetric):
        return value
    if isinstance(value, (int, float, np.integer, np.floating)) and value == 0:
        return ResourceMetric.zero()
    return ResourceMetric.delta(value)


def metric(values, category=DELTA):
    return ResourceMetric.basic(values, category)


def delta(values):
    return ResourceMetric.delta(values)


def zero():
    return ResourceMetric.zero()


def axis_max(left, right, axis):
    left_metric = ensure_metric(left)
    right_metric = ensure_metric(right)
    if left_metric[axis] >= right_metric[axis]:
        return left_metric.only(axis)
    return right_metric.only(axis)
