"""
Detect the fractionation method of an uploaded file from its file name.

Allowed methods: HLB, RP, SCX, QMA, MAX.
RP is detected by finding the uppercase sequence "RP" anywhere in the file name
(so "RPS_pH_10.csv" and "RP_pH10.csv" both count as RP).
If exactly two files are RP, the numbers in the file names (e.g. the pH) decide
which one is labelled "RP high pH" and which one "RP low pH".
"""
import os
import re

ALLOWED_METHODS = ["HLB", "RP", "SCX", "QMA", "MAX"]
RP_HIGH = "RP high pH"
RP_LOW = "RP low pH"

# Options offered in the manual-override dropdown
METHOD_OPTIONS = ["HLB", "RP", RP_HIGH, RP_LOW, "SCX", "QMA", "MAX"]
SKIP = "Skip file"

_NUM = r"\d+(?:\.\d+)?"


def _make_pattern(method):
    if method == "RP":
        # RP: the uppercase sequence "RP" anywhere in the name (case-sensitive)
        return re.compile("RP")
    # Other methods: must not be glued to other letters, case-insensitive
    # ("SCX_BSA.csv" matches SCX, "MAXQuant.csv" does not match MAX)
    return re.compile(rf"(?<![A-Za-z]){method}(?![A-Za-z])", re.IGNORECASE)


_PATTERNS = {m: _make_pattern(m) for m in ALLOWED_METHODS}


def _stem(filename):
    return os.path.splitext(os.path.basename(filename))[0]


def detect_methods(filename):
    """Return the list of allowed methods found in the file name."""
    stem = _stem(filename)
    return [m for m, pattern in _PATTERNS.items() if pattern.search(stem)]


def _ph_value(stem):
    """Number written next to 'pH' (pH10, pH_10, 10pH), or None."""
    m = (re.search(rf"pH\s*[_\-]?\s*({_NUM})", stem, re.IGNORECASE)
         or re.search(rf"({_NUM})\s*[_\-]?\s*pH", stem, re.IGNORECASE))
    return float(m.group(1)) if m else None


def _rank_rp_pair(name_a, name_b):
    """
    Decide which of two RP file names is the high-pH one.
    Returns 0 if name_a is higher, 1 if name_b is higher, None if undecidable.

    1) Compare the numbers written next to 'pH'.
    2) Otherwise compare the first number that differs between the two names
       (numbers shared by both names, e.g. a sample id, are ignored).
    """
    a, b = _stem(name_a), _stem(name_b)

    pa, pb = _ph_value(a), _ph_value(b)
    if pa is not None and pb is not None and pa != pb:
        return 0 if pa > pb else 1

    na = [float(x) for x in re.findall(_NUM, a)]
    nb = [float(x) for x in re.findall(_NUM, b)]
    for x, y in zip(na, nb):
        if x != y:
            return 0 if x > y else 1
    return None


def auto_assign_methods(filenames):
    """
    Returns (methods, notes):
      methods : list aligned with filenames; a value from METHOD_OPTIONS or None
                when no single method could be detected.
      notes   : human-readable messages about anything that needs attention.
    """
    methods, notes = [], []
    for fn in filenames:
        found = detect_methods(fn)
        if len(found) == 1:
            methods.append(found[0])
        else:
            methods.append(None)
            if not found:
                notes.append(f"No allowed method ({', '.join(ALLOWED_METHODS)}) found in '{fn}'.")
            else:
                notes.append(f"'{fn}' matches several methods ({', '.join(found)}) - please pick one.")

    rp_idx = [i for i, m in enumerate(methods) if m == "RP"]
    if len(rp_idx) == 2:
        i, j = rp_idx
        higher = _rank_rp_pair(filenames[i], filenames[j])
        if higher is None:
            notes.append("Two RP files found, but the numbers in their names don't show which "
                         "is high/low pH - please choose manually.")
        else:
            hi, lo = (i, j) if higher == 0 else (j, i)
            methods[hi], methods[lo] = RP_HIGH, RP_LOW
    elif len(rp_idx) > 2:
        notes.append(f"{len(rp_idx)} RP files found - high/low pH is only assigned automatically "
                     f"when there are exactly 2. Choose manually if needed.")
    return methods, notes


def make_unique_labels(methods):
    """['SCX', 'RP', 'SCX'] -> ['SCX (1)', 'RP', 'SCX (2)'] so every file has a distinct label."""
    totals = {}
    for m in methods:
        totals[m] = totals.get(m, 0) + 1
    seen, labels = {}, []
    for m in methods:
        if totals[m] > 1:
            seen[m] = seen.get(m, 0) + 1
            labels.append(f"{m} ({seen[m]})")
        else:
            labels.append(m)
    return labels