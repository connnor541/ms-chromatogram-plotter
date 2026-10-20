"""
Detect the fractionation method of an uploaded file from its file name.

Allowed methods: HLB, RP, SCX, QMA, MAX.
methods are detected by finding the respective uppercase name anywhere in the file name
If exactly two files are RP, the numbers in the file names decide
which one is labelled "RP high pH" and which one "RP low pH".
"""
import os
import re

ALLOWED_METHODS = ["HLB", "RP", "SCX", "QMA", "MAX"]
RP_HIGH = "RP high pH"
RP_LOW = "RP low pH"

# Options offered in the manual-override dropdown
METHOD_OPTIONS = ["HLB", "RP", RP_HIGH, RP_LOW, "SCX", "QMA", "MAX"]
#Add option to change method of file to skip file, to not analyze it anymore
SKIP = "Skip file"

#regex pattern to match numbers needed to decide between RP low pH and RP high pH
_NUM = r"\d+(?:\.\d+)?"

#return filename without path and extension
def _stem(filename):
    return os.path.splitext(os.path.basename(filename))[0]

def detect_methods(filename):
    #Returns the list of allowed methods found in the file name.
    stem = _stem(filename)
    return [m for m in ALLOWED_METHODS if m in stem] 
 
def _ph_value(stem):    
    #pH in filename is detected as number written after 'pH' (case-insensitive) or after 'RP'. If no match is found it is set to None
    #note that pH has priority of RP in the searching and first match gets picked
    m = (re.search(rf"pH\s*[_\-]?\s*({_NUM})", stem, re.IGNORECASE)
         or re.search(rf"RP\s*[_\-]?\s*({_NUM})", stem))
    return float(m.group(1)) if m else None
 
def _rank_rp_pair(name_a, name_b):
    #Decides which of two RP file names is the high-pH one. Returns 0 if name_a is higher, 1 if name_b is higher, None if undecidable.
    #There is still a fallback mechanism, if both values are the same or one doesn't exist it will look for other numbers in the filename    
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
    #takes a list of uploaded file names and tries to allocate a fractionation method to each one. It returns two lists:
    #methods: a list of the fractionation methods per file, notes: a list with warnings if things went wrong
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

