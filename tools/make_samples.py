"""Generate the sample .xyz files in samples/basics/ from idealised geometries.

Run with:  npm run samples   (or  python tools/make_samples.py)
No third-party packages needed.
"""

import json
import math
from pathlib import Path

OUT = Path(__file__).resolve().parent.parent / "samples" / "basics"

# ---------------------------------------------------------------------------
# tiny vector helpers
# ---------------------------------------------------------------------------


def add(a, b): return tuple(x + y for x, y in zip(a, b))
def sub(a, b): return tuple(x - y for x, y in zip(a, b))
def mul(a, s): return tuple(x * s for x in a)
def dot(a, b): return sum(x * y for x, y in zip(a, b))
def norm(a): return math.sqrt(dot(a, a))
def unit(a): return mul(a, 1 / norm(a))


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def zmat(a, b, c, bond, angle, dihedral):
    """Place atom D bonded to c, with angle b-c-D and dihedral a-b-c-D (degrees). NeRF method."""
    ang, dih = math.radians(angle), math.radians(dihedral)
    bc = unit(sub(c, b))
    n = unit(cross(sub(b, a), bc))
    m = cross(n, bc)
    d = (-bond * math.cos(ang), bond * math.sin(ang) * math.cos(dih), bond * math.sin(ang) * math.sin(dih))
    return add(c, add(mul(bc, d[0]), add(mul(m, d[1]), mul(n, d[2]))))


def two_h(c, n1, n2, length=1.09):
    """Two tetrahedral H's on atom c, which is already bonded to n1 and n2."""
    u1, u2 = unit(sub(n1, c)), unit(sub(n2, c))
    bis = unit(mul(add(u1, u2), -1))
    perp = unit(cross(u1, u2))
    half = math.radians(109.47 / 2)
    return [add(c, mul(add(mul(bis, math.cos(half)), mul(perp, s * math.sin(half))), length)) for s in (1, -1)]


# ---------------------------------------------------------------------------
# molecules: each returns [(symbol, (x, y, z)), ...]
# ---------------------------------------------------------------------------


def water():
    a = math.radians(104.5 / 2)
    return [("O", (0, 0, 0)), ("H", (0.9572 * math.sin(a), 0.9572 * math.cos(a), 0)),
            ("H", (-0.9572 * math.sin(a), 0.9572 * math.cos(a), 0))]


def methane():
    s = 1.09 / math.sqrt(3)
    return [("C", (0, 0, 0))] + [("H", (s * x, s * y, s * z)) for x, y, z in
                                 [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]]


def ethanol():
    c1 = (0.0, 0.0, 0.0)
    c2 = (1.52, 0.0, 0.0)
    t = math.radians(180 - 109.5)  # C1–C2–O angle of 109.5°
    o = add(c2, (1.43 * math.cos(t), 1.43 * math.sin(t), 0.0))
    atoms = [("C", c1), ("C", c2), ("O", o)]
    atoms.append(("H", zmat(c1, c2, o, 0.96, 108.5, 180)))
    for d in (180, 60, -60):
        atoms.append(("H", zmat(o, c2, c1, 1.09, 110, d)))
    for d in (120, -120):
        atoms.append(("H", zmat(o, c1, c2, 1.09, 110, d)))
    return atoms


def benzene():
    atoms = []
    for k in range(6):
        t = math.radians(60 * k)
        atoms.append(("C", (1.39 * math.cos(t), 1.39 * math.sin(t), 0)))
    for k in range(6):
        t = math.radians(60 * k)
        atoms.append(("H", (2.48 * math.cos(t), 2.48 * math.sin(t), 0)))
    return atoms


def cyclohexane():
    R = math.sqrt(1.53 ** 2 - 0.5 ** 2)  # puckering ±0.25 Å gives C–C = 1.53 Å
    cs = []
    for k in range(6):
        t = math.radians(60 * k)
        cs.append((R * math.cos(t), R * math.sin(t), 0.25 if k % 2 == 0 else -0.25))
    atoms = [("C", c) for c in cs]
    for k, c in enumerate(cs):
        for h in two_h(c, cs[k - 1], cs[(k + 1) % 6]):
            atoms.append(("H", h))
    return atoms


def sf6():
    d = 1.564
    return [("S", (0, 0, 0))] + [("F", v) for v in
                                 [(d, 0, 0), (-d, 0, 0), (0, d, 0), (0, -d, 0), (0, 0, d), (0, 0, -d)]]


def ferrocene():
    """Eclipsed (D5h) ferrocene: Fe between two parallel cyclopentadienyl rings."""
    r_c = 1.42 / (2 * math.sin(math.radians(36)))  # ring radius for C–C = 1.42 Å
    r_h = r_c + 1.08
    h = 1.66  # Fe to ring-plane distance
    atoms = [("Fe", (0, 0, 0))]
    for z in (h, -h):
        for k in range(5):
            t = math.radians(72 * k)
            atoms.append(("C", (r_c * math.cos(t), r_c * math.sin(t), z)))
        for k in range(5):
            t = math.radians(72 * k)
            atoms.append(("H", (r_h * math.cos(t), r_h * math.sin(t), z)))
    return atoms


def c60():
    """Buckminsterfullerene: truncated icosahedron, cyclic permutations of three point families."""
    phi = (1 + math.sqrt(5)) / 2
    scale = 1.43 / 2  # the canonical coordinates have edge length 2
    base = []
    for sy in (1, -1):
        for sz in (1, -1):
            base.append((0, sy * 1, sz * 3 * phi))
    for sx in (1, -1):
        for sy in (1, -1):
            for sz in (1, -1):
                base.append((sx * 1, sy * (2 + phi), sz * 2 * phi))
                base.append((sx * phi, sy * 2, sz * (2 * phi + 1)))
    pts = []
    for x, y, z in base:
        for p in ((x, y, z), (y, z, x), (z, x, y)):
            pts.append(mul(p, scale))
    assert len(pts) == 60
    return [("C", p) for p in pts]


def diamond(n=4):
    """A cube of n×n×n diamond unit cells (no hydrogen capping)."""
    a = 3.567
    basis = [(0, 0, 0), (0, .5, .5), (.5, 0, .5), (.5, .5, 0)]
    basis += [add(b, (.25, .25, .25)) for b in basis]
    seen, atoms = set(), []
    for i in range(n + 1):
        for j in range(n + 1):
            for k in range(n + 1):
                for b in basis:
                    f = add((i, j, k), b)
                    if max(f) > n + 1e-9:
                        continue
                    key = tuple(round(x, 4) for x in f)
                    if key not in seen:
                        seen.add(key)
                        atoms.append(("C", mul(f, a)))
    return atoms


def ethane_frame(torsion):
    """Ethane with the back methyl rotated by `torsion` degrees (0 = eclipsed)."""
    cc = 1.535
    ang = math.radians(180 - 111.2)
    r = 1.094 * math.sin(ang)
    dz = 1.094 * math.cos(ang)
    atoms = [("C", (0, 0, cc / 2)), ("C", (0, 0, -cc / 2))]
    for k in range(3):
        t = math.radians(120 * k)
        atoms.append(("H", (r * math.cos(t), r * math.sin(t), cc / 2 + dz)))
    for k in range(3):
        t = math.radians(120 * k + torsion)
        atoms.append(("H", (r * math.cos(t), r * math.sin(t), -cc / 2 - dz)))
    return atoms


# ---------------------------------------------------------------------------
# writing
# ---------------------------------------------------------------------------


def frame_text(atoms, comment):
    lines = [str(len(atoms)), comment]
    lines += [f"{s:<2} {x:12.6f} {y:12.6f} {z:12.6f}" for s, (x, y, z) in atoms]
    return "\n".join(lines) + "\n"


SAMPLES = [
    ("ferrocene.xyz", "Ferrocene", "Fe(C5H5)2, iron sandwiched between two rings", ferrocene),
    ("water.xyz", "Water", "H2O", water),
    ("methane.xyz", "Methane", "CH4, perfect tetrahedron", methane),
    ("ethanol.xyz", "Ethanol", "C2H5OH", ethanol),
    ("benzene.xyz", "Benzene", "C6H6 aromatic ring", benzene),
    ("cyclohexane.xyz", "Cyclohexane (chair)", "C6H12 in the chair conformation", cyclohexane),
    ("sf6.xyz", "Sulfur hexafluoride", "Octahedral SF6", sf6),
    ("c60.xyz", "Buckminsterfullerene", "C60 'buckyball'", c60),
    ("diamond.xyz", "Diamond lattice", "4x4x4 unit cells, ~600 carbon atoms", diamond),
]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    index = []
    for file, name, desc, fn in SAMPLES:
        (OUT / file).write_text(frame_text(fn(), f"{name}: {desc}"), encoding="utf-8")
        index.append({"file": file, "name": name, "description": desc})

    traj = "".join(frame_text(ethane_frame(t), f"Ethane, H-C-C-H torsion {t} deg") for t in range(0, 360, 5))
    (OUT / "ethane-rotation.xyz").write_text(traj, encoding="utf-8")
    index.append({"file": "ethane-rotation.xyz", "name": "Ethane rotation (trajectory)",
                  "description": "72-frame torsion scan around the C-C bond"})

    (OUT / "index.json").write_text(json.dumps(index, indent=2, ensure_ascii=False), encoding="utf-8")
    for f in index:
        print("wrote", f["file"])


if __name__ == "__main__":
    main()
