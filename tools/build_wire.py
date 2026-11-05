"""Build the repeat unit for the parametric "Molecular wire" object (Triplet molecules tab).

The wire is n 9,10-diphenylanthracene (DPA) units joined anthracene C2 -> 1,4-phenylene -> C6 of the next
anthracene, with H on both ends:

    H-[DPA-2,6-diyl]-(C6H4)-[DPA-2,6-diyl]-(C6H4)- ... -[DPA-2,6-diyl]-H        (n = 1 is plain DPA)

Both anthracene-2,6-diyl and p-phenylene are centrosymmetric, so in a straight wire every DPA + linker
"repeat unit" is the same set of atoms shifted by one translation vector t. This script finds that unit and t:

  1. builds an oligomer (RDKit ETKDGv3 + MMFF94s, lowest of several conformers),
  2. takes the middle DPA + linker as the unit and t as the shift between neighbouring DPAs,
  3. re-assembles a longer wire from the unit, relaxes it with MMFF94s and takes the middle unit again
     (repeated, so the unit is self-consistent with its neighbours on both sides).

The app then lays down n copies (src/parametric.js). This script validates assembled wires the same way
as the other examples (bond-length windows), plus: how far MMFF moves the assembled wire when relaxed
(RMSD), the closest non-bonded contact between units, and the formula.

Writes samples/triplet/dpa-wire.param.json and docs/triplet-report.md.
Run:  python tools/build_wire.py
"""

import json
import math
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolDescriptors

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_examples import check_bonds  # noqa: E402

RDLogger.DisableLog("rdApp.*")
sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "samples" / "triplet" / "dpa-wire.param.json"
REPORT = ROOT / "docs" / "triplet-report.md"
SEED = 20261105
C_H = 1.085  # aromatic C-H (Å), for the end caps

# DPA-2,6-diyl: [*:L] on anthracene C2, [*:R] on C6 (centrosymmetric partner). Atom order in this SMILES
# fixes the "local index" of each atom; the meso (9,10) phenyls carry R = H.
DPA = "[*:{L}]c1ccc2c(-c4ccccc4)c3cc([*:{R}])ccc3c(-c4ccccc4)c2c1"
LINKER = "[*:{L}]c1ccc([*:{R}])cc1"  # 1,4-phenylene
CAP = "[H][*:{M}]"


def frag(smiles, kind, k):
    m = Chem.MolFromSmiles(smiles, sanitize=False)
    for a in m.GetAtoms():
        if a.GetAtomicNum():
            a.SetProp("role", f"{kind}|{k}|{a.GetIdx()}")
    return m


def oligomer(n):
    """Wire with n DPA units as an RDKit mol (explicit H). Every atom gets a 'role' = kind|copy|local index."""
    parts = [frag(CAP.format(M=1), "capL", 0)]
    for k in range(n):
        parts.append(frag(DPA.format(L=2 * k + 1, R=2 * k + 2), "A", k))
        if k < n - 1:
            parts.append(frag(LINKER.format(L=2 * k + 2, R=2 * k + 3), "P", k))
    parts.append(frag(CAP.format(M=2 * n), "capR", 0))
    combo = parts[0]
    for p in parts[1:]:
        combo = Chem.CombineMols(combo, p)
    mol = Chem.molzip(combo)
    Chem.SanitizeMol(mol)
    mol = Chem.AddHs(mol)
    label_hydrogens(mol)
    return mol


def label_hydrogens(mol):
    """Added H atoms get their parent's role + "|H" (+ a number for the 2nd, 3rd H on the same atom, e.g. CH3)."""
    seen = {}
    for a in mol.GetAtoms():
        if not a.HasProp("role"):
            parent = a.GetNeighbors()[0].GetProp("role")
            j = seen[parent] = seen.get(parent, 0) + 1
            a.SetProp("role", parent + "|H" + (str(j) if j > 1 else ""))


def h_order(loc):
    """Stable atom order within a building block: heavy atoms by local index, then H atoms."""
    return ("|H" in loc, int(loc.split("|")[0]), loc)


def roles(mol):
    return [a.GetProp("role") for a in mol.GetAtoms()]


def local(role):
    """'A|3|7|H' -> ('A', 3, '7|H')"""
    kind, k, rest = role.split("|", 2)
    return kind, int(k), rest


def mmff(mol, conf_id=-1, iters=50000):
    props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
    ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=conf_id)
    ff.Minimize(maxIts=iters)
    return ff.CalcEnergy()


def embed_lowest(mol, n_confs=8):
    params = AllChem.ETKDGv3()
    params.randomSeed = SEED
    ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=params))
    res = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=20000, mmffVariant="MMFF94s")
    best = min(range(len(ids)), key=lambda i: res[i][1])
    return np.array(mol.GetConformer(ids[best]).GetPositions())


def extract_unit(mol, xyz, k):
    """DPA copy k + the linker after it, and the shift t from DPA k to DPA k+1."""
    unit = {"A": {}, "P": {}}
    lead = {}
    for i, r in enumerate(roles(mol)):
        kind, c, loc = local(r)
        if kind in ("A", "P") and c == k:
            unit[kind][loc] = xyz[i]
        if kind == "A" and c in (k, k + 1) and loc == "1":  # anthracene C2 (atom after the dummy)
            lead[c] = xyz[i]
    return unit, lead[k + 1] - lead[k]


def assemble(mol, unit, t, n):
    """Coordinates for oligomer(n) from the repeat unit - the same recipe as src/parametric.js."""
    caps = cap_positions(unit, t)
    xyz = np.zeros((mol.GetNumAtoms(), 3))
    for i, r in enumerate(roles(mol)):
        kind, k, loc = local(r)
        if kind in ("A", "P"):
            xyz[i] = unit[kind][loc] + k * t
        elif kind == "capL":
            xyz[i] = caps["L"]
        else:
            xyz[i] = caps["R"] + (n - 1) * t
    return xyz


def cap_positions(unit, t):
    left_c, right_c = unit["A"]["1"], unit["A"]["14"]  # C2 and C6 of the anthracene
    prev_far = unit["P"]["4"] - t  # para carbon of the previous unit's linker (bonded to C2)
    ipso = unit["P"]["1"]  # linker carbon bonded to C6
    unitv = lambda v: v / np.linalg.norm(v)  # noqa: E731
    return {"L": left_c + C_H * unitv(prev_far - left_c), "R": right_c + C_H * unitv(ipso - right_c)}


def set_coords(mol, xyz):
    mol.RemoveAllConformers()
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, p in enumerate(xyz):
        conf.SetAtomPosition(i, p.tolist())
    mol.AddConformer(conf, assignId=True)


def kabsch_rmsd(a, b):
    a0, b0 = a - a.mean(0), b - b.mean(0)
    u, _, vt = np.linalg.svd(a0.T @ b0)
    d = np.sign(np.linalg.det(u @ vt))
    r = u @ np.diag([1, 1, d]) @ vt
    return float(np.sqrt(((a0 @ r - b0) ** 2).sum(1).mean()))


def closest_contact(mol, xyz):
    """Shortest distance between atoms of different DPA/linker units that are not bonded (or 1-3 related)."""
    unit_of = []
    for r in roles(mol):
        kind, k, _ = local(r)
        unit_of.append(k if kind in ("A", "P") else -1)
    dm = Chem.GetDistanceMatrix(mol)
    best = (math.inf, None)
    for i in range(len(xyz)):
        for j in range(i + 1, len(xyz)):
            if unit_of[i] != unit_of[j] and dm[i, j] > 3:
                d = float(np.linalg.norm(xyz[i] - xyz[j]))
                if d < best[0]:
                    best = (d, (i, j))
    return best


def main():
    # 1. oligomer -> first guess of the repeat unit
    tri = oligomer(3)
    xyz = embed_lowest(tri)
    unit, t = extract_unit(tri, xyz, 1)
    print(f"trimer: {tri.GetNumAtoms()} atoms; first |t| = {np.linalg.norm(t):.3f} Å")

    # 2. self-consistency: assemble 5 units, relax, take the middle unit again
    penta = oligomer(5)
    for it in range(4):
        start = assemble(penta, unit, t, 5)
        set_coords(penta, start)
        mmff(penta)
        relaxed = penta.GetConformer().GetPositions()
        rmsd = kabsch_rmsd(start, relaxed)
        unit, t = extract_unit(penta, relaxed, 2)
        print(f"iteration {it + 1}: assembled vs relaxed pentamer RMSD {rmsd:.3f} Å, |t| = {np.linalg.norm(t):.3f} Å")

    # 3. validate assembled wires
    rows, worst = [], 0
    for n in (1, 2, 3, 6, 10):
        mol = oligomer(n)
        xyz = assemble(mol, unit, t, n)
        set_coords(mol, xyz)
        built = dict(mol=mol, coords=[tuple(p) for p in xyz], symbols=[a.GetSymbol() for a in mol.GetAtoms()],
                     extra_bonds=[], removed=set())
        per_type, outliers, _ = check_bonds(built)
        worst = max(worst, len(outliers))
        relax = mol.__copy__() if n <= 6 else None
        rmsd = None
        if relax is not None:
            mmff(relax)
            rmsd = kabsch_rmsd(xyz, relax.GetConformer().GetPositions())
        contact = closest_contact(mol, xyz) if n > 1 else (math.nan, None)
        caps = [i for i, r in enumerate(roles(mol)) if r.startswith("cap")]
        length = float(np.linalg.norm(xyz[caps[0]] - xyz[caps[1]]))
        formula = rdMolDescriptors.CalcMolFormula(mol)
        expect = f"C{26 * n + 6 * (n - 1)}H{18 * n + 2 * (n - 1)}"
        rows.append((n, mol.GetNumAtoms(), formula, expect, len(outliers), rmsd, contact[0], length, per_type))
        print(f"n={n}: {mol.GetNumAtoms()} atoms {formula} (expect {expect}), {len(outliers)} bonds outside range, "
              f"relax RMSD {rmsd if rmsd is None else round(rmsd, 3)}, closest contact {contact[0]:.2f} Å, "
              f"H···H length {length:.1f} Å")

    # geometry descriptors of the unit
    a, p = unit["A"], unit["P"]
    # local indices (see DPA / LINKER): A 1 = C2, 27 = C1, 14 = C6, 13 = C5, 19 = C9, 26 = C9a, 20/21 = meso
    # phenyl ipso/ortho; P 1 = ipso (bonded to C6), 2 = its ortho, 4 = para (bonded to the next C2), 3 = its ortho
    fold = lambda x: min(abs(x), 180 - abs(x))  # ring twist: 0-90°  # noqa: E731
    twist_left = fold(_dihedral(a["27"], a["1"], p["4"] - t, p["3"] - t))
    twist_right = fold(_dihedral(a["13"], a["14"], p["1"], p["2"]))
    meso_twist = fold(_dihedral(a["26"], a["19"], a["20"], a["21"]))

    data = {
        "builder": "dpa-wire",
        "description": "DPA-2,6-diyl + 1,4-phenylene repeat unit (R = H); see tools/build_wire.py",
        "method": "RDKit ETKDGv3 + MMFF94s oligomer, repeat unit made self-consistent in a relaxed pentamer",
        "t": [round(float(v), 4) for v in t],
        "dpa": [[s, *[round(float(v), 4) for v in xyz_]] for s, xyz_ in _atoms(unit["A"], tri, "A")],
        "linker": [[s, *[round(float(v), 4) for v in xyz_]] for s, xyz_ in _atoms(unit["P"], tri, "P")],
        "capLeft": [round(float(v), 4) for v in cap_positions(unit, t)["L"]],
        "capRight": [round(float(v), 4) for v in cap_positions(unit, t)["R"]],
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(data, separators=(",", ":")) + "\n", "utf-8")
    print("wrote", OUT.relative_to(ROOT), f"({len(data['dpa'])} + {len(data['linker'])} atoms per unit)")
    print(f"linker twist vs anthracene: {twist_left:.1f}° / {twist_right:.1f}°; meso phenyl twist {meso_twist:.1f}°")

    lines = ["# Triplet molecules: geometry report", "",
             "Generated by `tools/build_wire.py`.", "",
             "## Molecular wire (DPA-2,6-diyl / 1,4-phenylene, R = H)", "",
             f"- Repeat unit: {len(data['dpa'])} DPA atoms + {len(data['linker'])} linker atoms, "
             f"shift per unit |t| = {np.linalg.norm(t):.3f} Å",
             f"- Method: {data['method']}; the app stacks n copies and caps both ends with H (C–H {C_H} Å)",
             f"- Linker phenylene twisted {twist_left:.1f}° / {twist_right:.1f}° from the anthracene plane; "
             f"meso phenyls {meso_twist:.1f}°",
             "- *Relax RMSD*: how far MMFF94s moves the assembled wire when it is minimised freely "
             "(small = the stacked geometry is already at a minimum).",
             "- *Closest contact*: shortest distance between atoms of different units more than three bonds apart.",
             "",
             "| n | Atoms | Formula | Expected | Bonds outside range | Relax RMSD (Å) | Closest contact (Å) | H···H length (Å) |",
             "|---|---|---|---|---|---|---|---|"]
    for n, na, f, e, o, r, c, L, _ in rows:
        lines.append(f"| {n} | {na} | {f} | {e} | {o} | {'-' if r is None else f'{r:.3f}'} | "
                     f"{'-' if math.isnan(c) else f'{c:.2f}'} | {L:.1f} |")
    lines += ["", "Bond lengths in the n = 3 wire:", "", "| Bond | n | mean (Å) | min | max |", "|---|---|---|---|---|"]
    for key, ds in sorted(rows[2][8].items()):
        lines.append(f"| {key[0]}–{key[1]} {key[2].lower()} | {len(ds)} | {np.mean(ds):.3f} | {min(ds):.3f} | {max(ds):.3f} |")
    REPORT.write_text("\n".join(lines) + "\n", "utf-8")
    print("wrote", REPORT.relative_to(ROOT))
    return 1 if worst else 0


def _atoms(part, mol, kind):
    """(symbol, xyz) for a unit part, in a stable order (heavy atoms by local index, then their H)."""
    sym = {}
    for a in mol.GetAtoms():
        k, c, loc = local(a.GetProp("role"))
        if k == kind and c == 1:
            sym[loc] = a.GetSymbol()
    order = sorted(part, key=h_order)
    return [(sym[s], part[s]) for s in order]


def _dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - b0.dot(b1) * b1
    w = b2 - b2.dot(b1) * b1
    return math.degrees(math.atan2(np.cross(b1, v).dot(w), v.dot(w)))


if __name__ == "__main__":
    sys.exit(main())
