"""Build the end groups for the parametric "TTA Flora 4" object (Triplet molecules tab).

    sensitizer -(C6H4)- [DPA-2,6-diyl -(C6H4)-] x n - perylen-3-yl          (n = 0: panel a, one linker)

  sensitizer: thieno[2,3-b]quinoxaline 1,1-dioxide, 6,7-bis{4-[bis(4-methoxyphenyl)amino]phenyl}, linked at C2
  linkers:    1,4-phenylene
  mediators:  the DPA wire repeat unit (samples/triplet/dpa-wire.param.json, from tools/build_wire.py)
  emitter:    perylene, linked at C3

For n >= 1 the app builds the DPA wire exactly as for "Molecular wire" and adds two end groups:
  left  = sensitizer + its phenylene, in the frame of the first DPA;
  right = phenylene + perylene, in the frame of the last DPA (shifted by (n - 1) t in the app).
Each end group gets its own conformer search (RDKit ETKDGv3 + MMFF94s, end group + phenyl), its phenylene is
superimposed on the position a wire linker would take, and then the end group is relaxed with MMFF94s against
a fixed two-DPA wire. n = 0 is built directly (conformer search on the whole molecule).

Writes samples/triplet/tta-flora-4.param.json and appends to docs/triplet-report.md.
Run:  python tools/build_triplet.py   (after tools/build_wire.py)
"""

import json
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolDescriptors

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_examples import check_bonds  # noqa: E402
from build_wire import (CAP, DPA, LINKER, SEED, frag, h_order, kabsch_rmsd, label_hydrogens, local,  # noqa: E402
                        oligomer, roles, set_coords)

RDLogger.DisableLog("rdApp.*")
sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parent.parent
WIRE = ROOT / "samples" / "triplet" / "dpa-wire.param.json"
OUT = ROOT / "samples" / "triplet" / "tta-flora-4.param.json"
REPORT = ROOT / "docs" / "triplet-report.md"

ARYLAMINE = "c4ccc(N(c5ccc(OC)cc5)c5ccc(OC)cc5)cc4"  # 4-[bis(4-methoxyphenyl)amino]phenyl
SENS = "[*:{R}]C1=Cc2nc3cc(-" + ARYLAMINE + ")c(-" + ARYLAMINE + ")cc3nc2S1(=O)=O"
PERYLENE = "c1cc2cccc3c4cccc5cccc(c(c1)c23)c54"
RING = ["1", "2", "3", "4", "6", "7"]  # LINKER heavy atoms in ring order; 1 = ipso ([*:L]), 4 = para ([*:R])


def perylene_3yl():
    """Perylene with a dummy on C3 (peri: the CH whose quaternary neighbour also carries C4-H)."""
    m = Chem.MolFromSmiles(PERYLENE)
    assert rdMolDescriptors.CalcMolFormula(m) == "C20H12"
    for a in m.GetAtoms():
        quat = [n for n in a.GetNeighbors() if n.GetTotalNumHs() == 0]
        if a.GetTotalNumHs() == 1 and quat and sum(n.GetTotalNumHs() == 1 for n in quat[0].GetNeighbors()) == 2:
            rw = Chem.RWMol(m)
            d = rw.AddAtom(Chem.Atom(0))
            rw.GetAtomWithIdx(d).SetAtomMapNum(999)
            rw.AddBond(a.GetIdx(), d, Chem.BondType.SINGLE)
            Chem.SanitizeMol(rw)
            return Chem.MolToSmiles(rw).replace("[*:999]", "[*:{L}]")
    raise RuntimeError("no peri CH in perylene")


EMITTER = perylene_3yl()


def zipped(parts):
    combo = parts[0]
    for p in parts[1:]:
        combo = Chem.CombineMols(combo, p)
    mol = Chem.molzip(combo)
    Chem.SanitizeMol(mol)
    mol = Chem.AddHs(mol)
    label_hydrogens(mol)
    return mol


def flora(n):
    """The full molecule with n DPA mediators. Roles: S, PL, A|k, P|k, PR, E (and Q = the single linker if n = 0)."""
    if n == 0:
        return zipped([frag(SENS.format(R=1), "S", 0), frag(LINKER.format(L=1, R=2), "Q", 0),
                       frag(EMITTER.format(L=2), "E", 0)])
    parts = [frag(SENS.format(R=1), "S", 0), frag(LINKER.format(L=1, R=2), "PL", 0)]
    for k in range(n):
        parts.append(frag(DPA.format(L=2 * k + 2, R=2 * k + 3), "A", k))
        if k < n - 1:
            parts.append(frag(LINKER.format(L=2 * k + 3, R=2 * k + 4), "P", k))
    parts += [frag(LINKER.format(L=2 * n + 1, R=2 * n + 2), "PR", 0), frag(EMITTER.format(L=2 * n + 2), "E", 0)]
    return zipped(parts)


def search(mol, n_confs):
    params = AllChem.ETKDGv3()
    params.randomSeed = SEED
    params.pruneRmsThresh = 0.5
    ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=params))
    if len(ids) < n_confs // 2:  # ETKDG's default start fails for the n = 0 molecule; random starts work
        params.useRandomCoords = True
        ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=params))
    res = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=20000, mmffVariant="MMFF94s")
    best = min(range(len(ids)), key=lambda i: res[i][1])
    return np.array(mol.GetConformer(ids[best]).GetPositions()), len(ids)


def minimize(mol, fixed=()):
    props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
    ff = AllChem.MMFFGetMoleculeForceField(mol, props)
    for i in fixed:
        ff.AddFixedPoint(i)
    ff.Minimize(maxIts=100000)
    return np.array(mol.GetConformer().GetPositions())


def superpose(mobile, target):
    """Rotation r, translation such that mobile @ r + shift ~ target."""
    mc, tc = mobile.mean(0), target.mean(0)
    u, _, vt = np.linalg.svd((mobile - mc).T @ (target - tc))
    d = np.sign(np.linalg.det(u @ vt))
    r = u @ np.diag([1, 1, d]) @ vt
    return r, tc - mc @ r


def load_unit():
    """The wire repeat unit by role (same order as build_wire._atoms writes it)."""
    data = json.loads(WIRE.read_text("utf-8"))
    tri = oligomer(3)
    unit = {}
    for kind, key in (("A", "dpa"), ("P", "linker")):
        locs = [local(r)[2] for r in roles(tri) if local(r)[0] == kind and local(r)[1] == 1]
        order = sorted(locs, key=h_order)
        assert len(order) == len(data[key])
        unit[kind] = {loc: np.array(row[1:]) for loc, row in zip(order, data[key])}
    return data, unit, np.array(data["t"])


def end_group(parts, side, target_ring, n_confs):
    """Conformer search on an end group + phenyl, then place its phenylene on target_ring (dict loc -> xyz)."""
    mol = zipped(parts)
    xyz, n_emb = search(mol, n_confs)
    ring_idx = {local(r)[2]: i for i, r in enumerate(roles(mol)) if local(r)[0] == side}
    best = None
    for order in (RING, ["1", "7", "6", "4", "3", "2"]):  # both ways round the ring (flips the end group)
        mob = np.array([xyz[ring_idx[a]] for a in order])
        tgt = np.array([target_ring[a] for a in RING])
        r, shift = superpose(mob, tgt)
        placed = xyz @ r + shift
        rms = float(np.sqrt(((mob @ r + shift - tgt) ** 2).sum(1).mean()))
        if best is None or rms < best[0]:
            best = (rms, placed, order)
    return mol, best[1], best[0], n_emb


def assemble(mol, unit, t, n, ends):
    xyz = np.zeros((mol.GetNumAtoms(), 3))
    for i, r in enumerate(roles(mol)):
        kind, k, loc = local(r)
        if kind in ("A", "P"):
            xyz[i] = unit[kind][loc] + k * t
        elif kind in ("S", "PL"):
            xyz[i] = ends[(kind, loc)]
        else:  # PR, E
            xyz[i] = ends[(kind, loc)] + (n - 1) * t
    return xyz


def rows_of(mol, xyz, kinds):
    """[[symbol, x, y, z], ...] for atoms of the given kinds, plus their role keys (for core indices)."""
    out, keys = [], []
    for i, (a, r) in enumerate(zip(mol.GetAtoms(), roles(mol))):
        kind, _, loc = local(r)
        if kind in kinds:
            out.append([a.GetSymbol(), *[round(float(v), 4) for v in xyz[i]]])
            keys.append((kind, loc))
    return out, keys


SENS_CORE = "c1nc2ccccc2nc1"  # quinoxaline of the sensitizer


def expected_formula(n):
    # hand count: sensitizer-H C50H40N4O6S; panel a = C76H54N4O6S; each mediator adds DPA-diyl + C6H4 = C32H20
    return f"C{76 + 32 * n}H{54 + 20 * n}N4O6S"


def main():
    wire, unit, t = load_unit()
    print(f"wire unit loaded, |t| = {np.linalg.norm(t):.3f} Å; emitter fragment {EMITTER}")

    # 1. end groups: conformer search, then superpose their phenylene on the wire linker positions
    left_mol, left_xyz, left_rms, nl = end_group(
        [frag(SENS.format(R=1), "S", 0), frag(LINKER.format(L=1, R=2), "PL", 0), frag(CAP.format(M=2), "cap", 0)],
        "PL", {a: unit["P"][a] - t for a in RING}, 30)
    right_mol, right_xyz, right_rms, nr = end_group(
        [frag(CAP.format(M=1), "cap", 0), frag(LINKER.format(L=1, R=2), "PR", 0), frag(EMITTER.format(L=2), "E", 0)],
        "PR", {a: unit["P"][a] for a in RING}, 10)
    print(f"left end: {nl} conformers, phenylene fit {left_rms:.3f} Å; right end: {nr} conformers, fit {right_rms:.3f} Å")

    ends = {}
    for mol, xyz in ((left_mol, left_xyz), (right_mol, right_xyz)):
        for i, r in enumerate(roles(mol)):
            kind, _, loc = local(r)
            if kind != "cap":
                ends[(kind, loc)] = xyz[i]

    # 2. relax both end groups against a fixed two-DPA wire, keep them as the end-group geometry
    m2 = flora(2)
    xyz = assemble(m2, unit, t, 2, ends)
    set_coords(m2, xyz)
    wire_atoms = [i for i, r in enumerate(roles(m2)) if local(r)[0] in ("A", "P")]
    relaxed = minimize(m2, wire_atoms)
    moved = kabsch_rmsd(xyz, relaxed)
    for i, r in enumerate(roles(m2)):
        kind, _, loc = local(r)
        if kind in ("S", "PL"):
            ends[(kind, loc)] = relaxed[i]
        elif kind in ("PR", "E"):
            ends[(kind, loc)] = relaxed[i] - t
    print(f"end groups relaxed against the fixed wire (moved {moved:.3f} Å RMSD)")

    # 3. n = 0 (panel a): its own conformer search
    m0 = flora(0)
    zero_xyz, n0 = search(m0, 30)
    print(f"n = 0: {m0.GetNumAtoms()} atoms, {n0} conformers")

    # 4. validate
    report = []
    worst = 0
    for n in (0, 1, 2, 3, 5, 10):
        mol = flora(n)
        xyz = zero_xyz if n == 0 else assemble(mol, unit, t, n, ends)
        set_coords(mol, xyz)
        built = dict(mol=mol, coords=[tuple(p) for p in xyz], symbols=[a.GetSymbol() for a in mol.GetAtoms()],
                     extra_bonds=[], removed=set())
        per_type, outliers, _ = check_bonds(built)
        worst = max(worst, len(outliers))
        formula = rdMolDescriptors.CalcMolFormula(mol)
        # closest contact between different building blocks, > 3 bonds apart
        block = [(local(r)[0], local(r)[1]) for r in roles(mol)]
        dm = Chem.GetDistanceMatrix(mol)
        dd = np.linalg.norm(xyz[:, None] - xyz[None], axis=2)
        mask = (dm > 3) & np.array([[a != b for b in block] for a in block])
        contact = float(dd[mask].min())
        rmsd = None
        if n <= 2:
            free = Chem.Mol(mol)
            minimize(free)
            rmsd = kabsch_rmsd(xyz, np.array(free.GetConformer().GetPositions()))
        s_idx = [i for i in mol.GetSubstructMatch(Chem.MolFromSmarts(SENS_CORE))]
        e_idx = [i for i, r in enumerate(roles(mol)) if local(r)[0] == "E" and mol.GetAtomWithIdx(i).GetSymbol() == "C"]
        dist = float(np.linalg.norm(xyz[s_idx].mean(0) - xyz[e_idx].mean(0)))
        report.append((n, mol.GetNumAtoms(), formula, expected_formula(n), len(outliers), rmsd, contact, dist, outliers))
        print(f"n={n}: {mol.GetNumAtoms()} atoms {formula} (expect {expected_formula(n)}), "
              f"{len(outliers)} bonds outside range, relax RMSD {'-' if rmsd is None else f'{rmsd:.3f}'} Å, "
              f"closest contact {contact:.2f} Å, sensitizer-perylene {dist / 10:.2f} nm")
        for o in outliers:
            print("   outlier", o)

    # 5. write data for the app: left/right end groups, the n = 0 molecule, and core indices for the distance
    lm = flora(1)
    xyz1 = assemble(lm, unit, t, 1, ends)
    set_coords(lm, xyz1)
    left, lkeys = rows_of(lm, xyz1, ("S", "PL"))
    right, rkeys = rows_of(lm, xyz1, ("PR", "E"))
    s_core = {local(lm.GetAtomWithIdx(i).GetProp("role"))[2] for i in lm.GetSubstructMatch(Chem.MolFromSmarts(SENS_CORE))}
    zero, zkeys = rows_of(m0, zero_xyz, ("S", "Q", "E"))
    s_core0 = {local(m0.GetAtomWithIdx(i).GetProp("role"))[2] for i in m0.GetSubstructMatch(Chem.MolFromSmarts(SENS_CORE))}
    data = {
        "builder": "tta-flora",
        "wire": WIRE.name,
        "description": "sensitizer-(C6H4)-[DPA-(C6H4)]n-perylen-3-yl; see tools/build_triplet.py",
        "method": "RDKit ETKDGv3 + MMFF94s: end groups by conformer search, relaxed against the stacked DPA wire",
        "left": left,
        "right": right,
        "zero": zero,
        "leftCore": [j for j, (k, loc) in enumerate(lkeys) if k == "S" and loc in s_core],
        "rightCore": [j for j, (k, loc) in enumerate(rkeys) if k == "E" and right[j][0] == "C"],
        "zeroLeftCore": [j for j, (k, loc) in enumerate(zkeys) if k == "S" and loc in s_core0],
        "zeroRightCore": [j for j, (k, loc) in enumerate(zkeys) if k == "E" and zero[j][0] == "C"],
    }
    assert len(data["leftCore"]) == len(data["zeroLeftCore"]) == 10 and len(data["rightCore"]) == 20
    OUT.write_text(json.dumps(data, separators=(",", ":")) + "\n", "utf-8")
    print("wrote", OUT.relative_to(ROOT), f"(left {len(left)}, right {len(right)}, n=0 {len(zero)} atoms)")

    text = REPORT.read_text("utf-8").split("\n## TTA Flora 4")[0].rstrip("\n")
    lines = ["", "## TTA Flora 4 (sensitizer – phenylene – [DPA – phenylene]n – perylene)", "",
             "- Sensitizer: 6,7-bis{4-[bis(4-methoxyphenyl)amino]phenyl}thieno[2,3-b]quinoxaline 1,1-dioxide, "
             "linked at thiophene C2. Emitter: perylene, linked at C3. Read from the project figure.",
             f"- End groups: conformer search ({nl} / {nr} conformers, end group + phenyl), phenylene superimposed on "
             f"the wire linker position (fit {left_rms:.3f} / {right_rms:.3f} Å), then relaxed against a fixed "
             f"two-DPA wire (moved {moved:.3f} Å RMSD). n = 0 (panel a): own conformer search ({n0} conformers).",
             "- Expected formula: hand count, C(76 + 32n) H(54 + 20n) N4 O6 S.",
             "- *Relax RMSD*: how far a free MMFF94s minimisation moves the assembled molecule (the flexible "
             "methoxyphenyl arms contribute most). *Distance*: quinoxaline centroid to perylene centroid.", "",
             "| n | Atoms | Formula | Expected | Bonds outside range | Relax RMSD (Å) | Closest contact (Å) | Sensitizer–perylene (nm) |",
             "|---|---|---|---|---|---|---|---|"]
    for n, na, f, e, o, r, c, d, _ in report:
        lines.append(f"| {n} | {na} | {f} | {e} | {o} | {'-' if r is None else f'{r:.3f}'} | {c:.2f} | {d / 10:.2f} |")
    REPORT.write_text(text + "\n" + "\n".join(lines) + "\n", "utf-8")
    print("wrote", REPORT.relative_to(ROOT))
    return 1 if worst else 0


if __name__ == "__main__":
    sys.exit(main())
