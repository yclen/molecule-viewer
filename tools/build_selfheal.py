"""Build the parametric "Cy5–(phenylene)n–COT" object (Self-healing dyes tab): our rigid-linker self-healing dye.

    Cy5 (indolenine C5) -(C6H4)- x n - cycloocta-1,3,5,7-tetraen-1-yl          (n = 0: COT bonded straight to C5)

The published self-healing dyes join the triplet quencher (cyclooctatetraene, COT) to the dye through a flexible
alkyl/amide chain on the indolenine nitrogen, so triplet energy transfer (TTET) needs the two to touch. Here the
quencher sits on the dye's benzo ring (where sulfo-Cy5 carries a sulfonate) at the end of a rigid 1,4-phenylene
rod, the linker used in TTA Flora 4, so the coupling goes through bonds and the distance is fixed. The other
indolenine keeps its sulfonate and the N-hexanoic acid handle, as in the app's Cy5.

Geometry (RDKit ETKDGv3 + MMFF94s):
  1. the dye alone (H at C5) gets a conformer search; its lowest conformer is the dye pose for every n;
  2. for each n, the whole molecule is embedded with the dye's heavy atoms held at that pose, minimised with
     them fixed, and the lowest of the conformers kept (this searches the biaryl twists and the COT tub).
So the dye does not jump around when n changes: only the rod and the COT move.

Writes samples/selfheal/cy5-phenylene-cot.param.json and docs/selfheal-report.md.
Run:  python tools/build_selfheal.py   (needs rdkit)
"""

import json
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Geometry import Point3D

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_examples import check_bonds  # noqa: E402
from build_wire import LINKER, SEED, frag, label_hydrogens, local, roles  # noqa: E402

RDLogger.DisableLog("rdApp.*")
sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "samples" / "selfheal" / "cy5-phenylene-cot.param.json"
REPORT = ROOT / "docs" / "selfheal-report.md"
N_MAX = 4

# The app's sulfo-Cy5 (PubChem 17758493) with the left indolenine's 5-sulfonate replaced by the attachment point.
# One sulfonate is left, so the molecule is a neutral zwitterion.
DYE = ("CCN\\1C2=C(C=C(C=C2)[*:{R}])C(/C1=C/C=C/C=C/C3=[N+](C4=C(C3(C)C)C=C(C=C4)S(=O)(=O)[O-])CCCCCC(=O)O)(C)C")
COT = "[*:{L}]C1=CC=CC=CC=C1"  # cycloocta-1,3,5,7-tetraen-1-yl
CAP = "[H][*:{M}]"
POLYMETHINE = "[#7]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#7]"  # N1...N1': the Cy5 chromophore's spine


def zipped(parts):
    combo = parts[0]
    for p in parts[1:]:
        combo = Chem.CombineMols(combo, p)
    mol = Chem.molzip(combo)
    Chem.SanitizeMol(mol)
    mol = Chem.AddHs(mol)
    label_hydrogens(mol)
    return mol


def conjugate(n):
    """The molecule with n phenylenes. Roles: D (dye), P|k (phenylene k), Q (COT)."""
    if n == 0:
        return zipped([frag(DYE.format(R=1), "D", 0), frag(COT.format(L=1), "Q", 0)])
    parts = [frag(DYE.format(R=1), "D", 0)]
    for k in range(n):
        parts.append(frag(LINKER.format(L=k + 1, R=k + 2), "P", k))
    parts.append(frag(COT.format(L=n + 1), "Q", 0))
    return zipped(parts)


def mmff_confs(mol, n_confs, coord_map=None, fixed=()):
    params = AllChem.ETKDGv3()
    params.randomSeed = SEED
    params.pruneRmsThresh = 0.3
    if coord_map:
        params.SetCoordMap(coord_map)
    ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=params))
    if not ids:
        params.useRandomCoords = True
        ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=params))
    props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
    best = None
    for cid in ids:
        if coord_map:  # put the held atoms exactly back where they belong, then relax the rest around them
            xyz = np.array(mol.GetConformer(cid).GetPositions())
            keys = sorted(coord_map)
            r, shift = superpose(xyz[keys], np.array([[coord_map[k].x, coord_map[k].y, coord_map[k].z] for k in keys]))
            conf = mol.GetConformer(cid)
            for i, p in enumerate(xyz @ r + shift):
                conf.SetAtomPosition(i, Point3D(*p))
        ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=cid)
        for i in fixed:
            ff.AddFixedPoint(i)
        ff.Minimize(maxIts=50000)
        e = ff.CalcEnergy()
        if best is None or e < best[0]:
            best = (e, cid)
    return np.array(mol.GetConformer(best[1]).GetPositions()), len(ids)


def superpose(mobile, target):
    """Rotation r, translation such that mobile @ r + shift ~ target."""
    mc, tc = mobile.mean(0), target.mean(0)
    u, _, vt = np.linalg.svd((mobile - mc).T @ (target - tc))
    d = np.sign(np.linalg.det(u @ vt))
    r = u @ np.diag([1, 1, d]) @ vt
    return r, tc - mc @ r


def plane_normal(pts):
    pts = np.asarray(pts)
    return np.linalg.svd(pts - pts.mean(0))[2][2]


def twist(a, b):
    """Angle between two ring planes, 0-90°."""
    c = abs(float(np.dot(plane_normal(a), plane_normal(b))))
    return float(np.degrees(np.arccos(min(1.0, c))))


def main():
    # 1. the dye pose: conformer search on the dye with H at C5
    dye_mol = zipped([frag(DYE.format(R=1), "D", 0), frag(CAP.format(M=1), "cap", 0)])
    dye_xyz, n_dye = mmff_confs(dye_mol, 30)
    dye_pos = {r: dye_xyz[i] for i, r in enumerate(roles(dye_mol)) if local(r)[0] == "D" and "|H" not in r}
    print(f"dye: {dye_mol.GetNumAtoms()} atoms, lowest of {n_dye} conformers")

    frames, report = {}, []
    for n in range(N_MAX + 1):
        mol = conjugate(n)
        rl = roles(mol)
        held = {i: Point3D(*dye_pos[r]) for i, r in enumerate(rl) if r in dye_pos}
        xyz, n_conf = mmff_confs(mol, 20 if n else 30, held, fixed=list(held))

        # atom order for the app: dye, phenylenes in order, COT (each block followed by nothing else)
        def block(kind, k=None):
            return [i for i, r in enumerate(rl) if local(r)[0] == kind and (k is None or local(r)[1] == k)]
        order = block("D") + [i for k in range(n) for i in block("P", k)] + block("Q")
        assert sorted(order) == list(range(mol.GetNumAtoms()))
        atoms = [[mol.GetAtomWithIdx(i).GetSymbol(), *[round(float(v), 4) for v in xyz[i]]] for i in order]
        pos = {old: new for new, old in enumerate(order)}
        dye_sites = [pos[i] for i in block("D")]
        cot_sites = [pos[i] for i in block("Q")]

        # checks
        conf_mol = Chem.Mol(mol)
        conf_mol.RemoveAllConformers()
        conf = Chem.Conformer(mol.GetNumAtoms())
        for i, p in enumerate(xyz):
            conf.SetAtomPosition(i, Point3D(*p))
        conf_mol.AddConformer(conf)
        built = dict(mol=conf_mol, coords=[tuple(p) for p in xyz], symbols=[a.GetSymbol() for a in mol.GetAtoms()],
                     extra_bonds=[], removed=set())
        _, outliers, _ = check_bonds(built)
        formula = rdMolDescriptors.CalcMolFormula(mol)
        dm = Chem.GetDistanceMatrix(mol)
        dd = np.linalg.norm(xyz[:, None] - xyz[None], axis=2)
        blk = [(local(r)[0], local(r)[1]) for r in rl]
        mask = (dm > 3) & np.array([[a != b for b in blk] for a in blk])
        contact = float(dd[mask].min()) if mask.any() else float("nan")

        spine = mol.GetSubstructMatch(Chem.MolFromSmarts(POLYMETHINE))
        cot_c = [i for i in block("Q") if mol.GetAtomWithIdx(i).GetSymbol() == "C"]
        dist = float(np.linalg.norm(xyz[list(spine)].mean(0) - xyz[cot_c].mean(0)))
        # the bond the quencher hangs on, and the planes either side of every aryl-aryl / aryl-COT bond
        attach = next(i for i in block("D") if any(local(rl[nb.GetIdx()])[0] != "D"
                                                  for nb in mol.GetAtomWithIdx(i).GetNeighbors()))
        benzo = next(list(r) for r in mol.GetRingInfo().AtomRings() if len(r) == 6 and attach in r)
        rings = [benzo] + [[i for i in block("P", k) if mol.GetAtomWithIdx(i).GetSymbol() == "C"] for k in range(n)]
        q_ipso = next(i for i in cot_c if any(local(rl[nb.GetIdx()])[0] != "Q" for nb in mol.GetAtomWithIdx(i).GetNeighbors()))
        q_plane = [q_ipso] + [nb.GetIdx() for nb in mol.GetAtomWithIdx(q_ipso).GetNeighbors()]
        twists = [twist(xyz[a], xyz[b]) for a, b in zip(rings, rings[1:])] + [twist(xyz[rings[-1]], xyz[q_plane])]
        edge = float(np.linalg.norm(xyz[attach] - xyz[q_ipso]))  # along the rod: dye C5 to the COT carbon on it
        report.append((n, mol.GetNumAtoms(), formula, len(outliers), contact, dist, edge, twists, n_conf))
        tw = " / ".join(f"{t:.0f}°" for t in twists)
        print(f"n={n}: {mol.GetNumAtoms()} atoms {formula}, {len(outliers)} bonds outside range, closest contact "
              f"{contact:.2f} Å, Cy5 spine to COT {dist / 10:.2f} nm, C5 to COT {edge:.1f} Å, twists {tw}")
        for o in outliers:
            print("   outlier", o)

        what = "COT bonded straight to the dye" if n == 0 else f"{n} phenylene{'s' if n > 1 else ''}"
        frames[str(n)] = {
            "atoms": atoms,
            "sites": [{"label": "Cy5", "kind": "dye", "atoms": dye_sites},
                      {"label": "COT", "kind": "quencher", "atoms": cot_sites}],
            "comment": f"Cy5–(C6H4){n}–COT, {what} | {dist / 10:.2f} nm Cy5 spine to COT | twists {tw}",
            "note": f"Cy5 ↔ COT: {dist / 10:.2f} nm centre to centre, {edge:.1f} Å from the dye's C5 to the COT. "
                    f"Ring twists along the rod: {tw}",
        }

    data = {"builder": "frames",
            "description": "Cy5 (indolenine C5) -(C6H4)n- COT; see tools/build_selfheal.py",
            "method": "RDKit ETKDGv3 + MMFF94s, dye pose shared by every n",
            "frames": frames}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(data, separators=(",", ":")) + "\n", "utf-8")
    print("wrote", OUT.relative_to(ROOT))

    lines = ["# Self-healing dyes: geometry report", "",
             "Generated by `tools/build_selfheal.py`.", "",
             "## Cy5–(phenylene)n–COT", "",
             "- Cy5 = the app's sulfo-Cy5 with the left indolenine's 5-sulfonate replaced by the rod. COT = "
             "cycloocta-1,3,5,7-tetraen-1-yl.",
             f"- Dye pose: lowest of {n_dye} MMFF94s conformers of the dye alone, held fixed for every n; the rod and "
             "the COT are searched around it.",
             "- *Cy5–COT*: centroid of the Cy5 spine (N–C7–N) to the COT ring centroid; *C5–COT*: the rod, from the dye's "
             "C5 to the COT carbon bonded to it. *Twists*: angle between neighbouring ring planes along the rod, from the dye's benzo "
             "ring to the COT's ipso C=C (an angle near 90° cuts π overlap).",
             "- *Closest contact*: shortest distance between atoms of different blocks more than three bonds apart.", "",
             "| n | Atoms | Formula | Bonds outside range | Closest contact (Å) | Cy5–COT (nm) | C5–COT (Å) | Twists | Conformers |",
             "|---|---|---|---|---|---|---|---|---|"]
    for n, na, f, o, c, d, e, tws, nc in report:
        lines.append(f"| {n} | {na} | {f} | {o} | {c:.2f} | {d / 10:.2f} | {e:.1f} | "
                     f"{' / '.join(f'{t:.0f}°' for t in tws)} | {nc} |")
    REPORT.write_text("\n".join(lines) + "\n", "utf-8")
    print("wrote", REPORT.relative_to(ROOT))
    return 1 if any(r[3] for r in report) else 0


if __name__ == "__main__":
    sys.exit(main())
