"""Build the parametric "DNA FRET ruler: Cy3–Cy5" object (Interesting photophysics tab).

    Cy3 | n base pairs of B-DNA, poly(dA)·poly(dT) | Cy5          (the dyes stacked on the two blunt ends)

DNA: one base pair (A5·T20 with its sugars and phosphates) from the Dickerson dodecamer (PDB 1BNA, 1.9 Å) and the
screw (twist about and rise along the helix axis) that carries it onto the next pair, A6·T19, fitted by Kabsch. The
app stacks n copies with that screw, so the O3'–P bonds between copies close up. No H on the DNA (an X-ray model,
like the proteins).

Dyes: Cy3 and Cy5 (non-sulfonated, N-methyl; the three-carbon tethers to the 5′ phosphates are not drawn), built by
RDKit + MMFF94s, stacked flat on the end base pairs, 3.4 Å beyond them. Each dye's turn about the helix axis and its
slide across the end pair (at most 2 Å) are scanned for the fit with the most room from the DNA; the dye is then
locked to its end pair. That is the picture from experiments on terminally attached Cy3/Cy5 (Iqbal et al., PNAS 105,
11176, 2008): the dyes stack on the ends and turn with the helix, so the angle between their transition dipoles
(along their long axes) changes with n, and FRET with it.

Writes samples/photophysics/dna-fret-ruler.param.json. Run:  python tools/build_fret.py
"""

import json
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_examples import CACHE, fetch  # noqa: E402

RDLogger.DisableLog("rdApp.*")
sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "samples" / "photophysics" / "dna-fret-ruler.param.json"
SEED = 20261111
STACK = 3.4  # Å between a dye's plane and the end base pair
R0 = 5.4  # nm, Cy3–Cy5 Förster radius for κ² = 2/3 (a typical value on DNA)
DYES = {
    "Cy3": "CN1c2ccccc2C(C)(C)/C1=C\\C=C\\C1=[N+](C)c2ccccc2C1(C)C",
    "Cy5": "CN1c2ccccc2C(C)(C)/C1=C\\C=C\\C=C\\C1=[N+](C)c2ccccc2C1(C)C",
}


def residues(text):
    res = {}
    for line in text.splitlines():
        if line.startswith("ATOM"):
            res.setdefault((line[21], int(line[22:26])), []).append(
                (line[12:16].strip(), line[76:78].strip() or line[12:16].strip()[0],
                 np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])))
    return res


def screw(p, q):
    """Rotation matrix m, translation t with q ~ p @ m.T + t (Kabsch), and the fit RMSD."""
    pc, qc = p.mean(0), q.mean(0)
    u, _, vt = np.linalg.svd((p - pc).T @ (q - qc))
    d = np.sign(np.linalg.det(u @ vt))
    r = u @ np.diag([1, 1, d]) @ vt  # p @ r ~ q
    m = r.T
    t = qc - pc @ r
    return m, t, float(np.sqrt((((p - pc) @ r + qc - q) ** 2).sum(1).mean()))


def dye(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    params = AllChem.ETKDGv3()
    params.randomSeed = SEED
    ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=10, params=params))
    res = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=20000, mmffVariant="MMFF94s")
    best = min(range(len(ids)), key=lambda i: res[i][1])
    xyz = np.array(mol.GetConformer(ids[best]).GetPositions())
    n = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "N"]
    conj = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in "CN" and a.GetHybridization() != Chem.HybridizationType.SP3]
    return mol, xyz, n, conj


def main():
    res = residues(fetch("https://files.rcsb.org/download/1BNA.pdb", CACHE / "pdb" / "1BNA.pdb"))
    unit = res[("A", 5)] + res[("B", 20)]
    nxt = {(c, name): p for c, k in (("A", 6), ("B", 19)) for name, _, p in res[(c, k)]}
    keys = [("A", name) for name, _, _ in res[("A", 5)]] + [("B", name) for name, _, _ in res[("B", 20)]]
    p = np.array([q for _, _, q in unit])
    q = np.array([nxt[k] for k in keys])
    m, t, rmsd = screw(p, q)

    # helix frame: z along the screw axis, origin on the axis at the base pair's height
    w, v = np.linalg.eig(m)
    axis = np.real(v[:, np.argmin(abs(w - 1))])
    rise = float(t @ axis)
    if rise < 0:
        axis, rise = -axis, -rise
    perp = t - rise * axis
    centre = np.linalg.lstsq(np.eye(3) - m, perp, rcond=None)[0]  # a point on the axis
    centre = centre + axis * float((p.mean(0) - centre) @ axis)
    x = np.cross(axis, [1, 0, 0])
    x /= np.linalg.norm(x)
    y = np.cross(axis, x)
    frame = np.array([x, y, axis])
    local = (p - centre) @ frame.T
    mz = frame @ m @ frame.T  # the screw's rotation in the frame: about z
    twist = float(np.arctan2(mz[1, 0], mz[0, 0]))
    print(f"step fit {rmsd:.2f} Å RMSD, twist {np.degrees(twist):.1f}°, rise {rise:.2f} Å")

    def place(xyz, k):
        c, s = np.cos(k * twist), np.sin(k * twist)
        rot = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
        return xyz @ rot.T + [0, 0, k * rise]

    # check: backbone bonds between copies
    names = [f"{c}:{name}" for c, name in keys]
    i_o3a, i_pa = names.index("A:O3'"), names.index("A:P")
    i_o3b, i_pb = names.index("B:O3'"), names.index("B:P")
    nxt_copy = place(local, 1)
    o3a = float(np.linalg.norm(local[i_o3a] - nxt_copy[i_pa]))
    o3b = float(np.linalg.norm(nxt_copy[i_o3b] - local[i_pb]))
    print(f"O3'–P between copies: strand A {o3a:.2f} Å, strand B {o3b:.2f} Å (a real one is 1.60)")

    # the base pair's long axis (C1'–C1') and centre, for stacking the dyes like one more pair
    c1a, c1b = local[names.index("A:C1'")], local[names.index("B:C1'")]
    long_axis = (c1a - c1b) * [1, 1, 0]
    long_axis /= np.linalg.norm(long_axis)
    base = [i for i, k in enumerate(keys) if k[1][-1:] not in ("'", "P") and k[1] not in ("OP1", "OP2", "O1P", "O2P")]
    bp_centre = local[base].mean(0)

    def rot_z(a):
        c, s_ = np.cos(a), np.sin(a)
        return np.array([[c, -s_, 0], [s_, c, 0], [0, 0, 1]])

    # Cy3 sits below pair 0 (meeting pairs 0-2), Cy5 above the last pair (meeting it and the two before it, which in
    # the last pair's frame are pairs -2..0). Both are stored in their end pair's frame: the app turns and lifts them
    # with it (Cy3 with pair 0, Cy5 with pair n - 1)
    ends = {"Cy3": (-STACK, np.vstack([place(local, k) for k in range(3)])),
            "Cy5": (STACK, np.vstack([place(local, k) for k in (-2, -1, 0)]))}
    dyes = {}
    for label, smi in DYES.items():
        height, end3 = ends[label]
        mol, xyz, n_idx, conj = dye(smi)
        normal = np.linalg.svd(xyz[conj] - xyz[conj].mean(0))[2][2]
        nn = xyz[n_idx[1]] - xyz[n_idx[0]]
        nn -= normal * (nn @ normal)
        nn /= np.linalg.norm(nn)
        src = np.array([nn, np.cross(normal, nn), normal])
        dst = np.array([long_axis, np.cross([0, 0, 1], long_axis), [0, 0, 1]])
        flat = (xyz - xyz[[i for i in conj] + n_idx].mean(0)) @ (dst.T @ src).T  # flat in xy, long axis along the pair's
        heavy = [i for i, a in enumerate(mol.GetAtoms()) if a.GetSymbol() != "H"]
        best = None
        for flip in (1, -1):  # either face down
            f = flat * [1, flip, flip]
            for deg in range(0, 360, 4):
                r = f @ rot_z(np.radians(deg)).T
                for dx in np.arange(-2, 2.01, 0.5):
                    for dy in np.arange(-2, 2.01, 0.5):
                        if dx * dx + dy * dy > 4:
                            continue
                        pos = r + [bp_centre[0] + dx, bp_centre[1] + dy, height]
                        room = float(np.min(np.linalg.norm(pos[heavy][:, None] - end3[None], axis=2)))
                        score = room - 0.05 * np.hypot(dx, dy)  # most room; slightly prefer staying centred
                        if best is None or score > best[0]:
                            best = (score, room, pos, deg, flip, dx, dy)
        _, room, pos, deg, flip, dx, dy = best
        dyes[label] = dict(atoms=[[a.GetSymbol(), *[round(float(c), 4) for c in pos[a.GetIdx()]]]
                                  for a in mol.GetAtoms()], dipole=n_idx)
        print(f"{label}: turned {deg}°, slid {np.hypot(dx, dy):.1f} Å, closest heavy atom to the DNA {room:.2f} Å")

    data = {
        "builder": "helix",
        "description": "Cy3 | n bp B-DNA | Cy5, the dyes stacked on the ends; see tools/build_fret.py",
        "method": f"DNA: the A5·T20 pair of PDB 1BNA repeated by its own helical step ({np.degrees(twist):.1f}°, "
                  f"{rise:.2f} Å; fit {rmsd:.2f} Å), no H. Dyes: RDKit + MMFF94s, stacked {STACK} Å beyond the end pairs",
        "unit": [[e, *[round(float(c), 4) for c in xyz]] for (_, e, _), xyz in zip(unit, local)],
        "twist": twist, "rise": rise, "stack": STACK, "r0": R0, "o3p": [round(o3a, 2), round(o3b, 2)],
        "cy3": dyes["Cy3"], "cy5": dyes["Cy5"],
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(data, separators=(",", ":")) + "\n", "utf-8")
    print("wrote", OUT.relative_to(ROOT))


if __name__ == "__main__":
    main()
