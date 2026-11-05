"""Build the example molecules in tools/catalog.py into samples/<tab>/*.xyz, then check them.

Small molecules:  SMILES -> RDKit ETKDGv3 conformer search -> MMFF94s minimisation (UFF if MMFF lacks
                  parameters, e.g. boron) -> lowest-energy conformer.
Proteins:         X-ray structure from the RCSB PDB (chain A, no waters, first alternate location).

Validation (written to docs/geometry-report.md):
  1. every bond length is compared against typical crystallographic ranges for that element pair and
     bond order (Allen et al. 1987, J. Chem. Soc. Perkin Trans. 2, S1-S19);
  2. where PubChem has a 3D conformer for the same skeleton, bond lengths are compared bond-for-bond.

Run:  python tools/build_examples.py            (network needed the first time; results are cached)
      python tools/build_examples.py --only cy5  (rebuild one molecule)
"""

import argparse
import json
import math
import statistics
import sys
import urllib.request
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

sys.path.insert(0, str(Path(__file__).resolve().parent))
from catalog import MOLECULES, TABS  # noqa: E402

RDLogger.DisableLog("rdApp.*")
sys.stdout.reconfigure(encoding="utf-8")  # Windows consoles/pipes default to a code page without Δ, Å

ROOT = Path(__file__).resolve().parent.parent
SAMPLES = ROOT / "samples"
CACHE = Path(__file__).resolve().parent / "cache"
REPORT = ROOT / "docs" / "geometry-report.md"
SEED = 20261105

# Typical bond-length windows (Å) by (element pair, bond type). Upper C-C / C-N limits allow for
# quaternary and N+-alkyl carbons (~1.55-1.57 / ~1.51 Å), which are common in these dyes.
RANGES = {
    ("C", "C", "SINGLE"): (1.41, 1.58), ("C", "C", "DOUBLE"): (1.31, 1.44),
    ("C", "C", "AROMATIC"): (1.35, 1.46), ("C", "C", "TRIPLE"): (1.18, 1.23),
    ("C", "H", "SINGLE"): (1.07, 1.11),
    ("C", "N", "SINGLE"): (1.32, 1.52), ("C", "N", "DOUBLE"): (1.27, 1.40),
    ("C", "N", "AROMATIC"): (1.30, 1.42), ("C", "N", "TRIPLE"): (1.13, 1.18),
    ("C", "O", "SINGLE"): (1.24, 1.46), ("C", "O", "DOUBLE"): (1.19, 1.28), ("C", "O", "AROMATIC"): (1.33, 1.40),
    ("C", "S", "SINGLE"): (1.74, 1.84), ("O", "S", "SINGLE"): (1.42, 1.62), ("O", "S", "DOUBLE"): (1.40, 1.48),
    ("Cl", "S", "SINGLE"): (1.98, 2.08),
    ("C", "F", "SINGLE"): (1.32, 1.40), ("C", "Cl", "SINGLE"): (1.70, 1.77),
    ("C", "Br", "SINGLE"): (1.85, 1.93), ("C", "I", "SINGLE"): (2.04, 2.13),
    ("B", "F", "SINGLE"): (1.35, 1.43), ("B", "N", "SINGLE"): (1.52, 1.62),
    ("H", "N", "SINGLE"): (0.99, 1.04), ("H", "O", "SINGLE"): (0.95, 1.01),
    ("N", "Pt", "SINGLE"): (1.98, 2.08),  # Pt(II) porphyrins
    ("B", "C", "SINGLE"): (1.48, 1.60), ("B", "O", "SINGLE"): (1.44, 1.54),
    ("C", "Se", "SINGLE"): (1.88, 1.96), ("C", "Si", "SINGLE"): (1.80, 1.92),
    ("N", "S", "AROMATIC"): (1.58, 1.66), ("N", "N", "AROMATIC"): (1.33, 1.42),
    ("C", "O", "AROMATIC"): (1.33, 1.40),
}


# --------------------------------------------------------------------------- helpers

def fetch(url, path):
    """Download once, then reuse the cached copy."""
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        req = urllib.request.Request(url, headers={"User-Agent": "molecule-viewer-build/0.1"})
        with urllib.request.urlopen(req, timeout=60) as r:
            path.write_bytes(r.read())
    return path.read_text(encoding="utf-8", errors="replace")


def xyz_text(symbols, coords, comment):
    lines = [str(len(symbols)), comment]
    lines += [f"{s:<2} {x:12.6f} {y:12.6f} {z:12.6f}" for s, (x, y, z) in zip(symbols, coords)]
    return "\n".join(lines) + "\n"


def bond_key(a, b, btype):
    e1, e2 = sorted((a, b))
    return (e1, e2, btype)


# --------------------------------------------------------------------------- small molecules

# Known force-field errors, corrected wherever the pattern occurs (first and last SMARTS atoms are restrained).
GLOBAL_RESTRAINTS = [
    # MMFF94 puts Si–C ~0.1 Å too long; X-ray TIPS-ethynyl compounds: Si–C(sp) 1.84, Si–C(iPr) 1.88 Å
    ("smarts", "[Si]-[CX2]", 1.84, "Si–C(sp)"),
    ("smarts", "[Si]-[CX4]", 1.88, "Si–C(sp3)"),
    # MMFF types 2,1,3-benzothiadiazole as quinoid and stretches the ring-fusion C–C to ~1.51 Å (X-ray ~1.43)
    ("smarts", "c1nsnc1", 1.43, "benzothiadiazole C3a–C7a"),
]
# UFF makes sp3 C–H ~1.11 Å; neutron-diffraction mean ~1.09 Å
UFF_RESTRAINTS = [("smarts", "[CX4]-[#1]", 1.09, "C(sp3)–H")]

# Bump when the build method changes, to invalidate cached geometries.
BUILD_VERSION = 3


def build_small_cached(entry):
    """build_small() with an on-disk cache keyed on the catalog entry, so unchanged molecules are reused."""
    import hashlib
    import pickle
    key = hashlib.sha1(json.dumps([entry, BUILD_VERSION, GLOBAL_RESTRAINTS, UFF_RESTRAINTS],
                                  sort_keys=True, default=str).encode()).hexdigest()[:16]
    path = CACHE / "built" / f"{key}.pkl"
    if path.exists():
        data = pickle.loads(path.read_bytes())
        data["mol"] = Chem.Mol(data["mol"])
        data["cached"] = True
        return data
    built = build_small(entry)
    data = {k: v for k, v in built.items() if k != "conf"}
    data["mol"] = built["mol"].ToBinary()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(pickle.dumps(data))
    return built


def build_small(entry):
    mol = Chem.MolFromSmiles(entry["smiles"])
    if mol is None:
        raise ValueError(f"bad SMILES for {entry['name']}")
    mol = rdMolStandardize.LargestFragmentChooser().choose(mol)  # drop counter-ions
    mol = Chem.AddHs(mol)

    heavy = mol.GetNumHeavyAtoms()
    n_rot = rdMolDescriptors.CalcNumRotatableBonds(mol)
    n_confs = entry.get("confs") or (30 if heavy < 40 and n_rot < 10 else 20)

    params = AllChem.ETKDGv3()
    params.randomSeed = SEED
    params.pruneRmsThresh = 0.25
    conf_ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=params))
    if not conf_ids:
        params.useRandomCoords = True
        conf_ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=params))
    if not conf_ids:
        raise RuntimeError(f"could not embed {entry['name']}")

    if AllChem.MMFFHasAllMoleculeParams(mol):
        ff = "MMFF94s"
        results = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=20000, mmffVariant="MMFF94s")
    else:
        ff = "UFF"
        results = AllChem.UFFOptimizeMoleculeConfs(mol, maxIters=20000)

    # lowest energy among converged conformers (fall back to all if none converged)
    ok = [i for i, (nc, _) in enumerate(results) if nc == 0] or list(range(len(results)))
    best = min(ok, key=lambda i: results[i][1])
    conf_id = conf_ids[best]
    method = f"RDKit ETKDGv3 + {ff}, lowest of {len(conf_ids)} conformers"

    # Bonds the force field gets wrong are pulled to crystallographic means.
    rules = list(entry.get("restraints", [])) + GLOBAL_RESTRAINTS + (UFF_RESTRAINTS if ff == "UFF" else [])
    restraints, used = [], set()
    for r in rules:
        if r[0] == "smarts":
            _, pattern, d, label = r
            for match in mol.GetSubstructMatches(Chem.MolFromSmarts(pattern), uniquify=True, maxMatches=100000):
                restraints.append((match[0], match[-1], d))
                used.add(f"{label} {d:.3f} Å")
        else:
            for b in mol.GetBonds():
                e = sorted((b.GetBeginAtom().GetSymbol(), b.GetEndAtom().GetSymbol()))
                if sorted(r[:2]) == e:
                    restraints.append((b.GetBeginAtomIdx(), b.GetEndAtomIdx(), r[2]))
                    used.add(f"{r[0]}–{r[1]} {r[2]:.3f} Å")
    if restraints:
        field = _force_field(mol, conf_id, ff)
        for i, j, d in restraints:
            add = field.MMFFAddDistanceConstraint if ff.startswith("MMFF") else field.UFFAddDistanceConstraint
            add(i, j, False, d, d, 5000.0)
        field.Minimize(maxIts=20000)
        method += f"; restrained {', '.join(sorted(used))} (crystallographic means)"

    extra_bonds = []  # bonds added by post-processing (e.g. metal–N), for validation
    post = entry.get("post", "")
    if post.startswith("porphyrin:"):
        metal = post.split(":", 1)[1]
        extra_bonds = porphyrin_core(mol, conf_id, ff)
        method = (f"porphyrin core placed from metalloporphyrin X-ray averages (M–N 2.02 Å), "
                  f"periphery by RDKit {ff} with the core fixed; {metal} at the centre")
    else:
        metal = None

    conf = mol.GetConformer(conf_id)
    symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    coords = [tuple(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
    if metal:
        symbols.append(metal)
        coords.append((0.0, 0.0, 0.0))
        extra_bonds = [(n, len(symbols) - 1, "SINGLE") for n in extra_bonds]

    return dict(mol=mol, conf=conf, symbols=symbols, coords=coords, method=method,
                extra_bonds=extra_bonds, removed=set(), ff=ff)


def _force_field(mol, conf_id, ff):
    if ff.startswith("MMFF"):
        props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
        return AllChem.MMFFGetMoleculeForceField(mol, props, confId=conf_id)
    return AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)


def porphyrin_core(mol, conf_id, ff, m_n=2.02, n_ca=1.38, ca_n_ca=106.0, ca_cb=1.44, cb_cb=1.36, ca_cm=1.39):
    """
    Place the 24-atom porphyrin core of a dianionic porphyrin flat in the xy-plane, centred on the origin,
    using average metalloporphyrin bond lengths (Scheidt & Lee, Struct. Bonding 1987), then let the force
    field relax everything else (substituents, hydrogens) with the core held fixed.
    Returns the indices of the four N atoms (to bond to the metal placed at the origin).
    """
    # --- ideal core in the frame of the pyrrole on +x
    half = math.radians(ca_n_ca / 2)
    ca = (m_n + n_ca * math.cos(half), n_ca * math.sin(half))
    cb = (ca[0] + math.sqrt(ca_cb ** 2 - (ca[1] - cb_cb / 2) ** 2), cb_cb / 2)
    # meso carbon on the 45° diagonal (u, u), |Cm - Ca| = ca_cm: 2u² - 2u(ax+ay) + (ax²+ay²-ca_cm²) = 0
    s, q = ca[0] + ca[1], ca[0] ** 2 + ca[1] ** 2 - ca_cm ** 2
    u = (2 * s + math.sqrt(4 * s * s - 8 * q)) / 4
    local = {"N": (m_n, 0.0), "Ca+": ca, "Ca-": (ca[0], -ca[1]), "Cb+": cb, "Cb-": (cb[0], -cb[1]), "Cm": (u, u)}

    def rot(p, k):
        c, s_ = math.cos(k * math.pi / 2), math.sin(k * math.pi / 2)
        return (p[0] * c - p[1] * s_, p[0] * s_ + p[1] * c, 0.0)

    # --- walk the macrocycle to label atoms: N_k, Ca+_k (towards pyrrole k+1), Ca-_k, Cb±_k, meso_k
    ring = mol.GetRingInfo()
    nbrs = lambda i: [a.GetIdx() for a in mol.GetAtomWithIdx(i).GetNeighbors()]
    sym = lambda i: mol.GetAtomWithIdx(i).GetSymbol()
    ns = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "N" and ring.NumAtomRings(a.GetIdx())]
    alpha_of = {n: [c for c in nbrs(n) if sym(c) == "C"] for n in ns}
    alpha_all = {c for v in alpha_of.values() for c in v}
    n_of = {c: n for n, v in alpha_of.items() for c in v}
    meso_of = {c: next(m for m in nbrs(c) if m not in alpha_all and sym(m) == "C"
                       and sum(x in alpha_all for x in nbrs(m)) == 2) for c in alpha_all}
    beta_of = {c: next(b for b in nbrs(c) if b != meso_of[c] and sym(b) == "C") for c in alpha_all}

    target = {}
    n_cur, a_plus = ns[0], alpha_of[ns[0]][0]
    for k in range(4):
        a_minus = next(c for c in alpha_of[n_cur] if c != a_plus)
        for idx, key in ((n_cur, "N"), (a_plus, "Ca+"), (a_minus, "Ca-"),
                         (beta_of[a_plus], "Cb+"), (beta_of[a_minus], "Cb-"), (meso_of[a_plus], "Cm")):
            target[idx] = rot(local[key], k)
        nxt_alpha = next(c for c in nbrs(meso_of[a_plus]) if c in alpha_all and c != a_plus)
        n_cur = n_of[nxt_alpha]
        a_plus = next(c for c in alpha_of[n_cur] if c != nxt_alpha)
    if len(target) != 24:
        raise RuntimeError(f"porphyrin walk labelled {len(target)} atoms, expected 24")

    # --- align the conformer onto the ideal core, pin the core, relax the rest
    from rdkit.Chem import rdMolAlign
    from rdkit.Geometry import Point3D
    ref = Chem.RWMol(mol)
    ref.RemoveAllConformers()
    rc = Chem.Conformer(mol.GetNumAtoms())
    for i, p in target.items():
        rc.SetAtomPosition(i, Point3D(*p))
    ref.AddConformer(rc, assignId=True)
    rdMolAlign.AlignMol(mol, ref, prbCid=conf_id, atomMap=[(i, i) for i in target])
    conf = mol.GetConformer(conf_id)
    for i, p in target.items():
        conf.SetAtomPosition(i, Point3D(*p))
    field = _force_field(mol, conf_id, ff)
    for i in target:
        field.AddFixedPoint(i)
    field.Minimize(maxIts=50000)
    return ns


def check_bonds(built):
    """Compare every bond against RANGES. Returns (rows, outliers, unchecked)."""
    mol, coords, symbols = built["mol"], built["coords"], built["symbols"]
    removed = built["removed"]
    keep = [i for i in range(mol.GetNumAtoms()) if i not in removed]
    idx = {old: new for new, old in enumerate(keep)}

    bonds = []
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i in removed or j in removed:
            continue
        bonds.append((idx[i], idx[j], str(b.GetBondType())))
    bonds += built["extra_bonds"]

    per_type, outliers, unchecked = {}, [], 0
    for i, j, btype in bonds:
        d = math.dist(coords[i], coords[j])
        key = bond_key(symbols[i], symbols[j], btype)
        per_type.setdefault(key, []).append(d)
        rng = RANGES.get(key)
        if rng is None:
            unchecked += 1
        elif not rng[0] <= d <= rng[1]:
            outliers.append((key, i + 1, j + 1, d, rng))
    return per_type, outliers, unchecked


def generic_skeleton(m):
    """Heavy-atom graph with bond orders, aromaticity and charges erased, for matching across forms."""
    m = Chem.RWMol(Chem.RemoveHs(m, sanitize=False))
    for a in m.GetAtoms():
        a.SetIsAromatic(False)
        a.SetFormalCharge(0)
        a.SetNoImplicit(True)
        a.SetNumExplicitHs(0)
    for b in m.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE)
        b.SetIsAromatic(False)
    m.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(m)
    return m


def compare_pubchem(entry, built):
    """Bond-by-bond comparison with PubChem's own 3D conformer (MMFF94s via OMEGA), if one exists."""
    cid = entry.get("cid")
    if not cid or entry.get("post"):
        return None, "n/a"
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
    try:
        sdf = fetch(url, CACHE / "pubchem3d" / f"{cid}.sdf")
    except Exception:
        return None, "no PubChem 3D conformer"
    ref = Chem.MolFromMolBlock(sdf, removeHs=False, sanitize=False)
    if ref is None:
        return None, "no PubChem 3D conformer"

    ours = built["mol"]
    # Match heavy-atom skeletons (bond orders / charges ignored, since PubChem may draw another
    # resonance form). Positions in the H-stripped graphs map back via these lists.
    g_ref, g_ours = generic_skeleton(ref), generic_skeleton(ours)
    if g_ref.GetNumAtoms() != g_ours.GetNumAtoms() or g_ref.GetNumBonds() != g_ours.GetNumBonds():
        return None, "PubChem 3D record is a different form (e.g. closed lactone or another tautomer)"
    matches = g_ref.GetSubstructMatches(g_ours, uniquify=False, maxMatches=5000)
    if not matches:
        return None, "PubChem 3D record has a different skeleton"

    our_heavy = [a.GetIdx() for a in ours.GetAtoms() if a.GetAtomicNum() > 1]
    ref_heavy = [a.GetIdx() for a in ref.GetAtoms() if a.GetAtomicNum() > 1]
    pos = {i: k for k, i in enumerate(our_heavy)}
    ref_conf = ref.GetConformer()
    pairs = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in ours.GetBonds()
             if b.GetBeginAtomIdx() in pos and b.GetEndAtomIdx() in pos]
    ours_len = [math.dist(built["coords"][i], built["coords"][j]) for i, j in pairs]

    # Symmetric molecules have several valid atom mappings (e.g. which side carries the N+ in a
    # rhodamine); use the one that lines up best.
    best = None
    for match in matches:
        diffs = []
        for (i, j), d_ours in zip(pairs, ours_len):
            ri, rj = ref_heavy[match[pos[i]]], ref_heavy[match[pos[j]]]
            diffs.append(abs(d_ours - (ref_conf.GetAtomPosition(ri) - ref_conf.GetAtomPosition(rj)).Length()))
        if best is None or sum(diffs) < sum(best):
            best = diffs
    return best, "ok"


def describe_source(entry):
    src = entry.get("source", "pubchem")
    ref = f"; review ref {entry['ref']}" if entry.get("ref") else ""
    if src == "si":
        return f"original paper SI{ref}"
    if src == "drawing":
        return f"transcribed from the review's Fig. {entry.get('fig', '?')}{ref}"
    return f"PubChem CID {entry.get('cid', '-')}{ref}"


def formula_check(mol, entry):
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
    got = CalcMolFormula(mol).replace("+", "").replace("-", "")
    exp = entry.get("formula")
    if not exp:
        return f"{got} (no reference formula)"
    basis = {"si": "SI HRMS", "drawing": "hand count of the drawing"}.get(entry.get("source"), "PubChem")
    return f"{got} — {'matches' if got == exp else '**MISMATCH** vs'} {exp} ({basis})"


# --------------------------------------------------------------------------- proteins

def build_protein(entry):
    pdb_id = entry["pdb"].upper()
    text = fetch(f"https://files.rcsb.org/download/{pdb_id}.pdb", CACHE / "pdb" / f"{pdb_id}.pdb")
    symbols, coords, residues, hetero = [], [], set(), set()
    for line in text.splitlines():
        rec = line[:6].strip()
        if rec == "ENDMDL":
            break  # first model only
        if rec not in ("ATOM", "HETATM"):
            continue
        resname = line[17:20].strip()
        if rec == "HETATM" and resname not in entry.get("keep_het", []):
            continue  # waters, ions, crystallisation additives
        if line[21] != entry.get("chain", "A"):
            continue
        if line[16] not in (" ", "A"):  # first alternate location only
            continue
        element = line[76:78].strip() or line[12:16].strip()[0]
        element = element[0].upper() + element[1:].lower()
        symbols.append(element)
        coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        residues.add((line[22:27], resname))
        if rec == "HETATM":
            hetero.add(resname)
    n_h = symbols.count("H")
    h_note = f"{n_h} H from the deposited model" if n_h else "no H"
    method = f"RCSB PDB {pdb_id}, chain {entry.get('chain', 'A')}, X-ray ({h_note}; solvent/additives removed)"
    return dict(symbols=symbols, coords=coords, method=method, n_res=len(residues), hetero=sorted(hetero))


# --------------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", help="rebuild only entries whose file name contains this text")
    args = ap.parse_args()

    report = ["# Geometry report", "",
              "Generated by `tools/build_examples.py`. Bond lengths are checked against typical ",
              "crystallographic windows (Allen et al. 1987) and, where available, against PubChem's 3D conformer.", ""]
    summary = []

    for tab, _label in TABS:
        if tab == "basics":
            continue  # made by tools/make_samples.py
        out_dir = SAMPLES / tab
        out_dir.mkdir(parents=True, exist_ok=True)
        index = []
        if not all("param" in e for e in MOLECULES[tab]):  # parametric tabs: see docs/triplet-report.md
            report += [f"## {dict(TABS)[tab]}", ""]

        for entry in MOLECULES[tab]:
            index.append({"file": entry["file"], "name": entry["name"], "description": entry["desc"]})
            if "param" in entry:  # built in the app from a repeat unit (tools/build_wire.py)
                index[-1]["param"] = entry["param"]
                continue
            if args.only and args.only not in entry["file"]:
                continue
            print(f"[{tab}] {entry['name']} ...", end=" ", flush=True)

            if "pdb" in entry:
                built = build_protein(entry)
                comment = f"{entry['name']} | {built['method']}"
                (out_dir / entry["file"]).write_text(xyz_text(built["symbols"], built["coords"], comment), "utf-8")
                print(f"{len(built['symbols'])} atoms, {built['n_res']} residues")
                report += [f"### {entry['name']}", "",
                           f"- Source: {built['method']}",
                           f"- {len(built['symbols'])} atoms, {built['n_res']} residues; "
                           f"hetero groups: {', '.join(built['hetero']) or 'none'}",
                           "- Experimental coordinates, so no geometry check is applied.", ""]
                summary.append((entry["name"], len(built["symbols"]), "X-ray", 0, "-"))
                continue

            built = build_small_cached(entry)
            origin = (f"PubChem CID {entry['cid']}" if entry.get("source", "pubchem") == "pubchem" and entry.get("cid")
                      else f"from {entry.get('source')} (review ref {entry.get('ref', '-')})")
            comment = f"{entry['name']} | {origin} | {built['method']}"
            (out_dir / entry["file"]).write_text(xyz_text(built["symbols"], built["coords"], comment), "utf-8")

            per_type, outliers, unchecked = check_bonds(built)
            diffs, pc_status = compare_pubchem(entry, built)
            n_bonds = sum(len(v) for v in per_type.values())
            print(f"{len(built['symbols'])} atoms, {n_bonds} bonds, {len(outliers)} outside range, "
                  f"PubChem: {pc_status if diffs is None else f'mean |Δ| {statistics.mean(diffs):.3f} Å'}")

            report += [f"### {entry['name']}", "",
                       f"- Structure: {describe_source(entry)}"
                       + (f" ({entry['note']})" if entry.get("note") else ""),
                       f"- Formula: {formula_check(built['mol'], entry)}",
                       f"- Geometry: {built['method']}",
                       f"- {len(built['symbols'])} atoms, {n_bonds} bonds; "
                       f"**{len(outliers)} outside typical range**, {unchecked} without a reference window"]
            if diffs is None:
                report.append(f"- PubChem 3D comparison: {pc_status}")
                pc_cell = "-"
            else:
                pc_cell = f"{statistics.mean(diffs):.3f} / {max(diffs):.3f}"
                report.append(f"- PubChem 3D comparison ({len(diffs)} heavy-atom bonds): "
                              f"mean |Δ| = {statistics.mean(diffs):.3f} Å, max |Δ| = {max(diffs):.3f} Å")
            report += ["", "| Bond | n | mean (Å) | min | max | typical |", "|---|---|---|---|---|---|"]
            for key in sorted(per_type):
                ds = per_type[key]
                rng = RANGES.get(key)
                label = f"{key[0]}–{key[1]} {key[2].lower()}"
                report.append(f"| {label} | {len(ds)} | {statistics.mean(ds):.3f} | {min(ds):.3f} | {max(ds):.3f} | "
                              + (f"{rng[0]:.2f}–{rng[1]:.2f}" if rng else "—") + " |")
            for key, i, j, d, rng in outliers:
                report.append(f"\n> ⚠ atoms {i}–{j} ({key[0]}–{key[1]} {key[2].lower()}): {d:.3f} Å, "
                              f"expected {rng[0]:.2f}–{rng[1]:.2f}")
            report.append("")
            summary.append((entry["name"], len(built["symbols"]), built["ff"], len(outliers), pc_cell))

        (out_dir / "index.json").write_text(json.dumps(index, indent=2, ensure_ascii=False) + "\n", "utf-8")

    # tab order/labels for the app (single source of truth is catalog.TABS)
    tabs = [{"id": t, "label": label} for t, label in TABS]
    (SAMPLES / "tabs.json").write_text(json.dumps(tabs, indent=2) + "\n", "utf-8")

    if not args.only:
        table =["## Summary", "", "| Molecule | Atoms | Method | Bonds outside range | vs PubChem 3D: mean / max |Δ| (Å) |",
                 "|---|---|---|---|---|"]
        table += [f"| {n} | {a} | {m} | {o} | {p} |" for n, a, m, o, p in summary]
        report[5:5] = table + [""]
        REPORT.parent.mkdir(exist_ok=True)
        REPORT.write_text("\n".join(report) + "\n", encoding="utf-8")
        print("wrote", REPORT.relative_to(ROOT))


if __name__ == "__main__":
    main()
