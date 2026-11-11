"""Build "LHCII: chlorophylls and carotenoids" (Self-healing dyes tab): nature's version of a self-healing dye.

LHCII, the main light-harvesting complex of plants, holds 14 chlorophylls and 4 carotenoids in one protein. A
chlorophyll that crosses to its triplet would make singlet oxygen; the carotenoids next to it take that triplet by
triplet-triplet energy transfer (TTET) within nanoseconds and shed it as heat. There is no covalent linker: the protein
holds each carotenoid in van der Waals contact (~3.5 Å) with the chlorophylls it protects.

Source: RCSB PDB 1RWT (spinach LHCII, 2.72 Å, Liu et al., Nature 428, 287, 2004), chain A. Waters, lipids (LHG,
DGD), detergent and the Na+ are left out.

Hydrogens: X-ray maps at this resolution don't show H (one electron each), so the deposited model has none. The
pigments get theirs back: each pigment's bonds, bond orders and H count per atom come from its PDB Chemical Component
Dictionary entry (CLA, CHL, LUT, XAT, NEX), matched by atom name, and RDKit places the H on the crystal heavy atoms
(standard bond lengths and angles). Atoms missing from the crystal (disordered phytol tails) stay missing, and their
neighbours keep their own H count, so a cut tail is not capped with an extra H. The protein stays without H, like the
other proteins in the app.

File order: the protein (drawn faint, as `surroundings`), then each pigment: first its atoms outside the conjugated
system, then its conjugated system with its H, which is its site (where the excited-state glow goes):
  chlorophyll  Mg, the four N and the 20 carbons of the macrocycle (names MG, N?, C1?-C4?, CH?)
  carotenoid   the polyene chain: carbons outside the end rings with two or three carbon neighbours (no methyls)
"""

import shlex

import numpy as np
from rdkit import Chem
from rdkit.Geometry import Point3D

PIGMENTS = {"CLA": ("Chl a", "chlorophyll"), "CHL": ("Chl b", "chlorophyll"),
            "LUT": ("Lutein", "quencher"), "XAT": ("Violaxanthin", "quencher"), "NEX": ("Neoxanthin", "quencher")}
MACROCYCLE = {"MG", "CHA", "CHB", "CHC", "CHD"} | {f"N{r}" for r in "ABCD"} | {f"C{k}{r}" for k in "1234" for r in "ABCD"}
BOND = 1.9  # Å, heavy-atom bond cut-off within a pigment (C–C, C–O; the longest is C–C 1.55)


def _element(line):
    e = line[76:78].strip() or line[12:16].strip()[0]
    return e[0].upper() + e[1:].lower()


def _ccd(text):
    """Atoms {name: (element, charge, leaving)} and bonds [(name1, name2, order, metal)] of a Chemical Component Dictionary
    mmCIF entry."""
    lines = text.splitlines()

    def loop(prefix):
        i = next(k for k, line in enumerate(lines) if line.startswith(prefix))
        keys = []
        while lines[i].startswith(prefix):
            keys.append(lines[i].strip().split(".", 1)[1])
            i += 1
        rows = []
        while i < len(lines) and lines[i].strip() and not lines[i].startswith(("#", "_", "loop_")):
            rows.append(dict(zip(keys, shlex.split(lines[i]))))
            i += 1
        return rows

    atoms = {r["atom_id"]: (r["type_symbol"].capitalize(), int(r["charge"]), r.get("pdbx_leaving_atom_flag") == "Y")
             for r in loop("_chem_comp_atom.")}
    bonds = [(r["atom_id_1"], r["atom_id_2"], r["value_order"], r.get("pdbx_metal_coordination_flag") == "Y")
             for r in loop("_chem_comp_bond.")]
    return atoms, bonds


ORDER = {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE, "TRIP": Chem.BondType.TRIPLE}


def with_hydrogens(atoms, ccd_text, anchors=()):
    """Crystal heavy atoms [(name, element, xyz)] -> the same, then their H: [(name, element, xyz, parent or None)].

    H marked as leaving in the dictionary (lost when the residue joins a chain, e.g. a chromophore's N-terminal H2)
    are not added. anchors: [(atom name, element, xyz)] atoms outside the residue bonded to it (the peptide C and N
    of the neighbouring residues), used only to aim the H and then dropped."""
    ref_atoms, ref_bonds = _ccd(ccd_text)
    n_h = {}
    for a, b, _, _ in ref_bonds:
        for x, y in ((a, b), (b, a)):
            if ref_atoms[y][0] == "H" and not ref_atoms[y][2]:
                n_h[x] = n_h.get(x, 0) + 1
    idx = {name: i for i, (name, _, _) in enumerate(atoms)}
    unknown = [name for name in idx if name not in ref_atoms]
    assert not unknown, f"atoms not in the dictionary: {unknown}"
    mol = Chem.RWMol()
    for name, element, _ in atoms:
        a = Chem.Atom(element)
        a.SetFormalCharge(0 if element == "Mg" else ref_atoms[name][1])  # Mg: a bare ion here (its N bonds are dative)
        a.SetNoImplicit(True)
        a.SetNumExplicitHs(n_h.get(name, 0))
        mol.AddAtom(a)
    for a, b, order, metal in ref_bonds:
        if not metal and a in idx and b in idx:
            mol.AddBond(idx[a], idx[b], ORDER[order])
    extra = []
    for name, element, p in anchors:
        a = Chem.Atom(element)
        a.SetNoImplicit(True)
        k = mol.AddAtom(a)
        mol.AddBond(idx[name], k, Chem.BondType.SINGLE)
        extra.append(p)
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, p in enumerate([p for _, _, p in atoms] + extra):
        conf.SetAtomPosition(i, Point3D(*p))
    mol.AddConformer(conf)
    Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
    mol = Chem.AddHs(mol, addCoords=True)
    pos = mol.GetConformer().GetPositions()
    first_h = len(atoms) + len(extra)
    out = [(name, element, tuple(p), None) for name, element, p in atoms]
    for a in list(mol.GetAtoms())[first_h:]:
        out.append((f"H{a.GetIdx()}", "H", tuple(pos[a.GetIdx()]), a.GetNeighbors()[0].GetIdx()))
    return out


def _polyene(atoms):
    """Indices of a carotenoid's conjugated chain: non-ring carbons with >= 2 carbon neighbours."""
    xyz = np.array([p for _, _, p in atoms])
    is_c = [e == "C" for _, e, _ in atoms]
    d = np.linalg.norm(xyz[:, None] - xyz[None], axis=2)
    nbr = [[j for j in range(len(atoms)) if j != i and d[i, j] < BOND] for i in range(len(atoms))]

    def in_ring(i):
        """Does a bond of i lie on a ring: can a neighbour reach i again without that bond (within 7 bonds)?"""
        for b in nbr[i]:
            seen, frontier = {b}, [b]
            for _ in range(7):
                frontier = [c for a in frontier for c in nbr[a] if not (a == b and c == i) and c not in seen]
                if i in frontier:
                    return True
                seen.update(frontier)
        return False

    return [i for i in range(len(atoms)) if is_c[i] and not in_ring(i) and sum(is_c[j] for j in nbr[i]) >= 2]


def build_lhcii(text, ccd, chain="A", only=None, protein_too=True):
    """ccd(code): the Chemical Component Dictionary mmCIF of a pigment, for its H. only: pigment labels to keep (e.g.
    ["Chl a610", "Lutein 620"]), None for all; protein_too: keep the protein."""
    protein, pigments = [], {}
    for line in text.splitlines():
        rec = line[:6].strip()
        if rec == "ENDMDL":
            break
        if rec not in ("ATOM", "HETATM") or line[21] != chain or line[16] not in (" ", "A"):
            continue
        p = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        resname = line[17:20].strip()
        if rec == "ATOM":
            if protein_too:
                protein.append((_element(line), p))
        elif resname in PIGMENTS:
            pigments.setdefault((resname, int(line[22:26])), []).append((line[12:16].strip(), _element(line), p))

    if only:
        missing = set(only) - {_label(k) for k in pigments}
        assert not missing, f"no pigment {missing} in chain {chain}"
        pigments = {k: v for k, v in pigments.items() if _label(k) in only}
    symbols = [s for s, _ in protein]
    coords = [p for _, p in protein]
    blocks, chl_core = [], {}
    order = sorted(pigments, key=lambda k: (PIGMENTS[k[0]][1] != "chlorophyll", k[1]))  # chlorophylls, then carotenoids
    for key in order:
        atoms = pigments[key]
        name, kind = PIGMENTS[key[0]]
        if kind == "chlorophyll":
            core = [i for i, (n, _, _) in enumerate(atoms) if n in MACROCYCLE]
            assert len(core) == 25, (key, len(core))
            label = f"Chl {name[-1]}{key[1]}"  # "Chl a610": the numbering used in the LHCII literature
        else:
            core = _polyene(atoms)
            assert 16 <= len(core) <= 20, (key, len(core))
            label = f"{name} {key[1]}"
        full = with_hydrogens(atoms, ccd(key[0]))  # heavy atoms keep their indices; H follow, with their parent
        hs = {}
        for j, (_, _, _, parent) in enumerate(full):
            if parent is not None:
                hs.setdefault(parent, []).append(j)
        n_core = len(core) + sum(len(hs.get(i, [])) for i in core)
        for i in [i for i in range(len(atoms)) if i not in core] + core:
            for j in [i] + hs.get(i, []):
                symbols.append(full[j][1])
                coords.append(full[j][2])
        xyz = np.array([atoms[i][2] for i in core])
        if kind == "chlorophyll":
            chl_core[label] = xyz
        blocks.append((label, kind, n_core, len(symbols), xyz))  # len(symbols) = where this site ends

    # sites as the app reads them: the `last` atoms before the final `skip` atoms of the file
    total = len(symbols)
    sites = [dict(label=lb, kind=kd, last=n, **({"skip": total - end} if total - end else {}),
                  **({"on": True} if lb == "Chl a610" else {}))  # the lowest-energy chlorophyll starts lit
             for lb, kd, n, end, _ in blocks]

    # which chlorophylls each carotenoid touches (closest conjugated atoms), for the report and the description
    contacts = {}
    for label, kind, _, _, xyz in blocks:
        if kind != "quencher":
            continue
        d = sorted((float(np.min(np.linalg.norm(c[:, None] - xyz[None], axis=2))), lb) for lb, c in chl_core.items())
        contacts[label] = [(lb, round(x, 1)) for x, lb in d[:3]]

    n_chl = sum(PIGMENTS[k[0]][1] == "chlorophyll" for k in pigments)
    what = (f"pigment H placed from the PDB Chemical Component Dictionary; protein drawn faint (no H), {n_chl} chlorophylls and {len(pigments) - n_chl} carotenoids; waters, lipids, "
            f"detergent removed") if protein_too and not only else (
            f"only {' and '.join(_label(k) for k in order)}, at their crystal positions (protein and the rest removed); "
            f"H placed from the PDB Chemical Component Dictionary")
    method = f"RCSB PDB 1RWT (spinach LHCII, 2.72 Å, Liu et al. 2004), chain {chain}, X-ray: {what}"
    return dict(symbols=symbols, coords=coords, sites=sites, surroundings=len(protein), method=method,
                contacts=contacts)


def _label(key):
    """("CLA", 610) -> "Chl a610"; ("LUT", 620) -> "Lutein 620"."""
    name, kind = PIGMENTS[key[0]]
    return f"Chl {name[-1]}{key[1]}" if kind == "chlorophyll" else f"{name} {key[1]}"
