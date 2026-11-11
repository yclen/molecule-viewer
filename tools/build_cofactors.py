"""Proteins (or bare pigments) from the PDB with their photoactive cofactors as sites (Interesting photophysics tab).

A catalog entry's `cofactors` dict says what to take:
    pdb        PDB ID; assembly: the biological-assembly number, for complexes whose deposited file holds only part of
               the complex (e.g. LH2: 3 of its 9 subunit pairs). Copies of the assembly become separate chains.
    chains     chains to keep (default all)
    pigments   {het code: {name, kind, core, color}}: the cofactors that become sites. core says which atoms glow:
                 "macrocycle"  a (bacterio)chlorophyll's Mg, 4 N and 20 ring C
                 "polyene"     a carotenoid's conjugated chain
                 "all"         the whole residue (with its H, the glow finds the conjugated part by itself)
    keep       other het codes kept, drawn faint with the protein (Ca2+ ions, modified residues of the chain)
    protein    False: no protein at all (only the cofactors)
    labels     "lh2": name the LH2 bacteriochlorophylls B800 1-9 / B850 α1-9, β1-9 by ring and position

Every cofactor gets its H back from the PDB Chemical Component Dictionary (see build_lhcii.with_hydrogens): X-ray
models don't show H. A cofactor built into the protein chain (a fluorescent protein's chromophore) is bonded to its
neighbours; those bonds are found by distance and used to aim its H. The protein has no H, like the other proteins.

File order: protein + kept het (drawn faint: `surroundings`), then each cofactor: its non-glowing atoms, then its
glowing atoms, which are its site.
"""

import math

import numpy as np

from build_lhcii import MACROCYCLE, _element, _polyene, with_hydrogens

LINK = 1.75  # Å: a protein atom this close to a cofactor atom is covalently bonded to it (peptide C-N 1.33)


def _read(text, chains, wanted, keep):
    """Protein atoms [(element, xyz)] and residues {(chain, code, number): [(name, element, xyz)]} for the cofactors."""
    protein, residues, model = [], {}, 1
    for line in text.splitlines():
        rec = line[:6].strip()
        if rec == "MODEL":
            model = int(line.split()[1])
        if rec not in ("ATOM", "HETATM") or line[16] not in (" ", "A"):
            continue
        chain = line[21] + (str(model) if model > 1 else "")  # assembly copies: A, A2, A3
        if chains and line[21] not in chains:
            continue
        p = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        code = line[17:20].strip()
        if rec == "ATOM" or code in keep:
            protein.append((_element(line), p))
        elif code in wanted:
            residues.setdefault((chain, code, int(line[22:26])), []).append((line[12:16].strip(), _element(line), p))
    return protein, residues


def _lh2_labels(residues, spec):
    """B800 / B850 by ring (height along the 9-fold axis) and position around it; carotenoids numbered around it."""
    mg = {k: np.array(next(p for n, _, p in v if n == "MG")) for k, v in residues.items() if k[1] == "BCL"}
    pts = np.array(list(mg.values()))
    centre = pts.mean(0)
    axis = np.linalg.svd(pts - centre)[2][2]  # the rings are flat: the axis is the direction of least spread
    height = {k: float((p - centre) @ axis) for k, p in mg.items()}
    ref = np.linalg.svd(pts - centre)[2][0]
    side = np.cross(axis, ref)

    def angle(p):
        v = p - centre
        return math.atan2(float(v @ side), float(v @ ref)) % (2 * math.pi)

    hs = sorted(height.values())
    gap = max(range(len(hs) - 1), key=lambda i: hs[i + 1] - hs[i])
    low = set(hs[:gap + 1])
    b800_low = len(low) == 9
    labels = {}
    b800 = sorted((k for k in mg if (height[k] in low) == b800_low), key=lambda k: angle(mg[k]))
    for n, k in enumerate(b800, 1):
        labels[k] = f"B800 {n}"
    for sub in ("α", "β"):
        chains = spec["alpha"] if sub == "α" else spec["beta"]
        ring = sorted((k for k in mg if k not in labels and k[0][0] in chains), key=lambda k: angle(mg[k]))
        for n, k in enumerate(ring, 1):
            labels[k] = f"B850 {sub}{n}"
    cars = sorted((k for k in residues if k[1] != "BCL"),
                  key=lambda k: angle(np.mean([p for _, _, p in residues[k]], axis=0)))
    name = spec["pigments"][cars[0][1]]["name"] if cars else ""
    for n, k in enumerate(cars, 1):
        labels[k] = f"{name} {n}"
    order = b800 + [k for sub in "αβ" for k in sorted(mg, key=lambda k: angle(mg[k]))
                    if labels[k].startswith(f"B850 {sub}")] + cars
    return labels, order


def build_cofactors(text, ccd, spec):
    pig = spec["pigments"]
    protein, residues = _read(text, spec.get("chains"), set(pig), set(spec.get("keep", ())))
    if spec.get("protein") is False:
        protein = []
    if spec.get("labels") == "lh2":
        labels, order = _lh2_labels(residues, spec)
    else:
        multi = len({k[0] for k in residues}) > 1
        labels = {k: (f"{pig[k[1]]['name']} {k[0] if multi else ''}{k[2]}" if pig[k[1]].get("number", True)
                      else pig[k[1]]["name"]) for k in residues}
        kinds = list(dict.fromkeys(v["kind"] for v in pig.values()))
        order = sorted(residues, key=lambda k: (kinds.index(pig[k[1]]["kind"]), k[0], k[2]))
    if spec.get("only"):
        order = [k for k in order if labels[k] in spec["only"]]

    prot_xyz = np.array([p for _, p in protein]) if protein else np.zeros((0, 3))
    prot_el = [e for e, _ in protein]
    symbols = [e for e, _ in protein]
    coords = [p for _, p in protein]
    blocks, links = [], 0
    for key in order:
        atoms = residues[key]
        cfg = pig[key[1]]
        # covalent neighbours in the protein (a chromophore's peptide bonds), to aim its H
        anchors = []
        if len(prot_xyz):
            for name, _, p in atoms:
                d = np.linalg.norm(prot_xyz - np.array(p), axis=1)
                for j in np.where(d < LINK)[0]:
                    anchors.append((name, prot_el[j], tuple(prot_xyz[j])))
        links += len(anchors)
        full = with_hydrogens(atoms, ccd(key[1]), anchors)
        core = {"macrocycle": lambda: [i for i, (n, _, _) in enumerate(atoms) if n in MACROCYCLE],
                "polyene": lambda: _polyene(atoms),
                "all": lambda: list(range(len(atoms)))}[cfg.get("core", "all")]()
        hs = {}
        for j, (_, _, _, parent) in enumerate(full):
            if parent is not None:
                hs.setdefault(parent, []).append(j)
        n_core = len(core) + sum(len(hs.get(i, [])) for i in core)
        for i in [i for i in range(len(atoms)) if i not in core] + core:
            for j in [i] + hs.get(i, []):
                symbols.append(full[j][1])
                coords.append(full[j][2])
        blocks.append((labels[key], cfg["kind"], cfg.get("color"), n_core, len(symbols),
                       np.array([atoms[i][2] for i in core])))

    total = len(symbols)
    sites = []
    for label, kind, color, n, end, _ in blocks:
        site = dict(label=label, kind=kind, last=n)
        if total - end:
            site["skip"] = total - end
        if color:
            site["color"] = color
        if label in spec.get("on", ()):
            site["on"] = True
        sites.append(site)

    # closest pairs of glowing atoms between sites (for the report and the info text)
    contacts = []
    for a in range(len(blocks)):
        for b in range(a + 1, len(blocks)):
            d = float(np.min(np.linalg.norm(blocks[a][5][:, None] - blocks[b][5][None], axis=2)))
            contacts.append((round(d, 1), blocks[a][0], blocks[b][0]))
    contacts.sort()

    n_h = sum(1 for s in symbols[len(protein):] if s == "H")
    what = (f"{len(protein)} protein atoms drawn faint (no H)" if protein else "no protein") + \
        f"; {len(blocks)} cofactor sites with {n_h} H from the PDB Chemical Component Dictionary"
    method = f"RCSB PDB {spec['pdb']}" + (f" assembly {spec['assembly']}" if spec.get("assembly") else "") + \
        f", X-ray: {what}; waters, lipids, detergents and buffer removed"
    return dict(symbols=symbols, coords=coords, sites=sites, surroundings=len(protein), method=method,
                contacts=contacts, links=links)
