"""A protein model set into a patch of lipid bilayer. A static picture, not a simulation.

The bilayer is the last frame of a 270 ns MD run of 128 POPC lipids (Berger united-atom lipids, so no H on carbon,
like the protein models), from Ollila, Ferreira and Topgaard, Zenodo record 13279. Here it is:

  1. made whole again (the MD box wraps lipids across its edges) and tiled to cover the patch;
  2. turned so its normal is y and its midplane y = 0, the frame tools/build_voltron.py puts the protein in;
  3. trimmed to a disc around the protein, and stripped of every lipid with an atom too close to the protein.

Nothing is relaxed afterwards, so the lipids next to the protein leave small gaps where the deleted ones were. Water is
left out. Called by tools/build_examples.py for catalog entries with `membrane`. The lipids come first in the file, so
the protein's sites, counted from the end, are unchanged.
"""

import numpy as np
from scipy.spatial import cKDTree

CLASH = 2.8  # Å, a lipid with any atom closer than this to a protein atom is removed


def read_bilayer(gro_text, resname="PLA"):
    """[(symbols, xyz)] per lipid, made whole, and the box (Å). GROMACS .gro, positions in nm."""
    lines = gro_text.splitlines()
    n = int(lines[1])
    box = np.array([float(v) for v in lines[2 + n].split()[:3]]) * 10
    lipids, current, key = [], None, None
    for row in lines[2:2 + n]:
        if row[5:10].strip() != resname:
            continue
        k = int(row[0:5])
        if k != key:
            key, current = k, ([], [])
            lipids.append(current)
        current[0].append(row[10:15].strip()[0])
        current[1].append([float(row[20:28]) * 10, float(row[28:36]) * 10, float(row[36:44]) * 10])
    whole = []
    for symbols, xyz in lipids:
        xyz = np.array(xyz)
        anchor = xyz[symbols.index("P")]
        xyz -= box * np.round((xyz - anchor) / box)  # each atom to its image nearest the phosphate
        whole.append((symbols, xyz))
    return whole, box


def build_membrane(gro_text, protein_symbols, protein_coords, radius):
    """Lipids around a protein that sits in the y-normal, y = 0 midplane frame. Returns symbols, coords, counts."""
    lipids, box = read_bilayer(gro_text)
    midplane = np.mean([xyz[:, 2].mean() for _, xyz in lipids])

    protein = np.asarray(protein_coords, dtype=float)
    heavy = protein[[i for i, s in enumerate(protein_symbols) if s != "H"]]
    core = heavy[np.abs(heavy[:, 1]) < 25]  # the part in the membrane sets the disc's centre
    centre = core.mean(axis=0)
    tree = cKDTree(heavy)

    tiles = int(np.ceil(radius / min(box[:2]))) + 1
    symbols, coords, kept, clashed = [], [], 0, 0
    for i in range(-tiles, tiles + 1):
        for j in range(-tiles, tiles + 1):
            shift = np.array([(i - 0.5) * box[0], (j - 0.5) * box[1], -midplane])  # tile (0, 0) centred on 0
            for names, xyz in lipids:
                moved = xyz + shift
                placed = np.column_stack([moved[:, 0], moved[:, 2], -moved[:, 1]])  # bilayer normal z -> y
                placed[:, [0, 2]] += centre[[0, 2]]
                head = placed[names.index("P")]
                if np.hypot(head[0] - centre[0], head[2] - centre[2]) > radius:
                    continue
                if tree.query(placed, k=1)[0].min() < CLASH:
                    clashed += 1
                    continue
                kept += 1
                symbols += names
                coords += [tuple(p) for p in placed]
    return dict(symbols=symbols, coords=coords, lipids=kept, removed=clashed)
