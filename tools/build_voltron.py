"""Rhodopsin voltage indicators: a rhodopsin and a fluorophore domain fused into one chain. Composite models.

No structure of a whole indicator exists. Each domain here is a crystal structure, kept rigid; their relative pose is
invented, and only the residues at the junction are rebuilt:

  1. the rhodopsin (3AM6, chain A, with its retinal, Asp92 changed to Asn as in Ace2N) is placed as OPM places it in
     a membrane (Orientations of Proteins in Membranes), turned so the membrane normal is y, the midplane y = 0 and
     the cytoplasm below. The fluorophore domain is fused to its C-terminus, which microbial rhodopsins hold in the
     cytoplasm;
  2. random poses of the fluorophore domain are tried, keeping it below the lipid phosphates, clear of the rhodopsin and
     within reach of the junction, with a clear path between the chain ends. Poses are ranked by how near the
     chromophore comes to the retinal, since the design depends on FRET between them;
  3. RDKit regrows the junction residues between the crystal anchors (MMFF94s, the anchors held), going down the
     ranking until a junction fits without touching either domain. Their sequence is the patent's (US 11,598,792):
     Voltron joins ...IFK to IGTGF... with no linker at all, the five residues removed being KM from the rhodopsin
     and MAE from HaloTag.

Called by tools/build_examples.py for catalog entries with `fusion`. File order: rhodopsin, fluorophore domain,
junction, retinal, chromophore, so the two sites are found from the end of the file.
"""

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from rdkit.Geometry import Point3D
from scipy.spatial import cKDTree

SEED = 525
POSES = 20000
TRIES = 12  # best poses to try growing the junction in, and seeds per pose
CLEAR = 3.0  # Å, closest allowed heavy-atom distance between the two domains
PHOSPHATES = 3.0  # Å, how far outside OPM's hydrophobic boundary the lipid phosphates sit (POPC: ~18.5 vs 15.6)
REACH = 2.8  # Å per peptide bond the junction has to span; a fully stretched residue is 3.8
LAST_RHODOPSIN = 225  # crystal residues kept; the junction regrows the rest
CONTACT = 2.8  # Å, junction heavy atom to domain heavy atom closer than this counts as a close contact
LONGEST_BOND = 1.7  # Å, C/N/O bonds in the junction; a longer one means the chain was stretched to fit

# backbone-first residues, N to C, all L
RESIDUES = {
    "F": "N[C@@H](Cc1ccccc1)C(=O)", "G": "NCC(=O)", "I": "N[C@@H]([C@@H](C)CC)C(=O)", "K": "N[C@@H](CCCCN)C(=O)",
    "L": "N[C@@H](CC(C)C)C(=O)", "M": "N[C@@H](CCSC)C(=O)", "P": "N1[C@@H](CCC1)C(=O)",
    "R": "N[C@@H](CCCNC(=N)N)C(=O)", "S": "N[C@@H](CO)C(=O)", "T": "N[C@@H]([C@@H](C)O)C(=O)",
}


def read_chain(text, chain, keep_het=()):
    """Atoms of one chain in file order (first model, first alternate location, hydrogens dropped)."""
    atoms = []
    for line in text.splitlines():
        rec = line[:6].strip()
        if rec == "ENDMDL":
            break
        if rec not in ("ATOM", "HETATM") or line[21] != chain or line[16] not in (" ", "A"):
            continue
        res = line[17:20].strip()
        if rec == "HETATM" and res not in keep_het:
            continue
        name = line[12:16].strip()
        element = (line[76:78].strip() or name[0]).capitalize()
        if element in ("H", "D"):
            continue
        atoms.append(dict(name=name, res=res, seq=int(line[22:26]), element=element,
                          xyz=(float(line[30:38]), float(line[38:46]), float(line[46:54]))))
    return atoms


def xyz_of(atoms):
    return np.array([a["xyz"] for a in atoms], dtype=float)


def find(atoms, seq, name):
    return next(a for a in atoms if a["seq"] == seq and a["name"] == name)


def membrane_frame(rho, opm_text):
    """OPM's placement of this chain, as a function crystal xyz -> model xyz, and the hydrophobic half-thickness.

    OPM puts the membrane normal along z with the midplane at z = 0. Here the normal is y (the viewer's vertical), with
    the cytoplasm, where the C-terminus is, below.
    """
    half = next(float(line.split(":")[1]) for line in opm_text.splitlines() if "1/2 of bilayer thickness" in line)
    opm = {(a["seq"], a["name"]): a["xyz"] for a in read_chain(opm_text, "A") if a["name"] == "CA"}
    pairs = [(a["xyz"], opm[(a["seq"], a["name"])]) for a in rho if (a["seq"], a["name"]) in opm]
    src, dst = (np.array(p, dtype=float) for p in zip(*pairs))
    sc, dc = src.mean(axis=0), dst.mean(axis=0)
    u, _, vt = np.linalg.svd((src - sc).T @ (dst - dc))  # Kabsch: crystal -> OPM
    flip = np.diag([1.0, 1.0, np.sign(np.linalg.det(vt.T @ u.T))])
    rot = vt.T @ flip @ u.T
    if np.sqrt((((src - sc) @ rot.T + dc - dst) ** 2).sum(axis=1).mean()) > 0.5:
        raise RuntimeError("OPM's chain A does not superpose on the crystal's")
    to_y = np.array([[1.0, 0, 0], [0, 0, 1.0], [0, -1.0, 0]])  # (x, y, z) -> (x, z, -y)
    if ((np.array(rho[-1]["xyz"]) - sc) @ rot.T + dc)[2] > 0:
        to_y = to_y @ np.diag([1.0, -1.0, -1.0])  # C-terminus on the -y side
    return (lambda xyz: ((np.asarray(xyz, dtype=float) - sc) @ rot.T + dc) @ to_y.T), half


def random_rotation(rng):
    q = rng.normal(size=4)
    w, x, y, z = q / np.linalg.norm(q)
    return np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
                     [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
                     [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)]])


def dock(rho_xyz, rho_end, retinal_centre, tag_xyz, tag_start, chromophore, reach, floor):
    """Poses of the tag that fit, nearest the retinal first: [(distance, rotation, translation)], and how many fit.

    floor: the tag has to stay below this y, out of the membrane.
    """
    tree = cKDTree(rho_xyz)
    rng = np.random.default_rng(SEED)
    local = tag_xyz - tag_start
    local_chromophore = chromophore - tag_start
    path = np.linspace(0.25, 0.75, 5)[:, None]
    poses = []
    for _ in range(POSES):
        rot = random_rotation(rng)
        direction = rng.normal(size=3)
        direction /= np.linalg.norm(direction)
        target = rho_end + direction * rng.uniform(3.8, reach)
        moved = local @ rot.T + target
        if moved[:, 1].max() > floor:
            continue
        if tree.query(moved, k=1)[0].min() < CLEAR:
            continue
        line = rho_end + (target - rho_end) * path  # the junction needs a way through
        if tree.query(line, k=1)[0].min() < CLEAR or cKDTree(moved).query(line, k=1)[0].min() < CLEAR:
            continue
        d = float(np.linalg.norm(local_chromophore @ rot.T + target - retinal_centre))
        poses.append((d, rot, target - rot @ tag_start))
    if not poses:
        raise RuntimeError("no pose keeps the domains apart within reach of the junction")
    poses.sort(key=lambda pose: pose[0])
    return poses, len(poses)


def grow_junction(start, end, sequence, around):
    """Heavy atoms of `sequence` joining a chain end (CA, C, O) to the next chain start (N, CA).

    Returns ([(element, xyz)], close contacts with `around`), or None if every seed leaves the chain stretched.
    Seeds are ranked by contacts, then energy.
    """
    smiles = "CC(=O)" + "".join(RESIDUES[r] for r in sequence) + "NC"
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    ends = {heavy[0]: start["CA"], heavy[1]: start["C"], heavy[2]: start["O"],
            heavy[-2]: end["N"], heavy[-1]: end["CA"]}
    grown = [i for i in heavy if i not in ends]
    tree = cKDTree(around)
    best = None
    for seed in range(1, TRIES + 1):
        trial = Chem.Mol(mol)
        coord_map = {i: Point3D(*map(float, p)) for i, p in ends.items()}
        if AllChem.EmbedMolecule(trial, coordMap=coord_map, randomSeed=seed, useRandomCoords=True) < 0:
            continue
        ref = Chem.Mol(trial)
        for i, p in coord_map.items():
            ref.GetConformer().SetAtomPosition(i, p)
        rdMolAlign.AlignMol(trial, ref, atomMap=[(i, i) for i in ends])
        conf = trial.GetConformer()
        for i, p in coord_map.items():
            conf.SetAtomPosition(i, p)
        props = AllChem.MMFFGetMoleculeProperties(trial, mmffVariant="MMFF94s")
        field = AllChem.MMFFGetMoleculeForceField(trial, props)
        for i in ends:
            field.AddFixedPoint(i)
        field.Initialize()
        for _ in range(5):
            if field.Minimize(maxIts=2000) == 0:
                break
        pos = trial.GetConformer().GetPositions()
        bonds = [np.linalg.norm(pos[b.GetBeginAtomIdx()] - pos[b.GetEndAtomIdx()]) for b in trial.GetBonds()
                 if {b.GetBeginAtom().GetSymbol(), b.GetEndAtom().GetSymbol()} <= {"C", "N", "O"}]
        if max(bonds) > LONGEST_BOND:
            continue
        close = int((tree.query(pos[grown], k=1)[0] < CONTACT).sum())
        score = (close, field.CalcEnergy())
        if best is None or score < best[0]:
            best = (score, trial)
    if best is None:
        return None
    (close, _), trial = best
    pos = trial.GetConformer().GetPositions()
    return [(trial.GetAtomWithIdx(i).GetSymbol(), tuple(pos[i])) for i in grown], close


def build_fusion(rhodopsin_text, opm_text, entry, tag):
    """tag: dict(symbols, coords, chromophore=[indices], start=dict(N, CA), name), in the tag's own frame.

    opm_text: OPM's copy of the rhodopsin, which fixes where the membrane is.
    """
    sequence = entry["fusion"]["junction"]
    atoms = read_chain(rhodopsin_text, "A", keep_het=("RET",))
    whole = [a for a in atoms if a["res"] != "RET"]
    retinal = [a for a in atoms if a["res"] == "RET"]
    rho = [a for a in whole if a["seq"] <= LAST_RHODOPSIN]
    donor = find(rho, 92, "OD2")
    if donor["res"] != "ASP":
        raise RuntimeError("expected Asp92 in the rhodopsin")
    donor["element"] = "N"  # D92N: Asn's side chain is Asp's with one O swapped for N

    place, half = membrane_frame(whole, opm_text)
    rho_xyz, retinal_xyz = place(xyz_of(rho)), place(xyz_of(retinal))
    chain_end = {n: place(find(rho, LAST_RHODOPSIN, n)["xyz"]) for n in ("CA", "C", "O")}

    tag_xyz = np.array(tag["coords"], dtype=float)
    tag_symbols = tag["symbols"]
    heavy = [i for i, s in enumerate(tag_symbols) if s != "H"]
    chromophore = tag_xyz[[i for i in tag["chromophore"] if tag_symbols[i] != "H"]].mean(axis=0)
    reach = REACH * (len(sequence) + 1)
    poses, fits = dock(np.vstack([rho_xyz, retinal_xyz]), chain_end["C"], retinal_xyz.mean(axis=0), tag_xyz[heavy],
                       np.array(tag["start"]["N"]), chromophore, reach, -(half + PHOSPHATES))

    # the residues the junction is bonded to are next to it by design, so they don't count as contacts
    clear_rho = rho_xyz[[i for i, a in enumerate(rho) if a["seq"] != LAST_RHODOPSIN]]
    clear_tag = [i for i in heavy if i not in tag["anchor"]]
    chosen = None
    for rank, (distance, turn, shift) in enumerate(poses[:TRIES]):
        tag_placed = tag_xyz @ turn.T + shift
        start = {n: np.array(p) @ turn.T + shift for n, p in tag["start"].items()}
        grown = grow_junction(chain_end, start, sequence, np.vstack([clear_rho, tag_placed[clear_tag]]))
        if grown and (chosen is None or grown[1] < chosen[0][1]):
            chosen = (grown, distance, tag_placed, rank)
        if grown and grown[1] == 0:
            break
    if chosen is None:
        raise RuntimeError("could not grow the junction between the chain ends")
    (junction, close), distance, tag_placed, rank = chosen

    lit = set(tag["chromophore"])
    inner = [i for i in range(len(tag_symbols)) if i not in lit]
    symbols = [a["element"] for a in rho] + [tag_symbols[i] for i in inner] + [s for s, _ in junction]
    coords = list(rho_xyz) + [tag_placed[i] for i in inner] + [np.array(p) for _, p in junction]
    symbols += [a["element"] for a in retinal] + [tag_symbols[i] for i in tag["chromophore"]]
    coords += list(retinal_xyz) + [tag_placed[i] for i in tag["chromophore"]]

    method = (f"Composite model, not a structure of this protein. Rhodopsin: RCSB PDB 3AM6 chain A with retinal, "
              f"X-ray 3.20 Å, D92N, placed as OPM places it (membrane normal y, midplane y = 0, hydrophobic "
              f"half-thickness {half:.1f} Å). {tag['name']}. Domains rigid; their pose is one of {fits} that fit, "
              f"chromophore to retinal {distance:.1f} Å. Junction {sequence} (US 11,598,792) rebuilt by RDKit "
              f"MMFF94s between crystal anchors")
    return dict(symbols=symbols, coords=[tuple(p) for p in coords], method=method, distance=distance, fits=fits,
                contacts=close, n_chromophore=len(tag["chromophore"]))


def mneon_domain(text, first=5):
    """mNeonGreen (5LTR) from residue `first` on: the junction regrows what comes before."""
    atoms = [a for a in read_chain(text, "A", keep_het=("CR2", "KCX")) if a["seq"] >= first]
    return dict(symbols=[a["element"] for a in atoms], coords=[a["xyz"] for a in atoms],
                chromophore=[i for i, a in enumerate(atoms) if a["res"] == "CR2"],
                anchor={i for i, a in enumerate(atoms) if a["seq"] == first},
                start={n: find(atoms, first, n)["xyz"] for n in ("N", "CA")},
                name="mNeonGreen: RCSB PDB 5LTR chain A, X-ray 1.21 Å")


def halo_domain(conjugate, text, first=7):
    """HaloTag7 with JF525 (from build_conjugate), from residue `first` on: the junction regrows what comes before."""
    atoms = read_chain(text, "A")
    dye = conjugate["dye"]
    if len(atoms) != dye[0]:
        raise RuntimeError("the conjugate's protein atoms are not the chain's")
    keep = [i for i, a in enumerate(atoms) if a["seq"] >= first]
    order = keep + dye
    return dict(symbols=[conjugate["symbols"][i] for i in order], coords=[conjugate["coords"][i] for i in order],
                chromophore=list(range(len(keep), len(order))),
                anchor={k for k, i in enumerate(keep) if atoms[i]["seq"] == first},
                start={n: find(atoms, first, n)["xyz"] for n in ("N", "CA")},
                name="HaloTag7 with JF525: the HaloTag7–JF525 model in this tab")
