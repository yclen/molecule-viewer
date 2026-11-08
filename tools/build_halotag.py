"""HaloTag7 with JF525 covalently attached, modelled on the crystal structure of the TMR conjugate (PDB 6Y7A).

No crystal structure of HaloTag with JF525 exists, but JF525 and the TMR ligand in 6Y7A differ only in the two
amino groups: TMR has N(CH3)2, JF525 has 3,3-difluoroazetidine. Same xanthene, same pendant benzoate, same linker
at the same position. So:

  1. every protein atom and every shared dye atom keeps its crystal position, including the ester bond from the
     linker's last carbon to Asp106 (the chloride left when that bond formed, so the conjugate has no Cl);
  2. RDKit grows the azetidine rings and all dye hydrogens onto those fixed atoms (MMFF94s, the fixed atoms held);
  3. the seed with the fewest close contacts to the protein wins, and the contacts are reported.

Called by tools/build_examples.py for catalog entries with `conjugate`.
"""

import math

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from rdkit.Geometry import Point3D

# JF525 HaloTag ligand with the chloride replaced by an acetate, which stands in for the end of the Asp106 side chain
# (CB-CG(=OD1)-OD2). Open zwitterion: the fluorescent form, and the one bound to HaloTag.
ESTER_SMILES = ("[O-]C(=O)c1ccc(cc1C2=C3C=CC(=[N+]4CC(F)(F)C4)C=C3Oc5cc(N6CC(F)(F)C6)ccc25)"
                "C(=O)NCCOCCOCCCCCCOC(=O)C")
SEEDS = range(1, 41)
BOND_MAX = 1.9  # Å, heavy atoms in the ligand record
CONTACT = 2.8  # Å, new dye heavy atom to protein heavy atom closer than this counts as a close contact


def read_pdb(text, chain, ligand):
    """Protein atoms and ligand atoms of one chain (first model, first alternate location)."""
    protein, lig, link = [], {}, None
    for line in text.splitlines():
        rec = line[:6].strip()
        if rec == "ENDMDL":
            break
        if rec == "LINK" and ligand in line:
            link = line
        if rec not in ("ATOM", "HETATM") or line[21] != chain or line[16] not in (" ", "A"):
            continue
        name, resname, resseq = line[12:16].strip(), line[17:20].strip(), int(line[22:26])
        element = (line[76:78].strip() or name[0]).capitalize()
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        if rec == "ATOM":
            protein.append(dict(name=name, res=resname, seq=resseq, element=element, xyz=xyz))
        elif resname == ligand:
            lig[name] = dict(element=element, xyz=xyz)
    return protein, lig, link


def skeleton_from_coords(atoms):
    """Heavy-atom graph (all single bonds) of [(element, xyz)], bonds found from distances."""
    m = Chem.RWMol()
    for element, _ in atoms:
        a = Chem.Atom(element)
        a.SetNoImplicit(True)
        m.AddAtom(a)
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            if math.dist(atoms[i][1], atoms[j][1]) < BOND_MAX:
                m.AddBond(i, j, Chem.BondType.SINGLE)
    m.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(m)
    return m


def plain_skeleton(mol):
    """The same graph with bond orders, aromaticity and charges erased, keeping atom order (heavy atoms first)."""
    m = Chem.RWMol(mol)
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


def build_conjugate(pdb_text, entry):
    chain, ligand = entry.get("chain", "A"), entry["conjugate"]["ligand"]
    protein, lig, link = read_pdb(pdb_text, chain, ligand)
    if not link or "ASP A 106" not in link:
        raise RuntimeError("expected a LINK record from Asp106 to the ligand")
    asp = {a["name"]: a for a in protein if a["seq"] == 106}

    # crystal atoms the model keeps: the ligand, plus the end of Asp106 so the ester comes out right
    crystal = [(a["element"], a["xyz"]) for a in lig.values()]
    crystal += [(asp[n]["element"], asp[n]["xyz"]) for n in ("OD2", "CG", "OD1", "CB")]
    query = skeleton_from_coords(crystal)

    target = Chem.MolFromSmiles(ESTER_SMILES)
    match = plain_skeleton(target).GetSubstructMatch(query)
    if len(match) != len(crystal):
        raise RuntimeError("the crystal ligand does not fit the JF525 skeleton")

    # the N-methyl carbons become ring CH2 groups: they have to move, since a ring pinches the C-N-C angle
    ring_c = {i for i in match if target.GetAtomWithIdx(i).GetSymbol() == "C"
              and target.GetAtomWithIdx(i).IsInRingSize(4)}
    fixed = {t: crystal[q][1] for q, t in enumerate(match) if t not in ring_c}

    mol = Chem.AddHs(target)
    protein_heavy = np.array([a["xyz"] for a in protein if a["element"] != "H"])
    new_atoms = [a.GetIdx() for a in mol.GetAtoms() if a.GetIdx() not in fixed]

    best = None
    for seed in SEEDS:
        trial = Chem.Mol(mol)
        coord_map = {i: Point3D(*p) for i, p in fixed.items()}
        if AllChem.EmbedMolecule(trial, coordMap=coord_map, randomSeed=seed, useRandomCoords=True) < 0:
            continue
        ref = Chem.Mol(trial)
        conf = ref.GetConformer()
        for i, p in fixed.items():
            conf.SetAtomPosition(i, Point3D(*p))
        rdMolAlign.AlignMol(trial, ref, atomMap=[(i, i) for i in fixed])
        conf = trial.GetConformer()
        for i, p in fixed.items():
            conf.SetAtomPosition(i, Point3D(*p))
        props = AllChem.MMFFGetMoleculeProperties(trial, mmffVariant="MMFF94s")
        field = AllChem.MMFFGetMoleculeForceField(trial, props)
        for i in fixed:
            field.AddFixedPoint(i)
        field.Initialize()
        for _ in range(20):
            if field.Minimize(maxIts=5000) == 0:
                break
        pos = trial.GetConformer().GetPositions()
        heavy_new = [i for i in new_atoms if trial.GetAtomWithIdx(i).GetAtomicNum() > 1]
        d = np.linalg.norm(pos[heavy_new, None, :] - protein_heavy[None, :, :], axis=2)
        score = (int((d < CONTACT).sum()), field.CalcEnergy())
        if best is None or score < best[0]:
            best = (score, trial, float(d.min()))
    if best is None:
        raise RuntimeError("could not embed the dye on the crystal atoms")
    (n_close, energy), trial, closest = best

    heavy_pos = trial.GetConformer().GetPositions()
    protein_atoms = [a for a in protein if a["element"] != "H"]
    contacts = []
    for i in new_atoms:
        atom = trial.GetAtomWithIdx(i)
        if atom.GetAtomicNum() == 1:
            continue
        for a, p in zip(protein_atoms, protein_heavy):
            d = float(np.linalg.norm(heavy_pos[i] - p))
            if d < CONTACT:
                contacts.append(f"{atom.GetSymbol()} to {a['res']}{a['seq']} {a['name']} {d:.2f} Å")

    # drop the acetate stand-in (and its H): those atoms are Asp106, which the protein already has
    asp_atoms = {match[len(lig) + k] for k in range(4)}
    asp_atoms |= {n.GetIdx() for i in asp_atoms for n in trial.GetAtomWithIdx(i).GetNeighbors()
                  if n.GetAtomicNum() == 1}
    keep = [a.GetIdx() for a in trial.GetAtoms() if a.GetIdx() not in asp_atoms]
    pos = trial.GetConformer().GetPositions()

    symbols = [a["element"] for a in protein]
    coords = np.array([a["xyz"] for a in protein], dtype=float)
    dye_symbols = [trial.GetAtomWithIdx(i).GetSymbol() for i in keep]
    dye_coords = pos[keep]
    all_coords = np.vstack([coords, dye_coords])

    n_protein = len(symbols)
    residues = len({(a["seq"], a["res"]) for a in protein})
    method = (f"RCSB PDB {entry['pdb']} (HaloTag7 with TMR, X-ray 1.40 Å), chain {chain}. Crystal positions kept for "
              f"all protein atoms and the dye atoms shared with TMR; difluoroazetidines and dye H added by RDKit "
              f"MMFF94s with those atoms fixed")
    return dict(symbols=symbols + dye_symbols, coords=[tuple(p) for p in all_coords], method=method,
                n_res=residues, dye=list(range(n_protein, n_protein + len(dye_symbols))),
                contacts=contacts, closest=closest)
