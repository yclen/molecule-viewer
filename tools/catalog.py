"""Catalog of example molecules, grouped into the tabs shown in the app.

Each entry:
    file     - output name (inside samples/<tab>/)
    name     - display name
    desc     - one-line description (spectral values are approximate, solvent-dependent)
    smiles   - structure. Copied from the PubChem record `cid` unless `note` says otherwise.
    cid      - PubChem CID the structure was taken from / checked against
    note     - (optional) why the SMILES differs from the PubChem record
    restraints - (optional) [(el1, el2, Å)] bonds pulled to crystallographic means after minimisation
    post     - (optional) post-processing step, e.g. "porphyrin:Pt"

Counter-ions (Na+, K+, Cl-, Br-) are dropped: the coloured/fluorescent species is drawn.
"""

from catalog_uji2026 import EMITTERS as _UC_EMITTERS
from catalog_uji2026 import SENSITIZERS as _UC_SENSITIZERS

TABS = [
    ("basics", "Basics"),
    ("triplet", "Triplet molecules"),
    ("xanthenes", "Xanthenes"),
    ("cyanines", "Cyanines"),
    ("dyes", "Other dyes"),
    ("proteins", "Proteins"),
    ("tta", "TTA upconversion"),
    ("tadf", "TADF"),
    ("uc-sens", "Review: TADF sensitizers"),
    ("uc-emit", "Review: UC emitters"),
]

MOLECULES = {
    # ------------------------------------------------------------------ Triplet project (parametric objects)
    # `param` entries are built in the app by stacking a repeat unit made by tools/build_wire.py
    # (samples/triplet/<data>); build_examples.py only lists them.
    "triplet": [
        dict(file="molecular-wire", name="Molecular wire",
             desc="n DPA units joined anthracene C2 → 1,4-phenylene → C6, H on both ends (R = H). "
                  "Set n below; n = 1 is plain DPA.",
             param=dict(data="dpa-wire.param.json", label="DPA units", symbol="n", min=1, max=20, default=3)),
        dict(file="tta-flora-4", name="TTA Flora 4",
             desc="Sensitizer (bis-triarylamine thienoquinoxaline dioxide) – phenylene – n DPA mediators – "
                  "perylene. n = 0 is the direct sensitizer–linker–perylene dyad.",
             param=dict(data="tta-flora-4.param.json", label="DPA mediators", symbol="n", min=0, max=20,
                        default=2)),
    ],
    # ------------------------------------------------------------------ xanthene dyes
    "xanthenes": [
        dict(file="fluorescein.xyz", name="Fluorescein (dianion)",
             desc="Abs ~490 nm / Em ~515 nm (pH > 8). The open dianion is the bright species.",
             smiles="[O-]C(=O)c1ccccc1C1=C2C=CC(=O)C=C2Oc2cc([O-])ccc21", cid=16850,
             note="PubChem 16850 is the closed (non-fluorescent) lactone; drawn here as the open dianion."),
        dict(file="eosin-y.xyz", name="Eosin Y",
             desc="2',4',5',7'-Tetrabromofluorescein. Abs ~517 / Em ~540 nm. Photoredox catalyst.",
             smiles="C1=CC=C(C(=C1)C2=C3C=C(C(=O)C(=C3OC4=C(C(=C(C=C24)Br)[O-])Br)Br)Br)C(=O)[O-]", cid=11048),
        dict(file="erythrosin-b.xyz", name="Erythrosin B",
             desc="Tetraiodofluorescein. Abs ~526 / Em ~550 nm. Heavy-atom effect: strong triplet yield.",
             smiles="C1=CC=C(C(=C1)C2=C3C=C(C(=O)C(=C3OC4=C(C(=C(C=C24)I)[O-])I)I)I)C(=O)[O-]", cid=3259, restraints=[("C", "I", 2.095)],
             note="PubChem 3259 is the lactone; drawn as the open dianion (eosin Y skeleton with I for Br)."),
        dict(file="rose-bengal.xyz", name="Rose Bengal",
             desc="Tetrachloro-tetraiodofluorescein. Abs ~549 / Em ~570 nm. Efficient singlet-oxygen sensitizer.",
             smiles="C1=C2C(=C3C=C(C(=O)C(=C3OC2=C(C(=C1I)[O-])I)I)I)C4=C(C(=C(C(=C4Cl)Cl)Cl)Cl)C(=O)[O-]", cid=32343,
             restraints=[("C", "I", 2.095)]),
        dict(file="rhodamine-b.xyz", name="Rhodamine B",
             desc="Abs ~554 / Em ~577 nm (EtOH). Classic laser dye and tracer.",
             smiles="CCN(CC)C1=CC2=C(C=C1)C(=C3C=CC(=[N+](CC)CC)C=C3O2)C4=CC=CC=C4C(=O)O", cid=6694),
        dict(file="rhodamine-6g.xyz", name="Rhodamine 6G",
             desc="Abs ~530 / Em ~552 nm. Quantum yield ~0.95; the standard laser dye.",
             smiles="CCNC1=CC2=C(C=C1C)C(=C3C=C(C(=[NH+]CC)C=C3O2)C)C4=CC=CC=C4C(=O)OCC", cid=13806),
        dict(file="tmr.xyz", name="Tetramethylrhodamine (TMR)",
             desc="Abs ~550 / Em ~575 nm. Common label for proteins and antibodies.",
             smiles="CN(C)C1=CC2=C(C=C1)C(=C3C=CC(=[N+](C)C)C=C3O2)C4=CC=CC=C4C(=O)[O-]", cid=9952143),
        dict(file="texas-red.xyz", name="Texas Red",
             desc="Sulforhodamine 101 sulfonyl chloride. Abs ~589 / Em ~615 nm.",
             smiles="C1CC2=CC3=C(C4=C2N(C1)CCC4)OC5=C6CCC[N+]7=C6C(=CC5=C3C8=C(C=C(C=C8)S(=O)(=O)Cl)S(=O)(=O)[O-])CCC7",
             cid=452705),
    ],
    # ------------------------------------------------------------------ cyanines
    "cyanines": [
        dict(file="cy3.xyz", name="Cy3 (sulfo, acid)",
             desc="Trimethine indocyanine. Abs ~550 / Em ~570 nm.",
             smiles="CCN\\1C2=C(C=C(C=C2)S(=O)(=O)[O-])C(/C1=C/C=C/C3=[N+](C4=C(C3(C)C)C=C(C=C4)S(=O)(=O)O)CCCCCC(=O)O)(C)C",
             cid=91811029),
        dict(file="cy5.xyz", name="Cy5 (sulfo, acid)",
             desc="Pentamethine indocyanine. Abs ~649 / Em ~670 nm. Workhorse red dye for single-molecule work.",
             smiles="CCN\\1C2=C(C=C(C=C2)S(=O)(=O)[O-])C(/C1=C/C=C/C=C/C3=[N+](C4=C(C3(C)C)C=C(C=C4)S(=O)(=O)O)CCCCCC(=O)O)(C)C",
             cid=17758493, note="Polymethine chain written all-trans (the PubChem record mixes / and \\ markers)."),
        dict(file="cy7.xyz", name="Cy7 (sulfo, acid)",
             desc="Heptamethine indocyanine. Abs ~750 / Em ~773 nm (near-IR).",
             smiles="CCN\\1C2=C(C=C(C=C2)S(=O)(=O)[O-])C(/C1=C/C=C/C=C/C=C/C3=[N+](C4=C(C3(C)C)C=C(C=C4)S(=O)(=O)O)CCCCCC(=O)O)(C)C",
             cid=73554281),
        dict(file="icg.xyz", name="Indocyanine green (ICG)",
             desc="FDA-approved NIR dye. Abs ~780 / Em ~810 nm. Used in angiography and surgery.",
             smiles="CC1(C(=[N+](C2=C1C3=CC=CC=C3C=C2)CCCCS(=O)(=O)[O-])/C=C/C=C/C=C/C=C\\4/C(C5=C(N4CCCCS(=O)(=O)[O-])C=CC6=CC=CC=C65)(C)C)C",
             cid=5282412),
    ],
    # ------------------------------------------------------------------ other small-molecule dyes
    "dyes": [
        dict(file="coumarin-153.xyz", name="Coumarin 153",
             desc="Abs ~423 / Em ~530 nm (EtOH). Classic solvation-dynamics probe.",
             smiles="C1CC2=CC3=C(C4=C2N(C1)CCC4)OC(=O)C=C3C(F)(F)F", cid=72652),
        dict(file="bodipy.xyz", name="BODIPY core",
             desc="4,4-Difluoro-4-bora-3a,4a-diaza-s-indacene. Abs ~497 / Em ~505 nm; very narrow bands.",
             smiles="[B-]1(N2C=CC=C2C=C3[N+]1=CC=C3)(F)F", cid=25058173,
             restraints=[("B", "F", 1.39), ("B", "N", 1.55)],
             note="MMFF94 has no boron parameters, so UFF is used with B-F / B-N pulled to BODIPY X-ray means."),
        dict(file="nile-red.xyz", name="Nile Red",
             desc="Solvatochromic lipid stain. Em ~550-650 nm depending on environment.",
             smiles="CCN(CC)C1=CC2=C(C=C1)N=C3C4=CC=CC=C4C(=O)C=C3O2", cid=65182),
        dict(file="prodan.xyz", name="Prodan",
             desc="Push-pull naphthalene; strongly solvatochromic (Em ~400-530 nm).",
             smiles="CCC(=O)C1=CC2=C(C=C1)C=C(C=C2)N(C)C", cid=107729),
        dict(file="acridine-orange.xyz", name="Acridine Orange",
             desc="Nucleic-acid stain. Em ~525 nm (dsDNA) / ~650 nm (RNA).",
             smiles="CN(C)C1=CC2=C(C=C1)C=C3C=CC(=CC3=N2)N(C)C", cid=62344),
        dict(file="dapi.xyz", name="DAPI",
             desc="Blue DNA stain (minor groove). Abs ~358 / Em ~461 nm when bound.",
             smiles="C1=CC(=CC=C1C2=CC3=C(N2)C=C(C=C3)C(=N)N)C(=N)N", cid=2954),
        dict(file="hoechst-33258.xyz", name="Hoechst 33258",
             desc="Cell-permeant blue DNA stain. Abs ~352 / Em ~461 nm.",
             smiles="CN1CCN(CC1)C2=CC3=C(C=C2)N=C(N3)C4=CC5=C(C=C4)N=C(N5)C6=CC=C(C=C6)O", cid=2392),
        dict(file="ethidium.xyz", name="Ethidium",
             desc="Intercalating DNA stain. Em ~605 nm, ~20x brighter when bound.",
             smiles="CC[N+]1=C2C=C(C=CC2=C3C=CC(=CC3=C1C4=CC=CC=C4)N)N", cid=14710),
        dict(file="quinine.xyz", name="Quinine",
             desc="Fluorescence quantum-yield standard (in 0.5 M H2SO4). Abs ~350 / Em ~450 nm.",
             smiles="COC1=CC2=C(C=CN=C2C=C1)[C@H]([C@@H]3C[C@@H]4CCN3C[C@@H]4C=C)O", cid=3034034),
    ],
    # ------------------------------------------------------------------ fluorescent proteins (PDB) + chromophore
    "proteins": [
        dict(file="hbdi.xyz", name="GFP chromophore (HBDI)",
             desc="p-Hydroxybenzylidene-imidazolinone, the Z-isomer that forms inside GFP.",
             smiles="CC1=N/C(=C\\c2ccc(O)cc2)C(=O)N1C", cid=129849628,
             note="PubChem 129849628 has the carbonyl on the wrong ring carbon; this is 4-(4-hydroxybenzylidene)-"
                  "1,2-dimethyl-imidazol-5-one, (Z)."),
        # X-ray crystal structures from the RCSB Protein Data Bank (chain A, waters removed)
        # keep_het: the chromophore residue; other HETATM groups (solvent, ions, additives) are dropped
        dict(file="gfp.xyz", name="GFP (wild-type, 2WUR)", pdb="2WUR", chain="A", keep_het=["GYS"],
             desc="Aequorea victoria GFP, 0.9 Å X-ray structure. Abs 395/475 / Em 509 nm."),
        dict(file="egfp.xyz", name="EGFP (2Y0G)", pdb="2Y0G", chain="A", keep_het=["CRO"],
             desc="Enhanced GFP (F64L/S65T). Abs 488 / Em 507 nm. The standard green FP."),
        dict(file="mcherry.xyz", name="mCherry (2H5Q)", pdb="2H5Q", chain="A", keep_het=["CH6"],
             desc="Monomeric red FP derived from DsRed. Abs 587 / Em 610 nm."),
    ],
    # ------------------------------------------------------------------ triplet-triplet annihilation upconversion
    "tta": [
        dict(file="dpa.xyz", name="9,10-Diphenylanthracene (DPA)",
             desc="Benchmark annihilator. Abs ~375 / Em ~410-430 nm, QY ~0.95.",
             smiles="C1=CC=C(C=C1)C2=C3C=CC=CC3=C(C4=CC=CC=C42)C5=CC=CC=C5", cid=15159),
        dict(file="perylene.xyz", name="Perylene",
             desc="Blue annihilator. Abs ~436 / Em ~445-475 nm.",
             smiles="C1=CC2=C3C(=C1)C4=CC=CC5=C4C(=CC=C5)C3=CC=C2", cid=9142),
        dict(file="pyrene.xyz", name="Pyrene",
             desc="UV annihilator; famous for excimer emission (~480 nm).",
             smiles="C1=CC2=C3C(=C1)C=CC4=CC=CC(=C43)C=C2", cid=31423),
        dict(file="bpea.xyz", name="BPEA",
             desc="9,10-Bis(phenylethynyl)anthracene. Green-emitting annihilator, Em ~470-490 nm.",
             smiles="C1=CC=C(C=C1)C#CC2=C3C=CC=CC3=C(C4=CC=CC=C42)C#CC5=CC=CC=C5", cid=82338),
        dict(file="rubrene.xyz", name="Rubrene",
             desc="Tetraphenyltetracene. Yellow annihilator (Em ~560 nm); also a singlet-fission material.",
             smiles="C1=CC=C(C=C1)C2=C3C=CC=CC3=C(C4=C(C5=CC=CC=C5C(=C24)C6=CC=CC=C6)C7=CC=CC=C7)C8=CC=CC=C8",
             cid=68203),
        dict(file="ptoep.xyz", name="PtOEP (sensitizer)",
             desc="Pt(II) octaethylporphyrin. Absorbs ~535 nm, phosphoresces ~646 nm; pairs with DPA.",
             smiles="CCC1=C(C2=CC3=NC(=CC4=NC(=CC5=C(C(=C([N-]5)C=C1[N-]2)CC)CC)C(=C4CC)CC)C(=C3CC)CC)CC",
             cid=6537178, post="porphyrin:Pt",
             note="Porphyrin dianion from PubChem 6537178. MMFF puckers porphyrin cores, so the 24-atom core is "
                  "placed from metalloporphyrin X-ray averages and only the ethyls/H are force-field relaxed."),
    ],
    # ------------------------------------------------------------------ thermally activated delayed fluorescence
    "tadf": [
        dict(file="4czipn.xyz", name="4CzIPN",
             desc="Carbazole donors around a dicyanobenzene acceptor. Em ~507-520 nm; also a photocatalyst.",
             smiles="C1=CC=C2C(=C1)C3=CC=CC=C3N2C4=C(C(=C(C(=C4C#N)N5C6=CC=CC=C6C7=CC=CC=C75)"
                    "N8C9=CC=CC=C9C1=CC=CC=C18)N1C2=CC=CC=C2C2=CC=CC=C21)C#N",
             cid=102198498),
        dict(file="dmac-trz.xyz", name="DMAC-TRZ",
             desc="Acridan donor + triazine acceptor, nearly perpendicular. Sky-blue TADF emitter (~495 nm).",
             smiles="CC1(C2=CC=CC=C2N(C3=CC=CC=C31)C4=CC=C(C=C4)C5=NC(=NC(=N5)C6=CC=CC=C6)C7=CC=CC=C7)C",
             cid=118528399),
    ],
    # ------------------------------------------------------------------ Uji & Yanai 2026 review (catalog_uji2026.py)
    "uc-sens": _UC_SENSITIZERS,
    "uc-emit": _UC_EMITTERS,
}
