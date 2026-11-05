"""Molecules from: M. Uji, N. Yanai, "Photon upconversion based on triplet-triplet annihilation using
thermally activated delayed fluorescence sensitizers", J. Photochem. Photobiol. C 66 (2026) 100741,
doi:10.1016/j.jphotochemrev.2026.100741.

Fields as in catalog.py, plus:
    ref      - reference number in the review and the DOI of the original paper
    formula  - (optional) molecular formula the build checks the structure against. Independent for
               source "pubchem" (PubChem record) and "si" (HRMS in the SI); for "drawing" it is a hand atom
               count of the figure, so it only confirms the SMILES encodes what was read from the drawing.
    source   - where the structure came from:
               "pubchem"  - SMILES from the PubChem record `cid`
               "si"       - structure + formula confirmed from the original paper's SI
               "drawing"  - transcribed from the structure drawing in the review (figure given in `fig`)

Data in descriptions (ηUC, Ith, ΔEST, ΦISC) are from Table 1 / the text of the review.
"""

CZ = "n{a}c{b}ccccc{b}c{b}ccccc{b}{a}"  # carbazol-9-yl (ring-closure digits filled in)


def cz(a, b):
    return CZ.format(a=a, b=b)


TBU = "C(C)(C)C"
TIPS = "[Si](C(C)C)(C(C)C)C(C)C"
# 9,9-dimethylacridan-10-yl; ring digits 5/6 so it can sit inside rings 1-4 that are still open
DMAC = "N5c6ccccc6C(C)(C)c6ccccc65"
DMAC_BR = "N5c6ccc(Br)cc6C(C)(C)c6cc(Br)ccc65"  # 2,7-dibromo (para to N)
# 2,7-di-tert-butylpyren-4-yl, written from C4 (ring digits 4-7; caller must keep 1-3 free)
DBP4 = "c4cc5cc(C(C)(C)C)cc6ccc7cc(C(C)(C)C)cc4c7c65"
# carbazol-9-yl with a 2-(2,7-di-tBu-pyren-4-yl)ethyl arm at C3 (and optionally C6); uses digits 2,3 + DBP4's
CZ_DBP1 = "n2c3ccc(CC" + DBP4 + ")cc3c3ccccc32"
CZ_DBP2 = "n2c3ccc(CC" + DBP4 + ")cc3c3cc(CC" + DBP4 + ")ccc32"
NAR2_TBU = "N(c5cc(C(C)(C)C)cc(C(C)(C)C)c5)c5cc(C(C)(C)C)cc(C(C)(C)C)c5"  # bis(3,5-di-tBu-phenyl)amino
TPA_VINYL = "/C=C/c2ccc(N(c3ccccc3)c3ccccc3)cc2"  # (E)-2-[4-(diphenylamino)phenyl]vinyl
PYRAN_ARM = "/C=C/C{r}=CC(=C(C#N)C#N)C=C(C)O{r}"  # (E)-2-(6-methyl-4-dicyanomethylene-4H-pyran-2-yl)vinyl

B_C = ("B", "C", 1.52)  # triarylborane B–C (DABNA-type X-ray structures)
UFF_NOTE = "MMFF94 has no boron parameters, so UFF is used with B–C restrained to the X-ray mean."

SENSITIZERS = [
    dict(file="4czipn.xyz", name="4CzIPN", source="pubchem", cid=102198498, ref="58, 60, 61, 64, 66",
         smiles="C1=CC=C2C(=C1)C3=CC=CC=C3N2C4=C(C(=C(C(=C4C#N)N5C6=CC=CC=C6C7=CC=CC=C75)"
                "N8C9=CC=CC=C9C1=CC=CC=C18)N1C2=CC=CC=C2C2=CC=CC=C21)C#N", formula="C56H32N6",
         desc="D–A TADF benchmark (ΔEST 0.083 eV, ΦISC ~70%). Vis→UV: p-terphenyl ηUC 2.8%, "
              "p-TIPS-BP 8%, bTIPS-Bz (UVB) 2%. Refs 58, 64, 66."),
    dict(file="4czpn.xyz", name="4CzPN", source="pubchem", cid=101136487, ref="57",
         smiles="N#Cc1c(C#N)c(" + cz(2, 3) + ")c(" + cz(2, 3) + ")c(" + cz(2, 3) + ")c1" + cz(2, 3),
         formula="C56H32N6",
         desc="3,4,5,6-Tetrakis(carbazol-9-yl)phthalonitrile. Parent of the pyrene-tethered conjugates. Ref 57."),
    dict(file="4czpn-1dbp.xyz", name="4CzPN-1DBP", source="si", ref="57 (doi:10.1039/c5sc03245h)",
         smiles="N#Cc1c(C#N)c(" + cz(2, 3) + ")c(" + CZ_DBP1 + ")c(" + cz(2, 3) + ")c1" + cz(2, 3),
         formula="C82H60N6",
         note="Arm = –CH2CH2– from carbazole C3 to pyrene C4 (SI synthesis); the arm-bearing carbazole is at "
              "C4 of the phthalonitrile (fluorine position in the 3CzFPN X-ray structure, SI CIF).",
         desc="4CzPN covalently linked to 2,7-di-tBu-pyrene: intramolecular TET gives ~80× brighter UC in "
              "polyurethane films. Ref 57."),
    dict(file="4czpn-2dbp.xyz", name="4CzPN-2DBP", source="si", ref="57 (doi:10.1039/c5sc03245h)",
         smiles="N#Cc1c(C#N)c(" + cz(2, 3) + ")c(" + CZ_DBP2 + ")c(" + cz(2, 3) + ")c1" + cz(2, 3),
         formula="C108H88N6",
         desc="Two pyrene arms on one carbazole (3,6-positions). Ref 57."),
    dict(file="4czpn-8dbp.xyz", name="4CzPN-8DBP", source="si", ref="57 (doi:10.1039/c5sc03245h)",
         smiles="N#Cc1c(C#N)c(" + CZ_DBP2 + ")c(" + CZ_DBP2 + ")c(" + CZ_DBP2 + ")c1" + CZ_DBP2,
         formula="C264H256N6", confs=4,
         desc="Eight pyrene arms (all carbazoles 3,6-substituted); 520 atoms. Self-quenching between pyrenes "
              "limits UC. Ref 57."),
    dict(file="4czbn.xyz", name="4CzBN", source="drawing", fig="7b", ref="63 (doi:10.1021/jacs.1c13222)",
         smiles="N#Cc1c(" + cz(2, 3) + ")c(" + cz(2, 3) + ")cc(" + cz(2, 3) + ")c1" + cz(2, 3),
         formula="C55H33N5",
         desc="2,3,5,6-Tetra(carbazol-9-yl)benzonitrile. T1 2.71 eV, τTADF 62 µs. With 1,4-TIPS-Naph: "
              "ηUC 26.2%, a record for vis→UV. Ref 63."),
    dict(file="4cztpn-ph.xyz", name="4CzTPN-Ph", source="drawing", fig="5a", ref="45 (doi:10.1063/1.4926914)",
         smiles="N#Cc1c(n2c3ccc(-c4ccccc4)cc3c3cc(-c4ccccc4)ccc32)c(n2c3ccc(-c4ccccc4)cc3c3cc(-c4ccccc4)ccc32)"
                "c(C#N)c(n2c3ccc(-c4ccccc4)cc3c3cc(-c4ccccc4)ccc32)c1n1c2ccc(-c3ccccc3)cc2c2cc(-c3ccccc3)ccc21",
         formula="C104H64N6",
         desc="2,3,5,6-Tetrakis(3,6-diphenylcarbazol-9-yl)-1,4-dicyanobenzene. First TADF-sensitized TTA-UC "
              "(Baldo 2015): 532 nm → DPA in a solid bilayer, ηUC 0.28%. Ref 45."),
    dict(file="bn-2cz.xyz", name="BN-2Cz", source="drawing", fig="8a", ref="44 (doi:10.31635/ccschem.022.202101507)",
         smiles="B12c3cccc4c5ccccc5n(c34)c3c1c(n1c4ccccc4c4cccc2c41)c(" + cz(1, 4) + ")cc3" + cz(1, 4),
         formula="C54H31BN4", restraints=[B_C], note=UFF_NOTE,
         desc="MR-TADF, ε 2.4×10⁵ M⁻¹cm⁻¹ at 490 nm, ΦISC 75%. Green-to-UV with 1,4-TIPS-Naph: ηUC 7.5%, "
              "Ith 15.8 mW cm⁻². Ref 44."),
    dict(file="bn-2cz-tbu.xyz", name="BN-2Cz-tBu", source="drawing", fig="8a",
         ref="44 (doi:10.31635/ccschem.022.202101507)",
         smiles="B12c3cc(" + TBU + ")cc4c5cc(" + TBU + ")ccc5n(c34)c3c1c(n1c4ccc(" + TBU + ")cc4c4cc(" + TBU
                + ")cc2c41)c(n1c4ccc(" + TBU + ")cc4c4cc(" + TBU + ")ccc41)cc3n1c4ccc(" + TBU + ")cc4c4cc("
                + TBU + ")ccc41",
         formula="C86H95BN4", restraints=[B_C], note=UFF_NOTE,
         desc="All carbazoles 3,6-di-tBu. Absorbs at 504 nm (ε 2.95×10⁵); Ith 9.2 mW cm⁻² with 1,4-TIPS-Naph; "
              "drives photopolymerization. Ref 44."),
    dict(file="bn-se.xyz", name="BN-Se", source="drawing", fig="9a", ref="62 (doi:10.1021/jacs.2c09543)",
         smiles="CC(C)(C)c1cc2c3cc(C(C)(C)C)ccc3n3c2c(c1)B1c2cccc4c2N(c2ccccc2[Se]4)c2cccc3c21",
         formula="C38H33BN2Se", restraints=[B_C, ("C", "Se", 1.92)],
         note="UFF (no MMFF parameters for B or Se); B–C and C–Se restrained to X-ray means.",
         desc="Se-embedded N/B/N MR-TADF: ΦISC 88%, ΔEST 0.08 eV, kRISC 1.6×10⁶ s⁻¹. ηUC 21.4% with "
              "1,4-TIPS-Naph; Ith 1.3 mW cm⁻². Ref 62."),
    dict(file="bbcz-sb-br.xyz", name="BBCz-SB-Br", source="drawing", fig="10a", ref="68 (doi:10.1039/d4tc01820f)",
         smiles="CC(C)(C)c1cc2c3cc(C(C)(C)C)ccc3n3c2c(c1)B1c2cc(C(C)(C)C)cc4c5cc(C(C)(C)C)ccc5n(c24)"
                "c2cc(Br)cc3c21",
         formula="C46H48BBrN2", restraints=[B_C], note=UFF_NOTE,
         desc="Br-substituted B/N MR-TADF, ΦISC ~100%, ΔEST 0.17 eV. All-solid-state vis→UV with DPN "
              "emitters, ηUC up to 16.2% in solution. Ref 68."),
    dict(file="dabna-1.xyz", name="DABNA-1", source="pubchem", cid=123168845, ref="85 (doi:10.1002/adma.201505491)",
         smiles="B12C3=CC=CC=C3N(C4=C1C(=CC=C4)N(C5=CC=CC=C25)C6=CC=CC=C6)C7=CC=CC=C7",
         formula="C30H21BN2", restraints=[B_C], note=UFF_NOTE,
         desc="The first multiple-resonance TADF emitter (Hatakeyama 2016). Ref 85."),
    dict(file="dabna-2.xyz", name="DABNA-2", source="drawing", fig="3c", ref="85 (doi:10.1002/adma.201505491)",
         smiles="B12c3ccc(-c4ccccc4)cc3N(c3cccc(-c4ccccc4)c3)c3cccc(c31)N(c1ccccc1)c1cc(N(c3ccccc3)c3ccccc3)ccc12",
         formula="C54H38BN3", restraints=[B_C], note=UFF_NOTE,
         desc="DABNA-1 with Ph and NPh2 para to boron: f = 0.415 vs 0.205 for DABNA-1. Ref 85."),
    dict(file="mr-1.xyz", name="MR-1", source="drawing", fig="4", ref="50 (doi:10.1002/adom.202300504)",
         smiles="O=C1c2cc(" + NAR2_TBU + ")cc3c2N2c4c1cc(" + NAR2_TBU + ")cc4C(=O)c1cc(" + NAR2_TBU
                + ")cc(c12)C3=O",
         formula="C105H132N4O3", confs=6,
         desc="N-centred triketone MR-TADF core with three bis(3,5-di-tBu-phenyl)amino donors. "
              "532 nm → DPA, ηUC 14.0%. Ref 50."),
    dict(file="btz-dmac.xyz", name="BTZ-DMAC", source="drawing", fig="4 (parent of BTZ-DMAC-4Br)",
         ref="46 (doi:10.1039/c7tc04096b)",
         smiles="CC1(C)c2ccccc2N(c2ccc(" + DMAC + ")c3nsnc23)c2ccccc21",
         formula="C36H30N4S",
         desc="Red D–A–D TADF (benzothiadiazole acceptor, 9,9-dimethylacridan donors). 532 nm → DPA, "
              "ηUC 1.9%. Ref 46."),
    dict(file="btz-dmac-4br.xyz", name="BTZ-DMAC-4Br", source="drawing", fig="4", ref="48 (doi:10.1002/adfm.202104044)",
         smiles="CC1(C)c2cc(Br)ccc2N(c2ccc(" + DMAC_BR + ")c3nsnc23)c2ccc(Br)cc21",
         formula="C36H26Br4N4S",
         desc="Tetrabrominated BTZ-DMAC: heavy-atom effect lifts ISC (singlet-O2 yield 70%); 520 nm → DPA, "
              "ηUC 24.6%. Ref 48."),
    dict(file="dmacpdo.xyz", name="DMACPDO", source="drawing", fig="4", ref="51 (doi:10.1021/acs.jpcc.2c07906)",
         smiles="O=C(/C=C(\\O)c1ccc(" + DMAC + ")cc1)c1ccc(" + DMAC + ")cc1",
         formula="C45H38N2O2",
         note="Drawn, as in the review, as the H-bonded (Z)-enol tautomer of the 1,3-diketone.",
         desc="β-Diketone D–A–D TADF (ΔEST 0.11 eV, ΦISC 86%). 455 nm → DPA, ηUC 22.3%. Ref 51."),
    dict(file="dcf-mpym.xyz", name="DCF-MPYM", source="drawing", fig="14a", ref="47 (doi:10.1039/c9cc01868a)",
         smiles="OC(=O)c1ccccc1C1=C2C=C(Cl)C(=O)C(" + PYRAN_ARM.format(r=3) + ")=C2Oc2c1cc(Cl)c(O)c2"
                + PYRAN_ARM.format(r=1),
         formula="C42H22Cl2N4O7",
         desc="Dichlorofluorescein with two pyranylidene-malononitrile arms: red TADF (ΔEST 0.028 eV). "
              "635 nm → DPA 428 nm, ηUC 11.2%. Ref 47."),
    dict(file="s1-curcuminoid.xyz", name="S1 (BF2 curcuminoid)", source="drawing", fig="16a",
         ref="53 (doi:10.1038/s41566-017-0087-y)",
         smiles="F[B-]1(F)OC(" + TPA_VINYL + ")=C(C(=O)OCC)C(" + TPA_VINYL + ")=[O+]1",
         formula="C46H37BF2N2O4", restraints=[("B", "F", 1.37), ("B", "O", 1.49)],
         note="UFF (no MMFF boron parameters); B–F / B–O restrained to BF2-β-diketonate X-ray means.",
         desc="NIR TADF BF2 complex (Adachi), ΦISC only 10.2%; with TTBP ηUC 0.052%. Ref 53."),
    dict(file="s2-curcuminoid.xyz", name="S2 (Br-BF2 curcuminoid)", source="drawing", fig="16a",
         ref="56 (doi:10.1002/adma.202405509)",
         smiles="F[B-]1(F)OC(" + TPA_VINYL + ")=C(Br)C(" + TPA_VINYL + ")=[O+]1",
         formula="C43H32BBrF2N2O2", restraints=[("B", "F", 1.37), ("B", "O", 1.49)],
         note="UFF (no MMFF boron parameters); B–F / B–O restrained to BF2-β-diketonate X-ray means.",
         desc="Br raises ΦISC to 53.7%. 724 nm NIR → TTBP 485 nm (ηUC 4.9%; 13.5% in methyl oleate); "
              "used for in vivo optogenetics. Ref 56."),
    dict(file="hep.xyz", name="Hep (heptazine)", source="drawing", fig="13a", ref="69 (doi:10.1002/adom.202501279)",
         smiles="COc1ccc(-c2nc3nc(N4CCCCC4)nc4nc(N5CCCCC5)nc(n2)n34)cc1",
         formula="C23H27N9O",
         desc="Heptazine TADF sensitizer, T1 2.96 eV, ΦISC 83%. UVA/visible → UVC (~295 nm) with "
              "TIPS-benzenes. Ref 69."),
]

EMITTERS = [
    dict(file="dpa.xyz", name="DPA", source="pubchem", cid=15159, formula="C26H18",
         smiles="C1=CC=C(C=C1)C2=C3C=CC=CC3=C(C4=CC=CC=C42)C5=CC=CC=C5", ref="45–51",
         desc="9,10-Diphenylanthracene. Blue emitter used with 4CzTPN-Ph, MR-1, BTZ-DMAC(-4Br), DMACPDO, "
              "DCF-MPYM."),
    dict(file="perylene.xyz", name="Perylene", source="pubchem", cid=9142, formula="C20H12",
         smiles="C1=CC2=C3C(=C1)C4=CC=CC5=C4C(=CC=C5)C3=CC=C2", ref="47, 54",
         desc="Blue emitter with DCF-MPYM (635 nm excitation) and BNS (721 nm, ηUC 3.0%)."),
    dict(file="ttbp.xyz", name="TTBP", source="pubchem", cid=15082495, formula="C36H44",
         smiles="CC(C)(C)C1=CC2=C3C(=C1)C=C(C=C3C4=CC(=CC5=CC(=CC2=C54)C(C)(C)C)C(C)(C)C)C(C)(C)C", ref="53, 56",
         desc="2,5,8,11-Tetra-tert-butylperylene, T1 1.53 eV, ΦFL 94.8%. NIR→blue with S2. Ref 56."),
    dict(file="tips-bea.xyz", name="TIPS-BEA", source="pubchem", cid=71310996, formula="C36H50Si2",
         smiles="CC(C)[Si](C#CC1=C2C=CC=CC2=C(C3=CC=CC=C31)C#C[Si](C(C)C)(C(C)C)C(C)C)(C(C)C)C(C)C", ref="54, 55",
         desc="9,10-Bis[(triisopropylsilyl)ethynyl]anthracene. With BNS: 721 nm → 451 nm, ηUC 5.8%."),
    dict(file="dbp.xyz", name="DBP (2,7-di-tBu-pyrene)", source="pubchem", cid=185563, formula="C24H26",
         smiles="CC(C)(C)C1=CC2=C3C(=C1)C=CC4=CC(=CC(=C43)C=C2)C(C)(C)C", ref="57",
         desc="T1 2.5 eV, S1 3.3 eV. The emitter tethered to 4CzPN in the DBP conjugates. Ref 57."),
    dict(file="14-tips-naph.xyz", name="1,4-TIPS-Naph", source="drawing", fig="8a", formula="C32H48Si2",
         smiles="CC(C)[Si](C#Cc1ccc(C#C" + TIPS + ")c2ccccc12)(C(C)C)C(C)C",
         ref="44, 62, 63",
         desc="1,4-Bis(TIPS-ethynyl)naphthalene, T1 2.11 eV. Best UV emitter here: ηUC 26.2% with 4CzBN, 21.4% "
              "with BN-Se."),
    dict(file="15-tips-naph.xyz", name="1,5-TIPS-Naph", source="drawing", fig="8a", formula="C32H48Si2",
         smiles="CC(C)[Si](C#Cc1cccc2c(C#C" + TIPS + ")cccc12)(C(C)C)C(C)C", ref="44",
         desc="1,5-Bis(TIPS-ethynyl)naphthalene, T1 2.24 eV. With BN-2Cz-tBu (ηUC 2.5%). Ref 44."),
    dict(file="1-tna.xyz", name="1-TNA", source="drawing", fig="9a", formula="C21H28Si",
         smiles="CC(C)[Si](C#Cc1cccc2ccccc12)(C(C)C)C(C)C", ref="62",
         desc="1-(TIPS-ethynyl)naphthalene, T1 2.40 eV. With BN-Se: 0.93 eV anti-Stokes shift. Ref 62."),
    dict(file="p-terphenyl.xyz", name="p-Terphenyl (TP)", source="pubchem", cid=7115, formula="C18H14",
         smiles="C1=CC=C(C=C1)C2=CC=C(C=C2)C3=CC=CC=C3", ref="58, 63",
         desc="UV emitter (343 nm) with 4CzIPN: ηUC 2.8%, ΔE 0.83 eV. Ref 58."),
    dict(file="p-quaterphenyl.xyz", name="p-Quaterphenyl", source="pubchem", cid=8677, formula="C24H18",
         smiles="C1=CC=C(C=C1)C2=CC=C(C=C2)C3=CC=C(C=C3)C4=CC=CC=C4", ref="58",
         desc="Emits at 353 nm with 4CzIPN (ηUC 3.9%). Ref 58."),
    dict(file="r-s-tp.xyz", name="R/S-TP (chiral terphenyl dimer)", source="si", formula="C40H32",
         ref="61 (doi:10.1038/s41467-020-19479-1)",
         smiles="C1Cc2c(-c3ccc(-c4ccccc4)cc3)cc(cc2)CCc2c(-c3ccc(-c4ccccc4)cc3)cc1cc2",
         note="4,12-Bis(biphenyl-4-yl)[2.2]paracyclophane per the SI synthesis (MALDI m/z 512, C40H32). "
              "Planar-chiral: the 3D model is one enantiomer, not assigned as R or S.",
         desc="Two p-terphenyls stacked on a [2.2]paracyclophane. With 4CzIPN: circularly polarized UV "
              "upconversion at 385 nm, ηUC 7.9%. Ref 61."),
    dict(file="dpn.xyz", name="DPN", source="pubchem", cid=98868, formula="C22H16",
         smiles="C1=CC=C(C=C1)C2=CC=C(C3=CC=CC=C32)C4=CC=CC=C4", ref="68",
         desc="1,4-Diphenylnaphthalene. Emitter and host for all-solid-state UC with BBCz-SB-Br. Ref 68."),
    dict(file="sidpn-1.xyz", name="SiDPN-1", source="drawing", fig="10a", formula="C40H30Si",
         smiles="c1ccc([Si](c2ccccc2)(c2ccccc2)c2ccc(-c3ccc(-c4ccccc4)c4ccccc34)cc2)cc1", ref="68",
         desc="DPN with a triphenylsilyl group: amorphous glassy film, UC lifetime 2.8 ms. Ref 68."),
    dict(file="sidpn-2.xyz", name="SiDPN-2", source="drawing", fig="10a", formula="C56H40Si",
         smiles="c1ccc([Si](c2ccccc2)(c2ccc(-c3ccc(-c4ccccc4)c4ccccc34)cc2)c2ccc(-c3ccc(-c4ccccc4)c4ccccc34)cc2)cc1",
         ref="68", desc="Two DPN units on one silicon; Ith 44 mW cm⁻² in the solid state. Ref 68."),
    dict(file="p-tips-bp.xyz", name="p-TIPS-BP", source="drawing", fig="11a", formula="C23H30Si",
         smiles="CC(C)[Si](C#Cc1ccc(-c2ccccc2)cc1)(C(C)C)C(C)C", ref="64 (doi:10.1039/d1tc04782e)",
         desc="4-(TIPS-ethynyl)biphenyl, S1 4.04 eV. 447 nm → 328 nm with 4CzIPN, ΔE 1.01 eV. Ref 64."),
    dict(file="btips-bz.xyz", name="bTIPS-Bz", source="drawing", fig="12a", formula="C28H46Si2",
         smiles="CC(C)[Si](C#Cc1ccc(C#C" + TIPS + ")cc1)(C(C)C)C(C)C", ref="66 (doi:10.1002/anie.202215340)",
         desc="1,4-Bis(TIPS-ethynyl)benzene, S1 4.15 eV. First purely organic visible→UVB UC (318 nm); "
              "drives a Norrish type I reaction. Ref 66."),
    dict(file="tips-bz.xyz", name="TIPS-Bz", source="drawing", fig="13a", formula="C17H26Si",
         smiles="CC(C)[Si](C#Cc1ccccc1)(C(C)C)C(C)C", ref="69",
         desc="(Phenylethynyl)triisopropylsilane, S1 4.43 eV. UVC emitter (~295 nm) with Hep. Ref 69."),
    dict(file="bu-tips-bz.xyz", name="Bu-TIPS-Bz", source="drawing", fig="13a", formula="C21H34Si",
         smiles="CC(C)[Si](C#Cc1ccc(C(C)(C)C)cc1)(C(C)C)C(C)C", ref="69",
         desc="4-tert-Butyl TIPS-Bz; triplet lifetime 14.5 µs. Ref 69."),
    dict(file="bbu-tips-bz.xyz", name="bBu-TIPS-Bz", source="drawing", fig="13a", formula="C25H42Si",
         smiles="CC(C)[Si](C#Cc1cc(C(C)(C)C)cc(C(C)(C)C)c1)(C(C)C)C(C)C", ref="69",
         desc="3,5-Di-tert-butyl TIPS-Bz: bulk suppresses triplet excimers (τT 29.5 µs); best of the series. Ref 69."),
    dict(file="ppd.xyz", name="PPD", source="pubchem", cid=12888, formula="C14H10N2O",
         smiles="C1=CC=C(C=C1)C2=NN=C(O2)C3=CC=CC=C3", ref="63",
         desc="2,5-Diphenyl-1,3,4-oxadiazole; UV emitter screened with 4CzBN. Ref 63."),
    dict(file="ppo.xyz", name="PPO", source="pubchem", cid=7105, formula="C15H11NO",
         smiles="C1=CC=C(C=C1)C2=CN=C(O2)C3=CC=CC=C3", ref="63",
         desc="2,5-Diphenyloxazole (scintillator); ηUC ~14% with 4CzBN. Ref 63."),
    dict(file="ppf.xyz", name="PPF", source="pubchem", cid=70387, formula="C16H12O",
         smiles="C1=CC=C(C=C1)C2=CC=C(O2)C3=CC=CC=C3", ref="63",
         desc="2,5-Diphenylfuran; UV emitter screened with 4CzBN. Ref 63."),
    dict(file="2pi.xyz", name="2PI", source="pubchem", cid=353435, formula="C15H12",
         smiles="C1C2=CC=CC=C2C=C1C3=CC=CC=C3", ref="63",
         desc="2-Phenylindene; UV emitter screened with 4CzBN. Ref 63."),
]
