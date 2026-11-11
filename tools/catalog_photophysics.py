"""The "Interesting photophysics" tab: systems that each teach one piece of photophysics, with an About-panel text.

Grouped as energy transfer (FRET, BRET), light harvesting (LH2, PCP, FMO), triplets put to work (LOV2, miniSOG,
KillerRed), and fluorogenicity / sensing (JF525 open vs closed, GCaMP6m). Protein entries use `cofactors`
(tools/build_cofactors.py): the protein drawn faint, each photoactive cofactor a site with its H restored. Numbers in
the texts are approximate literature values; distances marked "here" are measured on the model.
"""

from textwrap import dedent


def info(text):
    return dedent(text).strip()


CHL = dict(name="BChl", kind="chlorophyll", core="macrocycle")

PHOTOPHYSICS = [
    # ------------------------------------------------------------------ energy transfer
    dict(file="dna-fret-ruler", name="DNA FRET ruler: Cy3–Cy5",
         desc="Cy3 and Cy5 stacked on the two ends of a B-DNA duplex of n base pairs. FRET falls with distance, and "
              "wobbles as the helix turns the dyes against each other (κ²).",
         param=dict(data="dna-fret-ruler.param.json", label="Base pairs", symbol="n", min=4, max=30, default=10),
         info=info("""
            ### What it is
            Two dyes on the ends of a DNA duplex: **Cy3** (the donor, absorbs ~550 nm, emits ~570 nm) and **Cy5** (the
            acceptor, absorbs ~650 nm, emits ~670 nm). Adding base pairs moves them apart by ~3.3 Å per pair, which
            makes DNA the standard ruler for calibrating FRET.

            ### The photophysics
            In **Förster resonance energy transfer**, an excited donor hands its energy to a nearby acceptor through
            the coupling of their transition dipoles, without emitting a photon. The efficiency is
            E = 1 / (1 + (r / R₀)⁶): it falls from ~100% to ~0% over about a factor of two in distance around the
            Förster radius R₀ (~5–6 nm for Cy3–Cy5).

            R₀⁶ is proportional to **κ²**, the orientation factor, κ = cos θ_DA − 3 cos θ_D cos θ_A. It runs from 0
            (dipoles at right angles) to 4 (in line). For dyes tumbling freely it averages to 2/3, the value almost
            every FRET calculation assumes. But dyes stacked on the ends of a helix can't tumble: they turn with the
            DNA, ~36° per base pair, so κ² and the FRET efficiency oscillate as n changes, on top of the r⁶ decay.
            That oscillation was measured for terminally stacked Cy3 and Cy5 (Iqbal et al. 2008).

            ### Why it matters for imaging
            Single-molecule FRET (smFRET) with Cy3/Cy5 measures distances and conformational changes of single
            proteins and nucleic acids. Its weak point is Cy5's photophysics: it blinks into its triplet and dark
            states and bleaches. That is exactly why the self-healing Cy5–COT dyes were made (Self-healing dyes tab).

            ### In this model
            - Move n and watch the panel note: distance, the angle between the two dipoles, κ², and E as drawn
              versus E for freely tumbling dyes. Light Cy3 (S) and imagine half of that energy jumping to Cy5.
            - The DNA repeats the A5·T20 base pair of a real crystal structure (PDB 1BNA) with its own helical step.
              That step is a little overwound (37.4° and 3.13 Å per pair, so ~9.6 pairs per turn instead of ~10.5),
              so the κ² oscillation here repeats a bit faster than in real DNA. No H on the DNA.
            - The dyes are stacked 3.4 Å beyond the end pairs, turned to fit without clashes; their tethers to the 5′
              phosphates aren't drawn. R₀ = 5.4 nm is an assumed typical value.
            - The dyes here are perfectly rigid, so κ² swings all the way between ~0 and ~1 (at 12 pairs FRET drops
              to 13% where tumbling dyes would give 83%). Real stacked dyes wobble and sometimes unstack, which
              averages κ² part of the way back toward 2/3: the measured oscillation is real but gentler.

            Sources: Förster, Ann. Phys. 1948. Iqbal et al., PNAS 105, 11176 (2008). Drew et al., PNAS 78, 2179
            (1981) (PDB 1BNA). Altman et al., Nat. Methods 9, 68 (2012).
            """)),
    dict(file="aequorin.xyz", name="Aequorin: calcium-triggered light (BRET to GFP)",
         desc="The jellyfish photoprotein that makes the blue light GFP turns green. Holds coelenterazine as a "
              "peroxide, ready to burst when Ca²⁺ binds. PDB 1EJ3.",
         cofactors=dict(pdb="1EJ3", chains=["A"], on=["Coelenterazine"],
                        pigments={"CZH": dict(name="Coelenterazine", kind="dye", color="#4aa3ff", number=False)}),
         info=info("""
            ### What it is
            Aequorin is the calcium-sensitive photoprotein of the jellyfish *Aequorea victoria*, the same animal GFP
            comes from. It holds its substrate, **coelenterazine**, already activated with oxygen as a
            2-hydroperoxide (the glowing site here), and waits.

            ### The photophysics
            This is **chemiluminescence**: the excited state is made by a chemical reaction, not by light. When Ca²⁺
            binds to aequorin's EF-hand sites, the protein changes shape and the peroxide attacks its own carbon,
            forming a strained four-membered dioxetanone. It falls apart into CO₂ and coelenteramide, and the
            energy released leaves the coelenteramide in its singlet excited state. Alone, it emits blue light,
            ~465 nm.

            In the jellyfish it never gets the chance: GFP sits close to aequorin, and the excited coelenteramide
            hands its energy to GFP's chromophore by resonance transfer (**BRET**, bioluminescence resonance energy
            transfer: FRET with a chemically made donor). GFP then emits green, ~509 nm. GFP's whole natural job is
            to be a BRET acceptor.

            ### Why it matters for imaging
            - Aequorin was one of the first calcium indicators: light output directly reports Ca²⁺, with no
              excitation light and so no background autofluorescence or phototoxicity.
            - BRET is the basis of luciferase-based sensors (e.g. NanoLuc fused to a fluorescent protein), which
              read out interactions or Ca²⁺ without illumination.
            - Shimomura found both aequorin and GFP in the 1960s; GFP's Nobel Prize (2008) grew from this.

            ### In this model
            Crystal structure PDB 1EJ3 (2.3 Å), protein drawn faint, coelenterazine-2-hydroperoxide with its H. There
            is no structure of the aequorin–GFP pair; GFP is in the Proteins tab.

            Sources: Head et al., Nature 405, 372 (2000) (PDB 1EJ3). Shimomura, Johnson & Saiga, J. Cell. Comp.
            Physiol. 59, 223 (1962).
            """)),
    # ------------------------------------------------------------------ light harvesting
    dict(file="lh2.xyz", name="LH2: two rings of bacteriochlorophyll",
         desc="Purple-bacterial antenna: 9 subunit pairs holding a loose B800 ring and a tight B850 ring of "
              "bacteriochlorophyll a, with carotenoids spanning both. PDB 1NKZ.",
         cofactors=dict(pdb="1NKZ", assembly=1, keep=["CXM"], labels="lh2", alpha="ACE", beta="BDF",
                        on=["B850 α1", "B850 β1", "B850 α2"],
                        pigments={"BCL": CHL,
                                  "RG1": dict(name="Rhodopin glucoside", kind="carotenoid", core="polyene")}),
         info=info("""
            ### What it is
            Light-harvesting complex 2 of the purple bacterium *Rhodoblastus* (formerly *Rhodopseudomonas*)
            *acidophilus*: nine pairs of small membrane helices (α and β) arranged in a ring, holding 27
            bacteriochlorophyll a (BChl) and 18 carotenoids (rhodopin glucoside).

            ### The photophysics
            The BChls form two rings named by where they absorb:
            - **B800**: 9 BChls lying flat, ~21 Å apart. Weakly coupled, they act roughly as separate molecules
              and absorb at 800 nm.
            - **B850**: 18 BChls standing up, packed so tightly (~9 Å Mg to Mg; their rings 3.5 Å apart here) that
              their excited states mix. The excitation is shared over several molecules at once (an **exciton**),
              which lowers its energy and shifts the absorption from ~800 to 850 nm.
            Energy flows downhill from B800 to B850 in under a picosecond, then on to the LH1 ring and the reaction
            centre. B800 → B850 is Förster-like hopping; within B850 it's coherent sharing. Comparing the two rings
            shows when the simple point-dipole FRET picture works and when it fails.

            The carotenoids touch both rings. They absorb blue-green light (their bright S₂ state) and pass it to
            the BChls, and they take BChl triplets before those can make singlet oxygen: the same TTET trick as the
            self-healing dyes.

            ### Why it matters for imaging
            Excitonic coupling is how dye aggregates (J- and H-aggregates) and some NIR imaging probes get their
            shifted, narrow, intense bands. LH2 is the textbook natural example.

            ### In this model
            The full ring from the PDB's biological assembly (the deposited file holds a third of it), protein
            drawn faint, pigments with H. BChls are labelled by ring (B800 1-9, B850 α1-9 and β1-9, numbered around
            the ring) and carotenoids 1-18. Three neighbouring B850 BChls start lit.

            Sources: Papiz et al., J. Mol. Biol. 326, 1523 (2003) (PDB 1NKZ). McDermott et al., Nature 374, 517
            (1995).
            """)),
    dict(file="pcp.xyz", name="Peridinin–chlorophyll protein (PCP)",
         desc="Dinoflagellate antenna with four carotenoids (peridinin) per chlorophyll. Here the carotenoids are the "
              "light harvesters, not just the protection. PDB 1PPR.",
         cofactors=dict(pdb="1PPR", chains=["M"], on=["Peridinin 612"],
                        pigments={"PID": dict(name="Peridinin", kind="carotenoid", core="polyene"),
                                  "CLA": dict(name="Chl a", kind="chlorophyll", core="macrocycle")}),
         info=info("""
            ### What it is
            A water-soluble antenna protein of dinoflagellates (marine algae). Each monomer holds 8 peridinins and
            2 chlorophyll a, a 4 : 1 carotenoid-to-chlorophyll ratio, in a protein shell. This is chain M of the
            trimer.

            ### The photophysics
            In the sea, red light is absorbed by the water within metres; what is left is blue-green, which
            chlorophyll absorbs poorly. Peridinin fills that gap: it absorbs strongly at ~470–550 nm and passes the
            energy to chlorophyll, which emits at ~670 nm.

            The catch, met with Self-healing JF525: a carotenoid's lowest singlet (S₁) is dark and usually short-lived,
            which makes it a poor energy donor. Peridinin is unusual: it carries a carbonyl in its conjugated chain,
            which gives its S₁ an **intramolecular charge-transfer (ICT)** character. That makes the S₁/ICT state
            more strongly coupled to chlorophyll, so peridinin hands over most of its energy (roughly 90% overall).
            The peridinins sit in direct contact with the chlorophylls (3.7 Å between conjugated atoms here), so
            they also take chlorophyll triplets.

            ### Why it matters for imaging
            It shows that the "quencher or harvester?" question for a carotenoid next to a dye is decided by the
            singlet energies and by coupling, and that a polar group in the chain changes the answer. The same
            charge-transfer idea underlies push–pull dyes and solvatochromic probes.

            ### In this model
            PDB 1PPR (2.0 Å), monomer M, protein drawn faint, pigments with H. Each peridinin's conjugated chain
            glows (orange); Peridinin 612, the one touching Chl a 601, starts lit.

            Sources: Hofmann et al., Science 272, 1788 (1996) (PDB 1PPR). Polívka & Sundström, Chem. Rev. 104, 2021
            (2004).
            """)),
    dict(file="fmo.xyz", name="FMO protein: an energy funnel",
         desc="Green-sulfur-bacterial antenna that wires the chlorosome to the reaction centre: 7 bacteriochlorophylls "
              "whose energies the protein tunes into a funnel. PDB 3EOJ, 1.3 Å.",
         cofactors=dict(pdb="3EOJ", chains=["A"], pigments={"BCL": CHL}),
         info=info("""
            ### What it is
            The Fenna–Matthews–Olson protein of green sulfur bacteria (here *Prosthecochloris aestuarii*), the first
            chlorophyll protein whose structure was solved (1975). It is a trimer; each monomer holds 7
            bacteriochlorophyll a (an 8th, loosely bound one is seen in some crystals, not in this one) in a
            protein "taco shell". It sits between the huge chlorosome antenna and the reaction centre.

            ### The photophysics
            All seven pigments are the same molecule, yet they absorb at different wavelengths (about 790–825 nm).
            The protein tunes each one's **site energy** with nearby charges, hydrogen bonds and ring distortions.
            Combined with their couplings, this makes an **energy funnel**: excitation entering near the chlorosome
            side flows down through the network to the lowest-energy pigments, next to the reaction centre, in a few
            picoseconds.

            FMO is also famous for the **quantum coherence** debate. In 2007, two-dimensional electronic
            spectroscopy saw oscillations ("quantum beats") lasting hundreds of femtoseconds, read as energy moving
            as a wave. Later work attributes the long-lived beats mainly to vibrations of the pigments. Electronic
            coherence between excitons is real but short-lived, and transport is well described as fast hopping
            between partly delocalised states.

            ### Why it matters for imaging
            The same molecule shows different colours in different protein pockets. That is how fluorescent proteins
            are tuned too: the chromophore stays the same and its surroundings shift its colour.

            ### In this model
            PDB 3EOJ at 1.3 Å, the sharpest structure in this tab (compare LHCII at 2.72 Å). Protein drawn faint;
            the BChls with H are labelled by PDB residue number (371–377). The literature numbers them BChl 1–7.

            Sources: Tronrud et al., Photosynth. Res. 100, 79 (2009) (PDB 3EOJ). Fenna & Matthews, Nature 258, 573
            (1975). Engel et al., Nature 446, 782 (2007).
            """)),
    # ------------------------------------------------------------------ triplets put to work
    dict(file="aslov2.xyz", name="LOV2 domain: a triplet that makes a bond",
         desc="Plant blue-light sensor (oat phototropin LOV2). Its flavin's triplet bonds to a cysteine, and that "
              "covalent bond is the signal. Dark state, PDB 2V0U.",
         cofactors=dict(pdb="2V0U", chains=["A"], on=["FMN"],
                        pigments={"FMN": dict(name="FMN", kind="flavin", number=False)}),
         info=info("""
            ### What it is
            The LOV2 (Light, Oxygen, Voltage) domain of oat phototropin 1, a plant blue-light receptor. It binds a
            flavin mononucleotide (FMN) non-covalently, with a cysteine (Cys450) next to the flavin.

            ### The photophysics
            Blue light (~450 nm) excites the flavin to S₁, which crosses efficiently to its **triplet** within
            nanoseconds. Usually a triplet is a liability (it makes singlet oxygen). Here the protein uses it:
            within microseconds the triplet reacts with the nearby cysteine (most likely by taking its hydrogen,
            then the two radicals joining) to form a **covalent bond between the cysteine's sulfur and the flavin's
            C4a**. The adduct no longer
            absorbs blue light or fluoresces (its absorption moves to ~390 nm). It holds for seconds to minutes, then
            breaks thermally and the domain resets.

            The adduct strains the protein and unwinds a helix (Jα) packed against it, which switches on the
            phototropin kinase: the photochemistry becomes a mechanical signal.

            ### Why it matters for imaging
            - LOV domains are the core of many **optogenetic tools**: photoactivatable proteins, light-induced
              dimerisers (iLID), LOV-based transcription switches. They work because the triplet makes a bond.
            - Remove the cysteine and the triplet has nowhere to go but oxygen: see miniSOG, next.
            - LOV-based fluorescent proteins (iLOV) emit green and work without oxygen to mature their chromophore.

            ### In this model
            Dark state, PDB 2V0U (1.4 Å), protein drawn faint, FMN with its H (only its conjugated isoalloxazine
            ring glows). Cys450's sulfur has two conformations in the crystal, 4.5 Å (the one drawn) and 3.7 Å from
            the flavin's C4a.

            Sources: Halavaty & Moffat, Biochemistry 46, 14001 (2007) (PDB 2V0U). Swartz et al., J. Biol. Chem. 276,
            36493 (2001).
            """)),
    dict(file="minisog.xyz", name="miniSOG: a flavin built to make singlet oxygen",
         desc="An engineered LOV domain without its cysteine, so the flavin triplet passes its energy to O₂. Used to "
              "label cells for electron microscopy and to kill them with light. PDB 6GPU.",
         cofactors=dict(pdb="6GPU", chains=["A"], on=["FMN"],
                        pigments={"FMN": dict(name="FMN", kind="flavin", number=False)}),
         info=info("""
            ### What it is
            miniSOG (mini Singlet Oxygen Generator) is a LOV domain from *Arabidopsis* phototropin 2, engineered so
            its flavin does the opposite of a self-healing dye: make as much singlet oxygen as possible. The key
            mutation replaces the adduct-forming cysteine with a glycine.

            ### The photophysics
            Blue light → flavin S₁ → intersystem crossing to T₁, as in LOV2. With no cysteine to react with, the
            triplet lives long enough (microseconds) for O₂ to reach it. **Triplet–triplet energy transfer to
            oxygen** then makes **singlet oxygen** (¹O₂, ~0.98 eV above ground). It is the same TTET as COT taking a
            dye's triplet, except the acceptor is O₂ and the product is a reactive oxidant.

            How much singlet oxygen is the interesting part. It was first reported at ~0.47 per photon absorbed,
            later measured nearer ~0.03: the protein shell shields the flavin from O₂ and quenches its triplet. The
            crystal-structure paper behind this model ("Tailing miniSOG") studies this complex photophysics,
            including how light alters the flavin while miniSOG works. Later variants (e.g. SOPP) reach much higher
            yields.

            ### Why it matters for imaging
            - **Correlative light and electron microscopy**: singlet oxygen polymerises diaminobenzidine into an
              electron-dense stain exactly where miniSOG is, so a genetically encoded tag shows up in the EM.
            - **CALI / photoablation**: kill or inactivate tagged proteins or cells with light.
            - The flip side of self-healing: every dye's triplet + O₂ makes some ¹O₂, which is the phototoxicity
              that the COT conjugates and Gentle Rhodamines reduce.

            ### In this model
            PDB 6GPU (1.17 Å), protein drawn faint, FMN with H. There is no cysteine by the flavin's C4a: the
            nearest residues here are Leu84 and Gly40.

            Sources: Torra et al., Sci. Rep. 9, 2428 (2019) (PDB 6GPU). Shu et al., PLoS Biol. 9, e1001041 (2011).
            Ruiz-González et al., J. Am. Chem. Soc. 135, 9564 (2013).
            """)),
    dict(file="killerred.xyz", name="KillerRed: a phototoxic fluorescent protein",
         desc="A red fluorescent protein that makes reactive oxygen under green light, through a water-filled "
              "channel to its chromophore. PDB 2WIQ.",
         cofactors=dict(pdb="2WIQ", chains=["A"], on=["Chromophore"],
                        pigments={"CRQ": dict(name="Chromophore", kind="fp", color="#ff5a4a", number=False)}),
         info=info("""
            ### What it is
            KillerRed is a red fluorescent protein (excitation ~585 nm, emission ~610 nm) evolved from a hydrozoan
            chromoprotein for phototoxicity. It is far more phototoxic than other fluorescent proteins.

            ### The photophysics
            Most fluorescent proteins wrap their chromophore tightly in an 11-stranded β-barrel, which keeps oxygen
            and water away and makes them bright and photostable. KillerRed's barrel has a **water-filled channel**
            running from the outside to the chromophore. Under green light, the excited chromophore (its singlet, and
            a triplet or radical formed from it) transfers an electron to O₂ arriving through the channel, making
            mainly **superoxide** (type I photosensitisation) rather than singlet oxygen. The same photochemistry
            bleaches KillerRed as it works.

            Compare miniSOG, whose flavin makes singlet oxygen by energy transfer (type II). Two different routes
            from an excited chromophore to damage, and either can be engineered in or out.

            ### Why it matters for imaging
            - Light-controlled cell killing and chromophore-assisted light inactivation (CALI) of tagged proteins.
            - A lesson for every fluorescent-protein label: the barrel's job is to isolate the chromophore, and a
              leaky barrel means phototoxicity and bleaching.

            ### In this model
            PDB 2WIQ (2.0 Å), protein drawn faint. The chromophore is made from Gln65–Tyr66–Gly67 (a DsRed-type red
            chromophore) and is part of the protein chain; its H were placed using its peptide bonds. Waters are
            not drawn, so the channel shows as a gap in the barrel wall.

            Sources: Carpentier et al., FEBS Lett. 583, 2839 (2009) (PDB 2WIQ). Bulina et al., Nat. Biotechnol. 24,
            95 (2006).
            """)),
    # ------------------------------------------------------------------ fluorogenicity and sensing
    dict(file="jf525-open.xyz", name="JF525, open: the fluorescent zwitterion",
         desc="Janelia Fluor 525 in its open form: a flat xanthene conjugated end to end, absorbing at 525 nm. "
              "Compare the closed spirolactone next.",
         smiles="[O-]C(=O)c1ccc(cc1C2=C3C=CC(=[N+]4CC(F)(F)C4)C=C3Oc5cc(N6CC(F)(F)C6)ccc25)C(=O)O",
         source="vendor", formula="C27H18F4N2O5",
         blocks=[dict(label="JF525 (open)", kind="dye", on=True)],
         info=info("""
            ### What it is
            Janelia Fluor 525, a rhodamine with 3,3-difluoroazetidine groups instead of dimethylamines. This is its
            **open zwitterion**: positive iminium on the xanthene, negative carboxylate on the pendant ring.

            ### The photophysics
            Rhodamines live in an equilibrium between two forms:
            - **open** (this one): the xanthene's central carbon (C9) is sp², so the π system runs across all three
              rings from one nitrogen to the other. Absorbs at 525 nm, emits at 549 nm, bright.
            - **closed spirolactone** (next entry): the carboxylate attacks C9, making it sp³. That breaks the π
              system into two separate aminophenyl halves, which absorb only in the UV. Colourless, not fluorescent.
            The balance is the lactone–zwitterion constant, K_L–Z. Polar, protic surroundings favour the open form;
            nonpolar ones favour the closed one.

            The azetidines matter twice. The small, strained rings resist twisting, which suppresses a twisted
            charge-transfer (TICT) decay path and makes the dye brighter than TMR. The fluorines pull electron
            density away from the nitrogens, which blue-shifts the dye (JF549 → JF525) and pushes K_L–Z toward the
            closed form.

            ### Why it matters for imaging
            - **Fluorogenicity**: tune K_L–Z so the dye stays mostly closed (and invisible) free in the cell, then
              opens and lights up when bound to its target, such as HaloTag. That gives less background without
              washing.
            - The closed, uncharged form also crosses membranes more easily, which helps cell permeability.
            - This is the dye in the Voltron525 and HaloTag models (Proteins tab).

            ### In this model
            RDKit + MMFF94s geometry. The glow sits across the xanthene, end to end.

            Sources: Grimm et al., Nat. Methods 12, 244 (2015); Grimm et al., Nat. Methods 14, 987 (2017).
            """)),
    dict(file="jf525-closed.xyz", name="JF525, closed: the colourless spirolactone",
         desc="The same atoms as the open form, with the carboxylate bonded to C9. The xanthene's conjugation is "
              "broken in two, so it no longer absorbs visible light.",
         smiles="O=C1OC2(c3ccc(N4CC(F)(F)C4)cc3Oc3cc(N4CC(F)(F)C4)ccc32)c2cc(C(=O)O)ccc21",
         source="drawing", formula="C27H18F4N2O5",
         blocks=[dict(label="JF525 (closed)", kind="dye", on=True)],
         note="The spirolactone isomer of the open zwitterion: the pendant ring's carboxylate (ortho to the xanthene) "
              "closed onto C9. Same formula.",
         info=info("""
            ### What it is
            JF525's other form: the carboxylate on the pendant ring has bonded to the xanthene's central carbon
            (C9), closing a five-membered lactone ring. Same atoms, same formula as the open form. It's an isomer,
            and the two interconvert in solution.

            ### The photophysics
            Making C9 sp³ cuts the π system in two. Each half is an aminophenyl ether, which absorbs only in the
            UV, so the closed form is **colourless and non-fluorescent**. C9 also turns the pendant ring
            perpendicular, so the molecule is no longer flat.

            Which form wins depends on the surroundings (the lactone–zwitterion equilibrium, K_L–Z). Water and
            hydrogen bonding stabilise the charges of the open form; membranes, nonpolar pockets and aggregation
            favour the closed form. Binding to a protein like HaloTag can pull the equilibrium open, and the dye
            turns on.

            ### Why it matters for imaging
            - **Fluorogenic labels** (JF dyes, MaP dyes, SiR): background stays dark until the dye binds.
            - **Spontaneous blinking** for super-resolution (e.g. HMSiR) uses a similar ring opening and closing,
              tuned so only a few molecules are open at any moment.
            - Cell permeability: the neutral closed form crosses membranes.

            ### In this model
            RDKit + MMFF94s geometry of the spirolactone, drawn from the open form's structure. Compare the glows.
            In the open form it sits across the xanthene, and the pendant ring stays dark. Here C9 is dark, and the
            glow spreads over three separate pieces: both aminophenyl halves and the lactone ring. With the π system
            cut, no single chromophore is left, and the simple Hückel picture behind the glow shares the excitation
            among pieces of similar energy, each of which would only absorb in the UV.

            Sources: Grimm et al., Nat. Methods 14, 987 (2017). Lukinavičius et al., Nat. Chem. 5, 132 (2013).
            """)),
    dict(file="gcamp6m.xyz", name="GCaMP6m: a calcium sensor that switches a proton",
         desc="Circularly permuted GFP fused to calmodulin and its binding peptide. Calcium binding closes a gap in "
              "the barrel and turns the chromophore bright. PDB 3WLD.",
         cofactors=dict(pdb="3WLD", chains=["A"], keep=["CA"], on=["Chromophore"],
                        pigments={"CRO": dict(name="Chromophore", kind="fp", number=False)}),
         info=info("""
            ### What it is
            GCaMP6m, the calcium indicator behind most neuronal imaging: a **circularly permuted EGFP** (its ends
            moved to the middle of the barrel) with **calmodulin** (CaM) on one side and the **M13 peptide** that
            CaM grabs on the other.

            ### The photophysics
            The chromophore (Thr65-Tyr66-Gly67, made by the protein itself) has two forms, set by the proton on its
            phenol:
            - **protonated (neutral)**: absorbs ~400 nm, barely excited by the usual 488 nm laser, dim;
            - **deprotonated (anionic)**: absorbs ~490 nm, bright green emission at ~510 nm.
            Circular permutation leaves an opening in the barrel right next to the phenol. Without Ca²⁺, water gets
            in and the chromophore stays mostly protonated, so it's dim. When Ca²⁺ binds, CaM wraps the M13 peptide
            and the complex packs against that opening, sealing it and reorganising the hydrogen bonds. The
            chromophore's pKa drops, it loses the proton, and the protein gets many times brighter.

            So GCaMP doesn't change how its chromophore emits. It changes **which ground state is populated**: a
            ground-state protonation equilibrium, the same idea as rhodamine's open/closed equilibrium.

            ### Why it matters for imaging
            ΔF/F from GCaMP is how most labs record neural activity. Its kinetics (rise and decay), affinity and
            dynamic range all come from this mechanism, and the GCaMP6/7/8 generations tuned the interfaces around
            the chromophore. Pair it in your head with Voltron (Proteins tab): a voltage sensor that works by FRET to
            a rhodopsin instead.

            ### In this model
            Ca²⁺-bound GCaMP6m, PDB 3WLD (2.7 Å). Protein and the 4 Ca²⁺ ions drawn faint, the chromophore glowing.
            It is part of the protein chain, and its H (drawn here with the phenol protonated, as the dictionary has
            it) were placed using its peptide bonds.

            Sources: Ding et al., Sci. China Life Sci. 57, 269 (2014) (PDB 3WLD). Chen et al., Nature 499, 295
            (2013). Akerboom et al., J. Biol. Chem. 284, 6455 (2009).
            """)),
]
