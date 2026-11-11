"""About-panel texts for the "Self-healing dyes" tab, keyed by catalog file name (merged into tools/catalog.py).

Numbers are approximate literature values (the same ones as README, "Self-healing dyes"); distances marked "here" are
measured on the model.
"""

from catalog_photophysics import info

SELFHEAL_INFO = {
    "cy5-cot.xyz": info("""
        ### What it is
        **Cy5** with **cyclooctatetraene (COT)** tied to one of its indolenine nitrogens by a flexible tether: the first
        "self-healing" dye (Altman et al. 2012). The same paper also tethered two other protective agents, nitrobenzyl
        alcohol and Trolox, to Cy5; COT is the one that works by taking the triplet.

        ### The photophysics
        Now and then an excited dye crosses into its **triplet state** (T₁) instead of fluorescing. A triplet can't
        emit quickly, so the dye goes dark for as long as it lasts (Cy5: ~60 µs without help), and it is reactive:
        with oxygen it makes singlet oxygen, and either route can end in bleaching.

        A quencher bonded to the dye fixes that by **triplet–triplet energy transfer** (TTET, a Dexter exchange that
        needs the two π systems to touch or be joined by bonds). Two conditions make it work:
        - **Triplet downhill:** COT's T₁ relaxes to ~0.95 eV by flattening, well below Cy5's 1.60 eV, so the triplet
          flows to COT and doesn't come back. COT then drops to its ground state as heat.
        - **Singlet uphill:** COT's S₁ is ~4 eV, far above Cy5's (~1.9 eV), so the fluorescence is left alone.

        On this flexible tether TTET happens when the tether lets the COT bump into the dye. Cy5's triplet lifetime falls
        from ~63 µs to ~1 µs with the shortest tether and ~13 µs with a 13-atom one (Zheng et al. 2012), and the number
        of photons per molecule rises severalfold, up to ~70× without oxygen.

        ### Why it matters for imaging
        The classic fix is millimolar COT (or Trolox, or an oxygen scavenger) in the buffer, which you can't put into a
        living cell. Bonding one COT to each dye gives it a quencher at a huge local concentration without adding
        anything to the sample. Cy5 is the workhorse of single-molecule FRET and super-resolution imaging, where blinking
        and bleaching set the limit.

        ### In this model
        - Cy5 (lit) and COT are separate sites. In simulation mode, the "Self-healing dye" preset counts triplets
          made and quenched; run it here and on plain Cy5 to compare.
        - The paper's exact 12-atom linker is in a supplementary figure we couldn't read, so the tether drawn here
          (hexanoyl amide to aminomethyl-COT) is illustrative. The geometry is one force-field conformer: the real tether
          flops around, which is the whole mechanism.

        Sources: Altman et al., Nat. Methods 9, 68 (2012). Zheng et al., J. Phys. Chem. Lett. 3, 2200 (2012).
        """),
    "gr555-mito.xyz": info("""
        ### What it is
        **Gentle Rhodamine GR555-mito**: tetramethylrhodamine (TMR) with COT attached as an ester on its pendant ring,
        COT–CH₂–O–C(=O)–Ar (Liu et al. 2024). It is a drop-in for TMRM, the methyl-ester dye used to stain mitochondria.

        ### The photophysics
        Same idea as Cy5–COT: the rhodamine's triplet sits well above COT's (~0.95 eV), so COT takes it by TTET; COT's
        singlet is far above the dye's, so the fluorescence is untouched. What the paper measured is the payoff for live
        cells: about **10× less singlet oxygen** than TMRM and **4–5× less phototoxicity**.

        Singlet oxygen needs ~0.98 eV from a triplet. A rhodamine triplet has plenty; COT's relaxed triplet has about
        that or less, so moving the triplet onto COT both shortens it and leaves less energy to hand to oxygen.

        ### Why it matters for imaging
        TMRM is a cation that piles up in mitochondria in proportion to their membrane potential, so it reports on
        mitochondrial health. But imaging it for long makes reactive oxygen that damages the same mitochondria.
        Capping the carboxyl with COT keeps the dye a lipophilic cation (it still goes to mitochondria) and builds the
        protection into the molecule.

        ### In this model
        - TMR (lit) and COT are the two sites; the ester and CH₂ between them are no site.
        - Linkage and position are from the paper; the 3D shape is one RDKit + MMFF94s conformer.

        Sources: Liu et al., ACS Cent. Sci. 10, 1933 (2024).
        """),
    "self-healing-jf525.xyz": info("""
        ### What it is
        Our design: **JF525**, a bright, cell-permeable rhodamine, carrying the triplet quencher photosynthesis uses, a
        **carotenoid**, at the end of a 1,4-phenylene (the TTA Flora linker). The phenylene sits on JF525's pendant ring
        para to the xanthene, so both carboxyls stay: the one that closes the lactone and the one that becomes the
        HaloTag handle. The carotenoid is 7′-apo-β-carotene (β-carotene cut at C7′ and bonded there, all-trans), the
        aryl-ended carotenoid Gust and Moore joined to porphyrins in their carotenoporphyrin dyads.

        ### The photophysics
        - **Triplet condition: met easily.** β-carotene's T₁ is ~0.9 eV, far below a rhodamine's (~1.8–2.0 eV; JF525's
          own isn't reported that we know of). It is also below the ~0.98 eV singlet oxygen needs, which is why
          carotenoids can't sensitize it and are used against it in nature.
        - **Singlet condition: probably not met.** A carotenoid's lowest singlet is the dark S₁ (2¹A_g⁻), which can't
          absorb or emit light but can accept energy. For this conjugation length it is ~1.8 eV, below JF525's S₁
          (~2.3 eV, from 525 / 549 nm). Through the same bonds that carry TTET, JF525's singlet could drop into it, so
          this version would likely be dim.

        In photosynthesis the gap is small because chlorophyll's S₁ is low (~1.85 eV), and even there one extra
        conjugated C=C turns a carotenoid into a chlorophyll singlet quencher (Berera et al. 2006). Shortening the
        polyene raises S₁ (and T₁), so a **shorter carotenoid** is the version that could satisfy both conditions.

        ### Why it matters for imaging
        If the levels can be tuned, a carotenoid is an attractive quencher: it is what nature already uses, it
        destroys singlet oxygen as well as preventing it, and it doesn't have COT's charge-transfer risk. JF525 is a
        common HaloTag dye, so a self-healing version would plug straight into live-cell labelling.

        ### In this model
        - JF525 (lit) and the carotenoid are sites; the phenylene is no site.
        - One RDKit + MMFF94s conformer. JF525's ring and the phenylene twist ~64° apart (a carboxyl is next door), the
          phenylene and polyene ~42°; the molecule is ~3 nm long. Energies are literature estimates, not computed.

        Sources: Gust, Moore et al., carotenoporphyrin dyads (e.g. Acc. Chem. Res. 1993). Berera et al., PNAS 103, 5343
        (2006). Grimm et al., Nat. Methods 14, 987 (2017) (JF525).
        """),
    "lhcii.xyz": info("""
        ### What it is
        **Light-harvesting complex II** (LHCII), the most abundant membrane protein on Earth and the main antenna of
        plant photosystem II: one monomer from the 2.72 Å X-ray structure of the spinach trimer (Liu et al. 2004). It
        holds **14 chlorophylls** (8 a, 6 b) and **4 carotenoids** (two luteins, violaxanthin, neoxanthin). The protein
        is drawn faint.

        ### The photophysics
        - **Light harvesting:** the chlorophylls absorb red and blue light and pass the excitation among themselves
          (Chl b → Chl a within picoseconds) and on to the reaction centre. Carotenoids help: their bright S₂ state
          absorbs blue-green light and hands it to chlorophyll.
        - **Triplets are the danger.** A chlorophyll that doesn't pass its excitation on crosses to its triplet often
          (~60% for Chl a in solution). Its T₁ (~1.3 eV) is above the ~0.98 eV needed to make singlet oxygen, which
          would wreck the protein.
        - **Carotenoids are nature's self-healing quencher.** Their T₁ (~0.9 eV) is below chlorophyll's, so they take
          the triplet by TTET within nanoseconds and turn it into heat, and it is below singlet oxygen's, so they can't
          make it themselves (they also quench any that forms).
        - There's **no linker**: TTET needs van der Waals contact, and the protein holds each lutein right on the
          chlorophylls it guards. Here Lutein 620 is 3.6 Å from Chl a610 and 3.7 Å from a612; Lutein 621 is 3.4 Å from
          a602 and 3.7 Å from a603 (closest conjugated atoms).

        ### Why it matters for imaging
        It is the design every self-healing dye copies: a bright pigment, a quencher whose triplet is lower and whose
        singlet is out of the way, held in contact. LHCII is also a reminder of the limit: carotenoids can quench
        chlorophyll singlets too, and plants are thought to use that to dump excess light (one proposed route of non-photochemical quenching).

        ### In this model
        - Each pigment is a site, but only its conjugated part glows (the macrocycle: Mg, 4 N, 20 C; or the polyene).
          Chl a610, the lowest-energy chlorophyll and the main triplet donor, starts lit.
        - X-ray maps at this resolution can't see hydrogens, so the deposited model has none. The pigments' H are put
          back from the PDB Chemical Component Dictionary; the protein keeps none. Phytol atoms too disordered to be in
          the crystal are missing. Waters, lipids and detergent are removed.
        - Checked against pea LHCII (PDB 2BHW, 2.5 Å): the same contacts within 0.2–0.4 Å.

        Sources: Liu et al., Nature 428, 287 (2004) (PDB 1RWT). Standfuss et al., EMBO J. 24, 919 (2005) (PDB 2BHW).
        """),
    "chl-a610-lutein-620.xyz": info("""
        ### What it is
        Just two pigments cut out of LHCII at their crystal positions: **Chl a610**, the chlorophyll that makes most
        of the complex's triplets, and **Lutein 620**, the carotenoid lying across it that takes them. Whole molecules
        with their H (C₅₅H₇₂MgN₄O₅ and C₄₀H₅₆O₂), phytol tail included; no protein, no linker.

        ### The photophysics
        This is the bare donor–quencher pair to compare with the self-healing dyes. Chl a's T₁ is ~1.3 eV and lutein's
        ~0.9 eV, so the triplet runs downhill to the lutein by TTET (Dexter exchange), which needs the two π systems to
        overlap. Here the lutein's polyene lies across the edge of the chlorophyll's ring, 3.6 Å away: van der Waals
        contact, close enough for electron exchange.

        ### Why it matters for imaging
        It shows what the synthetic dyes have to reproduce with bonds. Nature holds the pair in contact with a protein;
        Cy5–COT does it with a floppy tether that lets them collide; our phenylene designs try a rigid rod, which
        transfers through bonds instead of contact.

        ### In this model
        - Chl a610 (lit) and Lutein 620 are the two sites; only the macrocycle and the polyene glow.
        - Without the protein nothing holds them here: in solution the two would drift apart, and the triplet would have
          to wait for a collision.

        Sources: Liu et al., Nature 428, 287 (2004) (PDB 1RWT).
        """),
    "cy5-phenylene-cot": info("""
        ### What it is
        Our design: COT on Cy5's benzo ring (C5, where sulfo-Cy5 carries a sulfonate) at the end of a rigid rod of
        **n 1,4-phenylenes**, the linker of TTA Flora 4. n = 0 bonds COT straight to the dye. The other indolenine keeps
        its sulfonate and hexanoic-acid handle. Our literature search found no published dye–phenylene–COT (the closest
        are COT bonded to porphyrins, Ostapko et al. 2020).

        ### The photophysics
        The levels are those of Cy5–COT: Cy5's T₁ (1.60 eV) runs down to COT's (~0.95 eV), and COT's S₁ (~4 eV) is out of
        the way. What changes is how the triplet gets there. Estimates, not results:
        - **Through bonds, not by contact.** The rod holds COT 1.5 Å (n = 0) to 18.7 Å (n = 4) from C5, so the two can't
          touch. Through-bond triplet transfer falls off about e^(−0.45 Å⁻¹ · r) along conjugated bridges, so each
          phenylene (4.3 Å) costs ~7×. The flexible Cy5–COT reaches ~10⁶ s⁻¹; a rod has to beat that, so n = 0 or 1 is
          the range to test.
        - **Twists decouple.** Each biaryl here twists ~55° and the COT stands ~80° to its ring: that keeps Cy5's
          spectrum close to normal but also weakens the coupling TTET needs.
        - **Charge transfer is the risk.** COT accepts electrons. Calculations on COT bonded straight to dyes put a
          dye⁺/COT⁻ state below the dye's S₁, which would quench fluorescence. A phenylene weakens this as well as TTET,
          so n = 0 vs n = 1 is the real trade-off.
        - **Conjugation shifts levels.** Conjugating rings onto COT raises its T₁, and an aryl on C5 red-shifts Cy5
          somewhat. Both narrow the gaps, but Cy5's T₁ still clears COT's with room to spare.

        ### Why it matters for imaging
        A flexible tether only quenches when it happens to bring COT to the dye, and it can let COT fold onto the dye,
        where it may quench the fluorescence too. A rod fixes the distance, so the quenching rate would be set by design
        rather than by chance: one fixed rate instead of a spread.

        ### In this model
        - Move n: the note gives the dye–COT distance and the ring twists along the rod.
        - Cy5 (lit) and COT are sites; the phenylenes are no site. The dye's pose is the lowest of 30 MMFF94s conformers
          and is the same for every n; only the rod and COT are searched for each n (see docs/selfheal-report.md).

        Sources: Eng et al., J. Phys. Chem. B 110, 6483 (2006). Chanmungkalakul et al., Chin. Chem. Lett. 35, 109227
        (2024). El Bakouri et al., Chem. Sci. (2025). Ostapko et al., Chem. Eur. J. 26, 16666 (2020).
        """),
}
