# SARS-CoV-2

Currently, the world is being threatened by the SARS-CoV-2, of which the research on associated treatments/antibodies/drugs is also being conducted by many groups/companies, providing massive information about this coronavirus, including the structure of its main proteinase (**Mpro**; PDB Code: [6LU7](https://www.rcsb.org/structure/6lu7)) and its RNA-dependent RNA polymerase (**Rdrp**; [6M71](https://www.rcsb.org/structure/6M71)). Such efforts have facilitated the feasibility of drug screening research, but, of course, before these structures' presence, people had already used models derived from homology modelling, with the **Mpro** of the SARS-CoV as the template.

People havs chosen diverse compound libraries to screen, including FDA-approved drugs, natural products, drugs failing in clinical phases, etc.

## Fragment screening

The British synchrotron center, Diamond, has screened ~1500 fragments against the **Mpro** of SARS-CoV-2 and also provided the [data](https://www.diamond.ac.uk/covid-19/for-scientists/Main-protease-structure-and-XChem/Downloads.html) they obtained. In total, we have 91 crystal structures. As they reported, after their examination, there are 66 hits in the active site, of which, 22 are non-covalent hits and 44 are covalent hits.

## Binding mode analysis

Although the company has already analyzed them for us, I still want to do some little experiments. Here, I will extract the binding modes of the fragments and convert them to a interaction fingerprint-like fingerprint in order to find some patterns. The pseudocode is followed.
```
binding_site = {}

interactions = {'hbonds', 'halogenbonds',
                'pi_stacking', 'salt_bridges',
                'hydrophobic_contacts', 'pi_cation', 'pi_metal'}

for fragment in fragments:
    for interaction in interactions:
        contacted_residues = detect_interaction(fragment, protein, interaction)
        binding_site[interaction] = contacted_residues
```
(For complete implementation, please consult the `classifyLigand.py` and `classifyMols.py` in `src` folder.)

Well, the above code might be a little bit confusing.

- Basically, a binding site is an assemble of amino acid residues contacted by any fragments, which means it will stay at the residue level. But such a bindin site is divided by specific ligand-protein interactions. Therefore, if a residue is contacted by two fragments through two different interactions, this residue will appear twice in the so-called binding site, but in different interaction blocks.

- As for each fragment's fingerprint, it would be staying at the atom level (when detecting pi-stacking/cation/metal interaction, it is at the "ring" level), and still, it is divided by different interactions. But this time, each interaction block might contain at least two identical residues, if this residue interacts with the fragment through two different atoms.

- That is to say, the length of each fingerprint is dependent on the number of the amino acid residues of the binding site (so the length is identical). And each bit of the fingerprint represents the frequency of a ligand interacting a specific residue through a specific interaction. An illustration is followed.

<p align="center"><img src="https://github.com/XiaochangLiu-UoE/SARS-CoV-2/raw/master/Pic/InteractionFP_Illustration.png"></p>

In my implementation, the interaction detection is heavily dependent on [Open Drug Discovery Toolkit (oddt)](https://github.com/oddt/oddt).

Finally, we can use different metrics to evaluate the similarity between fingerprints, such as, Tanimoto and Dice. Also, because the length is identical, the euclidean distance can be suitable, either.

<p align="center"><img src="https://github.com/XiaochangLiu-UoE/SARS-CoV-2/raw/master/Pic/SARS_CoV_2_heatmap_tanimoto_corrected.png"></p>

Actually, in [ODDT](https://github.com/oddt/oddt), there are different interaction fingerprints we can use and evaluate similarity using Tanimoto or Dice.

OK. So what is the point of doing such things? Orginally, I was trying to use it to answer the frequent question "please analyse the ligand binding modes and find some patterns". I believe it would work for common lead screening result analysis, but, it might not perform well when used to deal with fragment screening results.

Fragment screening is a powerful and also popular tool in drug discovery, which may provide much more starting points for lead design than the conventional screening may do.
