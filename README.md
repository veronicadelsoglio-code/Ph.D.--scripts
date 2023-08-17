# Ph.D.--scripts
This folder contains all the Rosetta .xml scripts used for protein backbone and sequence design, for modeling and Python scripts for data analysis


STEPS OF THE DESIGN STRATEGY

******** A) TRUNCATION OF THE FULL-LENGTH PROTEIN BUNDLE ********

In this step the 3H5L protein bundle length was shortened in order to still enable the coordination with one b-type heme molecule, simplifying the overall design process. The blueprint file is a file that defines the secondary structure and amino acid identities for each residue. A blueprint file for the monomeric 3H5L bundle was generated by running getBluePrintFromCoords.pl script:
(Run command 1)
getBluePrintFromCoords.pl* -pdbfile [starting pdb] > [blueprint name]

The residues that needed to be removed from the bundle to create the truncated protein form were deleted from the blueprint file. Additionally, also the loop regions that con- nect the truncated helices were specified. Subsequently, the RosettaRemodel application, which handles fixed-backbone design, was run using the command:
(Run command 2)
remodel.linuxiccrelease -s input.pdb -remodel:blueprint blueprintfile.bp -ex1 -ex2 -num_trajectory 10

The Remodel application, available in Rosetta/main/source/bin directory, required the input pdb of the long 3H5L bundle structure and the blueprint file generated earlier to perform its design function. Rotamer re-sampling was favoured adding the flags -ex1 -ex2. Due to the stochastic Monte Carlo process used in sidechain selections, multiple design runs were created to ensure convergence in the search process. Setting the flag -num_trajectory 10 a reasonable amount of Monte Carlo sampling could have been per- formed, but only the few decoys with lowest energies were outputted. At the end of this process, the truncated version of the 3H5L bundle was obtained.


+++ (files used: blueprintfile.bp)










******** B) 3H5L TRUNCATED BUNDLE SEQUENCE REDESIGN AND OPTIMIZATION ********

The main goal of this step was to obtain an energetically optimised sequence folding as the initial model (T-3H5L bundle). A process of sequence design was conducted using the PackRotamer Mover. The packer is a powerful tool, the primary one in Rosetta, for designing amino acid sequences. The packer default behaviour is to design every position, considering all possible rotamers for each of the 20 canonical amino acids. For that reason the RestrictToRepacking TaskOperation - which favourites the packer to only pack existing side-chains without changing the residues identity - was not allowed in the designing process. A Rosetta XML script was set up, and the PackRotamer mover was executed. The script was run for 100 times in order to generate multiple outputs:
(Run command 3)
/Rosetta/.../rosetta_script.linuxrelease -s input.pdb -parser:protocol SequenceDesign.xml @enzflags

Subsequently, a Python script was used to analyse the 100 outputs both in terms of REU and of residue identity. This means for every position of the backbone, which residues the Packer chooses, how much frequently they occur considering all the outputs, their chemical nature and REU was investigated. The final sequence was build selecting the most frequent - low energy residue at each position. Also the residues ́s chemical identity was analysed - considering the homodimeric final structure.

+++ (files used: SequenceDesign.xml, enzflags, Mutation_Script.xml, Resfile)














******** C) HOMODIMERIZATION ********

This step aimed to generate a C2 symmetric homodimer with the rotation axis parallel to the long helical axis (TH-3H5L). Everything that needed to be provided to Rosetta software about the symmetry of the system, was encoded in the symmetry definition file (C2_symm). The latter was generated according to the instructions provided in the dedicated Rosetta Common section. An XML Rosetta script was run using SetupForSymmetry mover, which included the information of the symmetry definition file, generating the final homodimeric output structure. The script was run using the command:
(Run command 4)
/Rosetta/.../rosetta_script.linuxrelease -s input.pdb -parser:protocol Symm.XML @ enzflag






******** D) HEME INCLUSION ********

The b-type heme molecule was included inside the homodimeric scaffold in a bis-HIS fashion. To achieve this, the ligand molecule (HEM.pdb) was manually positioned within the protein scaffold, ensuring the desired orientation. The geometries between specific atoms of the Fe-porphyrin ring and the HIS coordinating residues were measured and saved in a CST file. The latter was divided into two blocks: one block specified the geometries from the Fe atom (heme) to the HIS residue of the first protein bundle, the second one specified the geometries from the Fe atom (heme) to the HIS of the second protein bundle. The measured geometries included the following:
1) Distance (Å) between the nitrogen atom of the respective HIS residues and the Fe atom of the b-type heme ligand;
2) Angle A, between three atoms: nitrogen (HIS) - Fe (heme) - nitrogen (one of the four planar nitrogen atoms coordinating the metal ion);
3) Angle B, between three atoms: Fe (heme) - nitrogen (HIS) - carbon (one of the two carbon atoms adjacent to the preceding nitrogen atom (HIS));
4) Torsion A, angle between four atoms: nitrogen (HIS) - Fe (heme) - nitrogen (one of the four nitrogen atoms in the pyrrole rings coordinating the metal ion) - carbon (one of the two carbon atoms adjacent to the preceding nitrogen atom);
5) Torsion B, angle between four atoms: Fe (heme) - nitrogen (HIS) - carbon (one of the two carbon atoms adjacent to the nitrogen atom (HIS)) - nitrogen (next to the preceding carbon atom);
6) Torsion AB, angle between four atoms: nitrogen (HIS) - Fe (heme) - nitrogen (one of the four nitrogen atoms coordinating the metal ion) - carbon (one of the two carbon atoms adjacent to the preceding nitrogen atom);

To validate the geometries specified in the CST file well describe the coordination event, ensuring no clashes occur with the side chain of the surrounding residues, the matching run was performed. The latter uses the matcher application, which scans for positions where a binding site for the ligand of interest within a protein target can be accomodated. The HIS 12 residue of the respective two bundles were found as a possible place where the described geometries could be placed. The files used to perform the run were: the CST file, the HEM.params file, the enzflag file. The HEM.params file was originated starting from the HEM.mol2 file, originated from the HEM.pdb one. The command used was:
(Run command 5) /ROSETTA/main/source/scripts/python/public/molfile_to_params.py HEM.mol2

The Matcher execution was performed with the command:
(Run command 6)
/rosetta/main/source/bin/match.linuxgccrelease -s PDBFILE.pdb -match:lig_name HEM -match:geometric_constraint_file CST_file.cst -match:scaffold_active_site_residues_for_ geomcsts Posfile.pos -match::enumerate_ligand_rotamers false @enzflags

Once the Matching run confirmed the HEM coordination can be feased, the position of the ligand was fixed. To this purpose, a RosettaScript script was run. The ligand coordinating HIS residues were mutated in not-protonated residues using the MutateResidue mover. Subsequently, the AddOrRemoverMatchCsts and EnzRepackMinimize movers favoured the incorporation of the geometric CST into the pose, followed by the design/repack and minimization of the protein-ligand interface. The files used to perform the run were: Rosetta.xml script, HEM.params file, CST file, enzflag file. The command used to run the RosettaScript script was:
(Run command 7) /ROSETTA/main/source/bin/rosetta.script.linuxgccrelease - s input.pdb -parser:protocol FinalDesign.xml -extra_res_fa HEM.params -match:geometric_constraint_file CST_file.cst @enzflags

The score of the output structure were evaluated.


+++ (files used: HEM.mol2, HEM.params, FinalDesign.xml, CST_file.cst, enzflags (used in the previous step), Posfile.pos)











******** E) HEME BINDING POCKET (BP) OPTIMIZATION ********

The aim of this step was to optimise the protein BP to enhance ligand binding. Since the heme is an highly hydrophobic molecule, the work was conducted with the goal of increasing the BP hydrophobicity. More specifically, three different RosettaScript scripts were written and executed, generating output structures showing BP with different characteristics. Their action scheme was the same; Rosetta FastDesign mover was used to redesign just specific residues constituting the BP. The selected BP positions were indicated within the <RESIDUESELECTOR/> tag and invoked through the RestrictToResidueProperties TaskOperations flag. The latter imposed the hydropho- bicity nature of the re-designed residues. The action of the FastDesign mover depended also on the AddCompositionContraint Mover. The latter was applied just to the residues of the BP - which through the <Comp subtag/> constrained the chemical identity of the new residues (e.g. level of aromatic residues, absence of polar charged residues). One of the parameter the FastDesign mover considered during the BP designing process was the aromaticity %. It defines the % amount of possible aromatic residues the mover could fix in the BP region. The aromatic residues, located in the BP positions close to the HIS, can help in its re-orienting, favouring the coordination with the ligand. The respective Rosetta XML scripts were assigned different aromaticity %, which were:

In 1.XML script: aromaticity % fraction of 0.1 (10 % of the residues redesigned by the FastDesign mover could be aromatic);
In 2.XML script: aromaticity % fraction of 0.5 (50 % of the redesigned residues could be aromatic);
In 3.XML script: aromaticity % fraction of of 0.8 (80 % of the new residues could be aromatic);

Not only the aromaticity % differentiated the three xml scripts; the selected BP positions was also considered. In the first two xml scripts the redesigned positions were ”9, 13, 16, 72, 75, 79”. In the third script, the BP was expanded to include position 82. Additionally, an ARG residue was enforced at position nine to facilitate coordination with the ligand’s isopropionate groups. The protein is an homodimer, hence the design process was favoured in a symmetrical way between chain A (one protein ́s monomer) and chain B (the other protein ́s monomer), using the SetUpForSequenceSymmetry Mover. Subse- quently, the InterfaceAnalyzer mover was utilized to evaluate the interface score between the protein monomers. The <SIMPLE_METRICS/> tag was used for data analysis and data reporting. Their execution was invoked through the RunSimpleMetrics mover and their different score values - Interaction energy of all the structure, Interaction energy of the protein interface, SASA metric (Solvent Accessible Surface Area) of the protein interface - were analysed.


+++ (files used: 1.xml, 2.xml, 3.xml, HEM.params - enzflags (used in the previous steps), MutationPosition9.xml, ResfilePosition9)










******** F) DESIGN OF PROTEIN VARIANTS: H12A, R9I, R9D, R9K ********

To generate TH-3H5L-4_3_3 and TH-3H5L-4_3_4 protein variants (H12A, R9D, R9I, R9K), specific Resfile file - one for each mutant - was generated. The Resfile is a file for specifying mutatable residues for a design run. The sequence position to be re-designed in which other residue was indicated in every generated file. The Resfile setting used was the same for all the generated variants:

NATAA start
[initial position] A (chain) PIKAA [replacing residue]
[initial position] B (chain) PIKAA [replacing residue]

The information from the Resfile was incorporated in a Rosetta XML script where PackRotamerMover was allowed to design sequence mutations in a symmetric man- ner for both protein bundles. Then a minimisation run (min_torsion, min_cartesian) was conducted to optimise the overall protein structures. The Abinitio folding prediction run was performed to confirm the bundle folding was preserved even after the single point mutations.














