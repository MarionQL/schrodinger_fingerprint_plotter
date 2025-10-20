# schrodinger_fingerprint_plotter
Visualizing the csv output from Schrodinger interaction fingerprints as a heatmap (multiple ligands) or bar graph (single ligand).

This code:
1. Parses Schrodinger fingerprint CSVs containing interaction data across multiple ligand poses and accross multiple ligands
2. Extracts ligand names from common naming patterns in the first column
3. Can filter for certain interaction types identified by schrodinger
4. Maps residue labels using an input PDB file (this is the protein that you docked your ligands to)
5. Counts interaction frequencies for each ligand-residue pair
6. Generates bar plots and heatmaps for selected interaction types

To run this script you need two things
1. The CSV from schrodinger
   - this csv can have as many ligands as you want
   - if you make a heatmap, all ligands will be included, if you make a bar graph, each ligand will be placed in its own bargraph
   - In order to identify ligand names, the first column of the CSV must be in the following format:
      1. `{protein}_{ligand}_{pose_number}`
      2. {ligand}_{pose_number}
      3. {ligand}
      4. {protein}_{ligand}
   * NOTE: Pose number must be an integer (not pose1, pose2, pose3)
   - this file should also be named {protein}_rest_of_csv.csv
      1. this is so that file outputs can be named with the protein name they are associated

2. The PDB of your protein

    - The script uses this protein file to extract residue names and is used to order the residue labels
    - The schrodinger output only includes the chain ID and the residue number, so that is why we need the pdb of the protein

Flags:

Required:

1. -i --interaction: this filters for different interactions, contact will give you all interactions, you can also get backbone or sidechain, or interaction type (polar, charged, donor, acceptor, hydrophobic, aromatic)
  - these interactions are all determined by schrodinger, so if there are errors in interaction type, there may be a problem with how you ran fingerprinting in Schrodinger
2. -g --graph: You can choose to generate a heatmap or a bar graph, examples of what these options look like are in the example_data folder

Optional:
1. -ic --ignore-chain: this flag will remove the chain letter from the residue name in the graph, this is not recommended if the protein has multiple chains
2. -s --show: this will show the graph in matplotlib as the graphs are created

In the example_data folder, there is a protein.pdb and a fingerprint.csv that work with this script and the heatmap and bar graphs that are output
These are the commands used to generate those graphs:

`python plot_schrodinger_fingerprints.py mfsd2b_fingerprint.csv mfsd2b.pdb -i contact -g bar -ic`

`python plot_schrodinger_fingerprints.py mfsd2b_fingerprint.csv mfsd2b.pdb -i contact -g heatmap -ic`

Aesthetic Issues:
If you have a long ligand name or many residue contacts, sometimes the size of the figure might not be wide enough. If this is an issue and labels are getting cut off, try increaseing the figure size in either the heatmap or bar graph functions depending on which graph you are having trouble with. The line you edit is this: plt.figure(figsize=(12, 8)) - to make the figure wider, change 12 to a larger number
