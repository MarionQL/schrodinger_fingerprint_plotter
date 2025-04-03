# schrodinger_fingerprint_plotter
Visualizing the csv output from Schrodinger interaction fingerprints as a heatmap (multiple ligands) or bar graph (single ligand).

This code:
1. Parses Schrodinger fingerprint CSVs containing interaction data across multiple ligand poses and accross multiple ligands
2. Extracts ligand names from common naming patterns in the first column
3. Can filter for certain interaction types identified by schrodinger
4. Maps residue labels using an input PDB file (this is the protein that you docked your ligands to)
5. Counts interaction frequencies for each ligand-residue pair
6. Generates bar plots and heatmaps for selected interaction types

To run this script you need to things
1. The CSV from schrodinger
   - this csv can have as many ligands as you want
   - if you make a heatmap, all ligands will be included, if you make a bar graph, each ligand will be placed in its own bargraph
   - In order to identify ligand names, the first column of the CSV must be in the following format:
      1. `{protein}_{ligand}_{pose_number}`
      2. {ligand}_{pose_number}
      3. {ligand}
      4. {protein}_{ligand}
   - this file should also be named {protein}_rest_of_csv.csv
      1. this is so that file outputs can be named with the protein name they are associated

2. The PDB of your protein

    - The script uses this protein file to extract residue names and is used to order the residue labels
    - The schrodinger output only includes the chain ID and the residue number, so that is why we need the pdb of the protein

