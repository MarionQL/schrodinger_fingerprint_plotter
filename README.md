# Schrodinger Fingerprint Plotter

Visualize Schrodinger interaction fingerprint CSVs as heatmaps or bar graphs across multiple ligands and protein conformations.

---

## What This Script Does

1. Parses one or more Schrodinger fingerprint CSVs containing interaction data across multiple ligand poses
2. Extracts ligand names from common naming patterns in the first column
3. Filters for specific interaction types identified by Schrodinger
4. Maps and aligns residue labels using input PDB files — supports multiple PDBs with different residue numbering via sequence-similarity alignment
5. Counts interaction frequencies for each ligand–residue pair across all input files
6. Generates heatmaps or bar graphs for selected interaction types

---

## Requirements

Create a conda environment with the required packages (only needs to be done once):

```bash
conda create -n fingerprints python=3.11 numpy=1.26 pandas seaborn matplotlib
```

Then activate it each time before running:

```bash
conda activate fingerprints
```

---

## Inputs

### CSV Files (`-c`)
One or more Schrodinger fingerprint CSV files. Each CSV is paired with a PDB at the same position — order matters.

- A CSV can contain as many ligands as you want
- The **first column** must follow one of these naming patterns:
  | Format | Example |
  |--------|---------|
  | `{protein}_{ligand}_{pose_number}` | `mfsd2b_compA_1` |
  | `{ligand}_{pose_number}` | `compA_1` |
  | `{ligand}` | `compA` |
  | `{protein}_{ligand}` | `mfsd2b_compA` |

  > **Note:** Pose number must be an integer (`1`, `2`, `3` — not `pose1`, `pose2`)  
  > **Note:** The ligand name must not be a plain integer, but it can contain integers

### PDB Files (`-p`)
One or more protein PDB files, paired to the CSVs by position.

- Used to extract residue names and ordering for the x-axis labels
- The **first PDB is the reference** — all other PDBs are aligned to it by residue-name sequence similarity, so shifted numbering or missing residues across conformations are handled automatically

---

## Flags

### Required

| Flag | Description |
|------|-------------|
| `-c`, `--csv-files` | One or more fingerprint CSV files |
| `-p`, `--pdb-files` | One or more PDB files (must match CSV count and order) |
| `-i`, `--interaction` | Interaction type to plot (see options below) |
| `-g`, `--graph` | Graph type: `heatmap` or `bar` |

**Interaction types** (`-i`):
`contact` · `backbone` · `sidechain` · `polar` · `hydrophobic` · `acceptor` · `donor` · `aromatic` · `charged` · `all`

> `all` runs every interaction type in sequence, producing a separate output file for each.

### Optional

| Flag | Description |
|------|-------------|
| `-ic`, `--ignore-chain` | Remove chain letter from residue labels — not recommended for multi-chain proteins |
| `-s`, `--show` | Display each plot interactively in matplotlib as it's generated |
| `-bar`, `--bar` | Add a bar graph of total interaction propensity above the heatmap |
| `-l`, `--ligands` | Filter to specific ligands by name (space-separated). If omitted, all ligands are plotted |
| `--prefix-ligand-with-source` | Prefix ligand labels with the CSV filename — useful when the same ligand appears in multiple CSVs and you want separate heatmap rows instead of pooled counts |
| `--order-ligands-by-name` | Order heatmap/bar rows alphabetically by ligand name instead of input order |
| `--min-residue-interactions` | Exclude residues with fewer than N total interactions across all ligands |

---

## Usage Examples

**Single CSV/PDB pair — heatmap:**
```bash
python plot_schrodinger_fingerprints.py \
  -c mfsd2b_fingerprint.csv \
  -p mfsd2b.pdb \
  -i contact -g heatmap -ic
```

**Single CSV/PDB pair — bar graph:**
```bash
python plot_schrodinger_fingerprints.py \
  -c mfsd2b_fingerprint.csv \
  -p mfsd2b.pdb \
  -i contact -g bar -ic
```

**Multiple CSV/PDB pairs — pooled heatmap with interaction bar on top:**
```bash
python plot_schrodinger_fingerprints.py \
  -c conf1_fp.csv conf2_fp.csv conf3_fp.csv \
  -p conf1.pdb conf2.pdb conf3.pdb \
  -i contact -g heatmap -bar
```

**Filter to specific ligands, minimum interaction threshold:**
```bash
python plot_schrodinger_fingerprints.py \
  -c protein_fp.csv \
  -p protein.pdb \
  -i polar -g heatmap \
  -l compA compB \
  --min-residue-interactions 5
```

---

## Output

Output files are named automatically:

- **Heatmap:** `{protein}_{interaction}_interaction_heatmap.png`
- **Bar graph:** `{protein}_{ligand}_{interaction}_bargraph.png`

The protein name is taken from the stem of the first CSV file (everything before the first `_`).

---

## Troubleshooting

**`ValueError: The number of FixedLocator locations (...) does not match the number of labels (...)`**  
or labels/axes are getting cut off:

Increase the figure size in the `heatmap` or `bar_graph` function. Find this line and increase the first number (width):

```python
plt.figure(figsize=(12, 8))  # change 12 to a larger value
```
