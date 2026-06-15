# Written by Marion Q. Lopresti
# Brown Lab 
# Last Update: 4/3/2025

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import re
import dataclasses import dataclass
from pathlib import Path

INTERACTIONS = [
    "contact", "backbone", "sidechain", "polar", "hydrophobic",
    "acceptor", "donor", "aromatic", "charged"
]


@dataclass(frozen=True)
class ResidueRecord:
    chain: str
    resname: str
    resnum: int
    icode: str = ""

    @property
    def csv_key(self) -> str:
        return f"{self.chain}{self.resnum}"

    def label(self, ignore_chain: bool = False) -> str:
        suffix = f"{self.resnum}{self.icode}".strip()
        if ignore_chain:
            return f"{self.resname}{suffix}"
        return f"{self.chain}_{self.resname}{suffix}"

    def order_key(self) -> tuple:
        return (self.chain, self.resnum, self.icode)
        
def load_data(file_path):
    '''
    arg: file_path (str) is path to CSV from user input
    returns: pandas dataframe of the CSV'''
    return pd.read_csv(file_path, low_memory=False)

def extract_ligand(df):
    '''extracts ligand from the first column, can handle these ligand formats in the csv:
    1. {protein}_{ligand}_{pose_number}
    2. {ligand}_{pose_number}
    3. {ligand}
    4. {protein}_{ligand}
    
    arg: pandas dataframe from load_csv function
    returns: df with modified first column to only contain ligand'''
    column_values = df.iloc[:, 0]  
    split_values = column_values.str.split('_')  
    has_pose_number = split_values.str[-1].str.isdigit() 
    ligand_values = split_values.str[-2].where(has_pose_number, split_values.str[-1]) 
    df['Ligand'] = ligand_values 
    return df

def filter_ligands(df: pd.DataFrame, ligands_to_plot: list[str] | None) -> pd.DataFrame:
    if not ligands_to_plot:
        return extract_ligand(df)

    filtered_parts = []

    for ligand in ligands_to_plot:
        ligand = ligand.strip()
        ligand_rows = df[
            df.iloc[:, 0].astype(str).str.contains(
                re.escape(ligand),
                case=False,
                na=False,
            )
        ].copy()

        if ligand_rows.empty:
            print(f"Warning: no rows found for requested ligand: {ligand}")
            continue

        ligand_rows["Ligand"] = ligand
        filtered_parts.append(ligand_rows)

    if not filtered_parts:
        raise ValueError(f"No matching rows found for requested ligands: {ligands_to_plot}")

    return pd.concat(filtered_parts, ignore_index=True)
    
def filter_columns(df, interaction):
    '''can handle all the different interactions based on user input
    arg: pandas dataframe from extract_ligand function
    arg: interaction (str) from user input
    returns: list of column names that match the interaction type'''
    return [col for col in df.columns if f'_{interaction}' in col]

def parse_pdb_residues(pdb_file: str | Path) -> list[ResidueRecord]:
    residues = []
    seen = set()

    with open(pdb_file, "r") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue

            chain = line[21].strip() or "_"
            resname = line[17:20].strip()
            resnum_text = line[22:26].strip()
            icode = line[26].strip()

            if not resnum_text:
                continue

            try:
                resnum = int(resnum_text)
            except ValueError:
                continue

            key = (chain, resnum, icode, resname)
            if key in seen:
                continue

            seen.add(key)
            residues.append(
                ResidueRecord(
                    chain=chain,
                    resname=resname,
                    resnum=resnum,
                    icode=icode,
                )
            )

    return residues


def group_residues_by_chain(residues: list[ResidueRecord]) -> dict[str, list[ResidueRecord]]:
    grouped = {}
    for residue in residues:
        grouped.setdefault(residue.chain, []).append(residue)
    return grouped


def residue_seq_similarity(ref_residues, target_residues) -> float:
    ref_seq = [r.resname for r in ref_residues]
    target_seq = [r.resname for r in target_residues]

    if not ref_seq or not target_seq:
        return 0.0

    return SequenceMatcher(a=ref_seq, b=target_seq, autojunk=False).ratio()


def pair_target_chains_to_reference_chains(
    reference_residues: list[ResidueRecord],
    target_residues: list[ResidueRecord],
) -> dict[str, str]:
    """
    Returns target_chain -> reference_chain.

    First prefers exact chain-ID matches.
    If chain IDs differ, pairs chains by residue-name sequence similarity.
    """
    ref_groups = group_residues_by_chain(reference_residues)
    target_groups = group_residues_by_chain(target_residues)

    chain_map = {}
    used_ref_chains = set()

    # First pass: exact chain ID matches
    for target_chain in target_groups:
        if target_chain in ref_groups:
            chain_map[target_chain] = target_chain
            used_ref_chains.add(target_chain)

    # Second pass: sequence-similarity matching
    for target_chain, target_chain_residues in target_groups.items():
        if target_chain in chain_map:
            continue

        best_ref_chain = None
        best_score = -1.0

        for ref_chain, ref_chain_residues in ref_groups.items():
            if ref_chain in used_ref_chains:
                continue

            score = residue_seq_similarity(ref_chain_residues, target_chain_residues)

            if score > best_score:
                best_score = score
                best_ref_chain = ref_chain

        if best_ref_chain is not None:
            chain_map[target_chain] = best_ref_chain
            used_ref_chains.add(best_ref_chain)

            print(
                f"Info: target chain {target_chain} mapped to reference chain "
                f"{best_ref_chain} by sequence similarity = {best_score:.3f}"
            )

    return chain_map


def align_residue_records(
    reference_residues: list[ResidueRecord],
    target_residues: list[ResidueRecord],
    ignore_chain: bool = False,
) -> dict[str, ResidueRecord]:
    """
    Maps target CSV residue keys, e.g. C135, onto reference ResidueRecord objects.

    Alignment is based on residue-name sequence similarity, not residue number.
    Handles:
    - shifted residue numbering
    - missing crystal residues
    - different chain IDs across structures
    """
    mapping = {}

    if ignore_chain:
        ref_chain_residues = reference_residues
        target_chain_residues = target_residues

        ref_seq = [r.resname for r in ref_chain_residues]
        target_seq = [r.resname for r in target_chain_residues]

        matcher = SequenceMatcher(a=ref_seq, b=target_seq, autojunk=False)

        for block in matcher.get_matching_blocks():
            if block.size == 0:
                continue

            for offset in range(block.size):
                ref_res = ref_chain_residues[block.a + offset]
                target_res = target_chain_residues[block.b + offset]
                mapping[target_res.csv_key] = ref_res

        return mapping

    ref_groups = group_residues_by_chain(reference_residues)
    target_groups = group_residues_by_chain(target_residues)

    target_to_ref_chain = pair_target_chains_to_reference_chains(
        reference_residues=reference_residues,
        target_residues=target_residues,
    )

    for target_chain, ref_chain in target_to_ref_chain.items():
        target_chain_residues = target_groups[target_chain]
        ref_chain_residues = ref_groups[ref_chain]

        ref_seq = [r.resname for r in ref_chain_residues]
        target_seq = [r.resname for r in target_chain_residues]

        matcher = SequenceMatcher(a=ref_seq, b=target_seq, autojunk=False)

        for block in matcher.get_matching_blocks():
            if block.size == 0:
                continue

            for offset in range(block.size):
                ref_res = ref_chain_residues[block.a + offset]
                target_res = target_chain_residues[block.b + offset]

                mapping[target_res.csv_key] = ref_res

    return mapping

def build_reference_maps(
    pdb_files: list[str | Path],
    ignore_chain: bool = False,
) -> tuple[list[str], dict[str, tuple], list[dict[str, str]]]:
    """
    Returns:
    1. reference residue labels in reference PDB order
    2. residue_order dict for sorting
    3. one CSV-key-to-reference-label map per PDB
    """
    all_residue_sets = [parse_pdb_residues(pdb) for pdb in pdb_files]

    reference_residues = all_residue_sets[0]

    reference_labels = [res.label(ignore_chain=ignore_chain) for res in reference_residues]
    residue_order = {
        res.label(ignore_chain=ignore_chain): res.order_key()
        for res in reference_residues
    }

    all_maps = []

    for i, target_residues in enumerate(all_residue_sets):
        target_to_ref = align_residue_records(
            reference_residues=reference_residues,
            target_residues=target_residues,
            ignore_chain=ignore_chain,
        )

        csv_key_to_ref_label = {
            target_key: ref_res.label(ignore_chain=ignore_chain)
            for target_key, ref_res in target_to_ref.items()
        }

        all_maps.append(csv_key_to_ref_label)

        mapped = len(csv_key_to_ref_label)
        total = len(target_residues)
        print(f"PDB {i + 1}: mapped {mapped}/{total} residues onto reference numbering")

    return reference_labels, residue_order, all_maps


def get_column_residue_key(column_name: str) -> str:
    """
    Assumes fingerprint columns start with residue keys such as:
    A135_contact, B42_polar, A135_sidechain, etc.
    """
    return column_name.split("_")[0]


def value_is_interaction(value: object, binary_format: bool) -> bool:
    if pd.isna(value):
        return False

    if binary_format:
        try:
            return float(value) == 1.0
        except ValueError:
            return False

    return str(value).strip() != ""


def detect_binary_format(df: pd.DataFrame, contact_columns: list[str]) -> bool:
    for col in contact_columns:
        nonnull = df[col].dropna()
        if nonnull.empty:
            continue

        first_val = nonnull.iloc[0]
        try:
            return float(first_val) in [0.0, 1.0]
        except ValueError:
            return False

    return False

def count_interactions(df, contact_columns, residue_mapping):
    '''
    Backwards compatible with all versions of Schrodinger
    counts the number of interactions between each {protein}_{ligand} 
    and each {chain ID}{residue_number}, replacing with formatted labels.
    arg: df with interaction data
    arg: list of columns with the correct interaction type
    arg: dictionary of residue with formatted labels
    returns: dictionary with (ligand, residue) keys and interaction counts
    '''
    interaction_counts = {}
    total_residue_counts = {}

    first_col = contact_columns[0]
    first_val = df[first_col].iloc[1]
    binary_format = isinstance(float(first_val), (int, float)) and float(first_val) in [0, 1]
    for _, row in df.iterrows():
        ligand = row['Ligand']
        for col in contact_columns:
            value = row[col]
            if binary_format:
                if value == 1:
                    residue_key = col.split('_')[0]
                    formatted_residue = residue_mapping.get(residue_key, residue_key)
                    interaction_counts[(ligand, formatted_residue)] = interaction_counts.get((ligand, formatted_residue), 0) + 1
                    total_residue_counts[formatted_residue] = total_residue_counts.get(formatted_residue, 0) +1
            else:
                if pd.notna(value) and str(value).strip() != '':
                    residue_key = col.split('_')[0]
                    formatted_residue = residue_mapping.get(residue_key, residue_key)
                    interaction_counts[(ligand, formatted_residue)] = interaction_counts.get((ligand, formatted_residue), 0) + 1
                    total_residue_counts[formatted_residue] = total_residue_counts.get(formatted_residue, 0) +1
    return interaction_counts, total_residue_counts

def merge_count_dicts(
    total_counts: dict,
    new_counts: dict,
) -> dict:
    for key, value in new_counts.items():
        total_counts[key] = total_counts.get(key, 0) + value
    return total_counts
    
def create_interaction_dataframe(interaction_counts):
    '''takes interaction counts dictionary and turns it into a dataframe'''
    return pd.DataFrame(
        [(ligand, residue, count) for (ligand, residue), count in interaction_counts.items()],
        columns=['Ligand', 'Residue', 'Count']
    )

def sort_residues(
    interaction_df: pd.DataFrame,
    residue_order: dict[str, tuple],
    reference_labels: list[str] | None = None,
) -> list[str]:
    if reference_labels is not None:
        present = set(interaction_df["Residue"].dropna().unique())
        return [label for label in reference_labels if label in present]

    heatmap_data = interaction_df.pivot(index="Ligand", columns="Residue", values="Count").fillna(0)

    return sorted(
        heatmap_data.columns,
        key=lambda label: residue_order.get(label, ("", float("inf"), "")),
    )

def heatmap(interaction_df, sorted_residues, ligands, protein, show, total_residue_counts, bar, interaction):
    '''creates and saves a heatmap from the interaction df with correctly sorted x-ticks
    arg: interaction_df is the dataframe with the interaction counts
    arg: sorted_residues (list) is the residue labels in order
    arg: ligands (list) is the ligand names
    arg: protein (str) has the protein name to label the output files
    arg: show (bool) specified by user and if true it will display a heatmap'''
    
    heatmap_data = interaction_df.pivot(index='Ligand', columns='Residue', values='Count').fillna(0)
    heatmap_data = heatmap_data.reindex(index=ligands).fillna(0)
    heatmap_data = heatmap_data.reindex(columns=sorted_residues).fillna(0)
    plt.figure(figsize=(max(12, len(sorted_residues) * 0.35), max(8, len(ligands) * 0.35)))
    ax = sns.heatmap(
        heatmap_data, 
        cmap='Blues', 
        annot=False, 
        fmt='.0f', 
        square=True, 
        cbar_kws={
            'orientation': 'horizontal', 
            'shrink': 0.6, 
            'location': 'bottom', 
            'pad': 0.25
        }
    )
    cbar = ax.collections[0].colorbar
    cbar.outline.set_edgecolor('black')
    cbar.outline.set_linewidth(1.25)
    
    # set black border around heatmap
    for _, spine in ax.spines.items():
        spine.set_visible(True)  
        spine.set_color('black') 
        spine.set_linewidth(1.25)
        
    plt.xticks(rotation=90, fontsize=16)
    ax.set_xticklabels(sorted_residues, fontsize=16)  
    ax.set_yticklabels(ligands, fontsize=16)
    plt.ylabel('Ligand', fontsize=18, fontweight='bold')
    plt.xlabel('Residue', fontsize=18, fontweight='bold')
    # adds bar graph to top of heatmap
    if bar:
        residue_counts = [total_residue_counts.get(res, 0) for res in sorted_residues]
        bar_ax = ax.inset_axes([0, 1.02, 1, 0.25])
        x_pos = range(len(sorted_residues))
        bar_ax.bar(x_pos, residue_counts, color='grey', edgecolor='black', align='center')
        bar_ax.set_xticks(x_pos)
        bar_ax.set_xticklabels(sorted_residues, rotation=90, fontsize=8)
        bar_ax.set_xticks([])
        bar_ax.set_yticks([])
        bar_ax.set_xlim(-0.5, len(sorted_residues)-0.5)

    plt.tight_layout()
    plt.savefig(f'{protein}_{interaction}_interaction_heatmap.png', dpi=300)
    if show:
        plt.show()

def bar_graph(interaction_df, sorted_residues, ligands, protein, show, interaction):
    '''its a bar graph
    arg: interaction_df is the dataframe with the interaction counts
    arg: sorted_residues (list) is the residue labels in order
    arg: ligands (list) is the ligand names
    arg: protein (str) has the protein name to label the output files
    arg: show (bool) specified by user and if true it will display a heatmap'''
    for ligand in ligands:
        ligand_df = interaction_df[interaction_df['Ligand'] == ligand]
        ligand_df = ligand_df.set_index('Residue').reindex(sorted_residues).reset_index()
        muted_blue = sns.color_palette("Blues")[2] 
        plt.figure(figsize=(12, 8))
        ax = sns.barplot(
            x=ligand_df['Residue'], 
            y=ligand_df['Count'], 
            color=muted_blue,
            edgecolor='black'
        )
        plt.xticks(rotation=90, fontsize=16)
        plt.yticks(fontsize=16)
        ax.set_xticklabels(sorted_residues, fontsize=16)  
        plt.ylabel('Interaction Count', fontsize=18, fontweight='bold')
        plt.xlabel('Residue', fontsize=18, fontweight='bold')
        plt.title(f'Interaction Count per Residue - {ligand}', fontsize=20, fontweight='bold')
        for _, spine in ax.spines.items():
            spine.set_visible(True)  
            spine.set_color('black') 
            spine.set_linewidth(1.25)
        plt.tight_layout()
        plt.savefig(f'{protein}_{ligand}_{interaction}_bargraph.png', dpi=300)
        if show:
            plt.show()
def sort_ligands_by_name(ligands: list[str]) -> list[str]:
    def ligand_sort_key(value: str):
        text = str(value)

        if "::" in text:
            source, ligand_part = text.split("::", 1)
        else:
            source, ligand_part = "", text

        try:
            ligand_key = int(ligand_part)
        except ValueError:
            ligand_key = ligand_part

        return (ligand_key, source)

    return sorted(ligands, key=ligand_sort_key)

def process_single_interaction_for_pairs(
    csv_files: list[str],
    pdb_files: list[str],
    interaction: str,
    ignore_chain: bool,
    show: bool,
    graph: str,
    bar: bool,
    ligands_to_plot: list[str] | None = None,
    prefix_ligand_with_source: bool = False,
    order_ligands_by_name: bool = False,
    min_residue_interactions: int = 0,
):
    if len(csv_files) != len(pdb_files):
        raise ValueError(
            "You must provide the same number of CSV files and PDB files. "
            "Each CSV is paired with the PDB at the same position."
        )

    reference_labels, residue_order, residue_maps = build_reference_maps(
        pdb_files=pdb_files,
        ignore_chain=ignore_chain,
    )

    all_interaction_counts = {}
    all_total_residue_counts = {}
    ligand_order = []

    for pair_i, (csv_file, residue_map) in enumerate(zip(csv_files, residue_maps), start=1):
        df = load_data(csv_file)
        df = filter_ligands(df, ligands_to_plot)

        if prefix_ligand_with_source:
            source_label = Path(csv_file).stem
            df["Ligand"] = source_label + "::" + df["Ligand"].astype(str)

        for ligand in df["Ligand"].dropna().unique():
            if ligand not in ligand_order:
                ligand_order.append(ligand)

        contact_columns = filter_columns(df, interaction)

        if not contact_columns:
            print(f"Warning: no columns found for interaction '{interaction}' in {csv_file}")
            continue

        interaction_counts, total_residue_counts = count_interactions(
            df=df,
            contact_columns=contact_columns,
            residue_mapping=residue_map,
        )

        merge_count_dicts(all_interaction_counts, interaction_counts)
        merge_count_dicts(all_total_residue_counts, total_residue_counts)

    if not all_interaction_counts:
        print(f"No {interaction} interactions found")
        return

    interaction_df = create_interaction_dataframe(all_interaction_counts)

    sorted_residues = sort_residues(
        interaction_df=interaction_df,
        residue_order=residue_order,
        reference_labels=reference_labels,
    )
    if min_residue_interactions > 0:
        sorted_residues = [
            res for res in sorted_residues
            if all_total_residue_counts.get(res, 0) >= min_residue_interactions
        ]

        if not sorted_residues:
            print(
                f"No residues passed --min-residue-interactions "
                f"{min_residue_interactions} for {interaction}"
            )
            return
    if order_ligands_by_name:
        ligand_order = sort_ligands_by_name(ligand_order)

    protein = Path(csv_files[0]).stem.split("_")[0]

    if graph == "heatmap":
        heatmap(
            interaction_df=interaction_df,
            sorted_residues=sorted_residues,
            ligands=ligand_order,
            protein=protein,
            show=show,
            total_residue_counts=all_total_residue_counts,
            bar=bar,
            interaction=interaction,
        )

    elif graph == "bar":
        bar_graph(
            interaction_df=interaction_df,
            sorted_residues=sorted_residues,
            ligands=ligand_order,
            protein=protein,
            show=show,
            interaction=interaction,
        )

def process_interaction_data(
    csv_files: list[str],
    pdb_files: list[str],
    interaction: str,
    ignore_chain: bool,
    show: bool,
    graph: str,
    bar: bool,
    ligands_to_plot: list[str] | None = None,
    prefix_ligand_with_source: bool = False,
    order_ligands_by_name: bool = False,
    min_residue_interactions: int = 0,
):
    if interaction == "all" or interaction is None:
        for inter in INTERACTIONS:
            print(f"\nProcessing {inter} interactions")
            process_single_interaction_for_pairs(
                csv_files=csv_files,
                pdb_files=pdb_files,
                interaction=inter,
                ignore_chain=ignore_chain,
                show=show,
                graph=graph,
                bar=bar,
                ligands_to_plot=ligands_to_plot,
                prefix_ligand_with_source=prefix_ligand_with_source,
                order_ligands_by_name=order_ligands_by_name,
                min_residue_interactions=min_residue_interactions,
            )
        return

    process_single_interaction_for_pairs(
        csv_files=csv_files,
        pdb_files=pdb_files,
        interaction=interaction,
        ignore_chain=ignore_chain,
        show=show,
        graph=graph,
        bar=bar,
        ligands_to_plot=ligands_to_plot,
        prefix_ligand_with_source=prefix_ligand_with_source,
        order_ligands_by_name=order_ligands_by_name,
        min_residue_interactions=min_residue_interactions,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Generate heatmaps or bar graphs from Schrodinger fingerprint CSV files.\n"
            "Supports one or multiple CSV/PDB pairs.\n\n"
            "The first CSV/PDB pair defines the reference residue numbering.\n"
            "Additional PDBs are aligned to the first PDB by residue-name sequence similarity.\n\n"
            "Example single pair:\n"
            "  python plot_schrodinger_fingerprints.py -c protein1_fp.csv -p protein1.pdb -i contact -g heatmap\n\n"
            "Example multiple pairs:\n"
            "  python plot_schrodinger_fingerprints.py \\\n"
            "    -c conf1_fp.csv conf2_fp.csv conf3_fp.csv \\\n"
            "    -p conf1.pdb conf2.pdb conf3.pdb \\\n"
            "    -i contact -g heatmap -bar\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument(
        "-c", "--csv-files",
        nargs="+",
        required=True,
        help="One or more fingerprint CSV files. Order must match --pdb-files.",
    )

    parser.add_argument(
        "-p", "--pdb-files",
        nargs="+",
        required=True,
        help="One or more PDB files. Order must match --csv-files. First PDB is the reference.",
    )

    parser.add_argument(
        "-i", "--interaction",
        type=str,
        choices=INTERACTIONS + ["all"],
        default="contact",
        help="Type of interaction to plot.",
    )

    parser.add_argument(
        "-g", "--graph",
        type=str,
        choices=["bar", "heatmap"],
        default="heatmap",
        help="Type of graph.",
    )

    parser.add_argument(
        "-ic", "--ignore-chain",
        action="store_true",
        help="Ignore chain IDs in residue labels.",
    )

    parser.add_argument(
        "-s", "--show",
        action="store_true",
        help="Show plots interactively.",
    )

    parser.add_argument(
        "-bar", "--bar",
        action="store_true",
        help="Adds a bar graph of total interaction propensity above the heatmap.",
    )

    parser.add_argument(
        "-l", "--ligands",
        nargs="+",
        default=None,
        help="Optional ligand names to visualize. If omitted, all ligands are visualized.",
    )

    parser.add_argument(
        "--prefix-ligand-with-source",
        action="store_true",
        help=(
            "Prefix ligand labels with the CSV filename. Use this if the same ligand appears "
            "in multiple CSVs and you want separate heatmap rows instead of pooled counts."
        ),
    )
    parser.add_argument(
        "--order-ligands-by-name",
        action="store_true",
        help="Order heatmap/bar rows by ligand name instead of input protein/CSV order.",
    )
    parser.add_argument(
        "--min-residue-interactions",
        type=int,
        default=0,
        help=(
            "Exclude residues with fewer than this many total interactions "
            "across all plotted ligands/proteins."
        ),
    )
    args = parser.parse_args()

    process_interaction_data(
        csv_files=args.csv_files,
        pdb_files=args.pdb_files,
        interaction=args.interaction,
        ignore_chain=args.ignore_chain,
        show=args.show,
        graph=args.graph,
        bar=args.bar,
        ligands_to_plot=args.ligands,
        prefix_ligand_with_source=args.prefix_ligand_with_source,
        order_ligands_by_name=args.order_ligands_by_name,
        min_residue_interactions=args.min_residue_interactions,
    )
