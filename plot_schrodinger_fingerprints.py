# Written by Marion Q. Lopresti
# Brown Lab 
# Last Update: 4/3/2025

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import os

def load_data(file_path):
    '''
    arg: file_path (str) is path to CSV from user input
    returns: pandas dataframe of the CSV'''
    return pd.read_csv(file_path, dtype=str)

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


def filter_columns(df, interaction):
    '''can handle all the different interactions based on user input
    arg: pandas dataframe from extract_ligand function
    arg: interaction (str) from user input
    returns: list of column names that match the interaction type'''
    return [col for col in df.columns if f'_{interaction}' in col]

def extract_pdb_residues(pdb_file, ignore_chain):
    '''
    arg: pdb_file (str) is path to pdb file from user input
    arg" ignore_chain (bool) is from user input, if -ic flag used then True
    returns a tuple with two dictionaries:
    1. residue_mapping: formatted label based on --ignore-chain
    2. residue_order: just chain ID and residue number to use for sorting residues in graphs later
    '''
    residue_mapping = {}
    residue_order = {}  # store for sorting later because some resn have numbers at start (eg. 2MA = methylated adenosine)

    with open(pdb_file, 'r') as pdb:
        for line in pdb:
            if line.startswith('ATOM') or line.startswith('HETATM'):  
                chain_id = line[21]  
                res_name = line[17:20].strip()  
                res_num = int(line[22:26].strip()) 
                key = f'{chain_id}{res_num}' 
                # to ignore chain in graph
                if ignore_chain:
                    formatted_label = f'{res_name}{res_num}'
                else:
                    formatted_label = f'{chain_id}_{res_name}{res_num}' 
                residue_mapping[key] = formatted_label  
                residue_order[formatted_label] = (chain_id, res_num) # will need chain ID to order the residues in the proper order even if they ignore chain
    return residue_mapping, residue_order

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
            else:
                if pd.notna(value) and str(value).strip() != '':
                    residue_key = col.split('_')[0]
                    formatted_residue = residue_mapping.get(residue_key, residue_key)
                    interaction_counts[(ligand, formatted_residue)] = interaction_counts.get((ligand, formatted_residue), 0) + 1
    return interaction_counts


def create_interaction_dataframe(interaction_counts):
    '''takes interaction counts dictionary and turns it into a dataframe'''
    return pd.DataFrame(
        [(ligand, residue, count) for (ligand, residue), count in interaction_counts.items()],
        columns=['Ligand', 'Residue', 'Count']
    )

def sort_residues(interaction_df, residue_order):
    '''
    arg: interaction_df is the pandas dataframe from create_interaction_dataframe with interaction counts
    arg: residue_order (dict) is from extract_pdb_residues function with chain ID and residue number
    returns: list of sorted residues'''
    heatmap_data = interaction_df.pivot(index='Ligand', columns='Residue', values='Count').fillna(0)
    sorted_residues = sorted(
        heatmap_data.columns, 
        key=lambda label: residue_order.get(label, ('', float('inf')))
    )
    return sorted_residues

def heatmap(interaction_df, sorted_residues, ligands, protein, show):
    '''creates and saves a heatmap from the interaction df with correctly sorted x-ticks
    arg: interaction_df is the dataframe with the interaction counts
    arg: sorted_residues (list) is the residue labels in order
    arg: ligands (list) is the ligand names
    arg: protein (str) has the protein name to label the output files
    arg: show (bool) specified by user and if true it will display a heatmap'''
    heatmap_data = interaction_df.pivot(index='Ligand', columns='Residue', values='Count').fillna(0)
    heatmap_data = heatmap_data.reindex(index=ligands)
    heatmap_data = heatmap_data[sorted_residues]
    plt.figure(figsize=(12, 8))
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
    plt.tight_layout()
    plt.savefig(f'{protein}_interaction_heatmap.png', dpi=300)
    if show:
        plt.show()

def bar_graph(interaction_df, sorted_residues, ligands, protein, show):
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
        plt.savefig(f'{protein}_{ligand}_bargraph.png', dpi=300)
        if show:
            plt.show()

def process_interaction_data(csv_file, pdb_file, interaction, ignore_chain, show, graph):
    '''
    csv_file: csv file path (str)
    pdb_file: pdb_file path (str)
    interaction: interaction type to filter (str)
    ignore_chain: exclude chain IDs when true (bool)
    show: will show plots in matplotlib when true (bool)
    '''
    df = load_data(csv_file)
    df = extract_ligand(df)
    contact_columns = filter_columns(df, interaction)
    filename = os.path.basename(csv_file)
    protein = filename.split('_')[0]  # assumes format: protein_ligand_data.csv
    residue_mapping, residue_order = extract_pdb_residues(pdb_file, ignore_chain)
    interaction_counts = count_interactions(df, contact_columns, residue_mapping)
    if not interaction_counts:
        print(f'No {interaction} interactions found')
        return
    interaction_df = create_interaction_dataframe(interaction_counts)
    sorted_residues = sort_residues(interaction_df, residue_order)
    if graph == 'heatmap':
        heatmap(interaction_df, sorted_residues, df['Ligand'].unique(), protein, show)
    elif graph == 'bar':
        bar_graph(interaction_df, sorted_residues, df['Ligand'].unique(), protein, show)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "This script will generate a heatmap and bar graphs from Schrodinger fingerprint CSV files.\n"
            "Supports multiple ligands per file!\n\n"
            "Example Command:\n"
            "  python plot_schrodinger_fingerprints.py protein_fingerprint.csv protein.pdb -i contact -g bar\n\n"
            "This script assumes that the CSV file is formatted with the protein name first:\n"
            "  {protein}_rest_of_name.csv\n\n"
            "Additionally, the first column in your CSV must be one of the following formats:\n"
            "  1. {protein}_{ligand}_{pose_number}\n"
            "  2. {ligand}_{pose_number}\n"
            "  3. {ligand}\n"
            "  4. {protein}_{ligand}"
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("csv_file", type=str, help="Path to CSV")
    parser.add_argument("pdb_file", type=str, help="Path to PDB file with protein")
    parser.add_argument("-i", "--interaction", type=str, choices=["contact", "backbone", "sidechain", "polar", "hydrophobic", 
                                                                  "acceptor", "donor", "aromatic", "charged"],
        help="Type of interaction to plot")
    parser.add_argument("-g", "--graph", type=str, choices=["bar", "heatmap"], help="Type of graph")
    parser.add_argument('-ic', "--ignore-chain", action="store_true", help="Ignore chain IDs in residue labels.")
    parser.add_argument("-s", "--show", action="store_true", help="Shows the plots as they are plotted in Matplotlib")   
    args = parser.parse_args()
    process_interaction_data(args.csv_file, args.pdb_file, args.interaction, args.ignore_chain, args.show, args.graph)
