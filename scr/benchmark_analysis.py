# Script for the analysis of the results and plot generation

######### IMPORTS #########

import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import ttest_rel, pearsonr

########################### SINGLE RNA - FIG 1A ###########################

'''
CSV file names are expected to be named as tm_{method-name}.csv for single RNA TM-scores
and {method}.csv for chain TM-scores
'''

def plot_tm_scores(
    input_dir,
    chain_tm_dir,
    output_dir,
    output_filename="Fig-1A.png"
):
    """
    Generate a TM-score plot and save to file.

    Args:
        input_dir (str): Directory containing original TM-score CSVs.
        chain_tm_dir (str): Directory containing RNA/protein-specific new TM-score CSVs.
        output_dir (str): Directory to save the output plot.
        output_filename (str): Name of the output image file.
    """

    os.makedirs(output_dir, exist_ok=True)

    custom_labels = {
        'tm_af3': 'AF3',
        'tm_chai': 'Chai',
        'tm_hf3': 'HF3',
        'tm_drf': 'DRFold',
        'tm_rf2na': 'RF2NA',
        'tm_boltz': 'Boltz-1',
        'tm_rhofold': 'RhoFold+',
        'tm_nufold': 'NuFold',
        'tm_trrosetta': 'trRosettaRNA'
    }

    # === Load original TM-score data ===
    dfs = []
    csv_files = [f for f in os.listdir(input_dir) if f.endswith('.csv')]
    common_ids = None

    for file in csv_files:
        file_path = os.path.join(input_dir, file)
        df = pd.read_csv(file_path)
        base_name = file.replace('.csv', '')
        model_name = custom_labels.get(base_name, base_name)
        df = df[['PDB_ID', 'TM_Score']].drop_duplicates()
        df['Model'] = model_name
        df['Type'] = 'original'
        dfs.append(df)
        if common_ids is None:
            common_ids = set(df['PDB_ID'])
        else:
            common_ids.intersection_update(set(df['PDB_ID']))

    joint_df = pd.concat([df[df['PDB_ID'].isin(common_ids)] for df in dfs], ignore_index=True)

    # === Load RNA/Protein-specific TM-score data ===
    file_mapping = {
        'af3.csv': 'AF3',
        'boltz.csv': 'Boltz-1',
        'rf2na.csv': 'RF2NA',
        'hf3.csv': 'HF3'
    }

    new_dfs = []
    for filename, model_name in file_mapping.items():
        file_path = os.path.join(chain_tm_dir, filename)
        if not os.path.exists(file_path):
            continue  # skip missing files
        df = pd.read_csv(file_path)
        df = df[['PDB_ID', 'TM_Score', 'Type']].copy()
        df['Model'] = model_name
        new_dfs.append(df)

    combined_df = pd.concat([joint_df] + new_dfs, ignore_index=True)

    # === Map types and define colors ===
    type_mapping = {
        'rna': 'RNA Chains',
        'protein': 'Protein Chains',
        'original': 'RNA'
    }

    color_palette = {
        'RNA Chains': '#ee82ee',
        'RNA': 'skyblue',
        'Protein Chains': '#1f77b4'
    }

    combined_df['Type'] = combined_df['Type'].map(type_mapping)

    # === Get model order ===
    models_with_data = combined_df['Model'].unique()
    original_means = (
        combined_df[combined_df['Type'] == 'RNA']
        .groupby('Model')['TM_Score']
        .mean()
        .sort_values(ascending=False)
    )
    model_order = [m for m in original_means.index if m in models_with_data]

    # === Manual plotting with fixed offsets ===
    fig, ax = plt.subplots(figsize=(2 * len(model_order), 5))

    x_ticks = []
    x_labels = []
    legend_handles = {}
    box_width = 0.30

    type_offsets = {
        'RNA Chains': -box_width,
        'RNA': 0,
        'Protein Chains': box_width
    }

    for i, model in enumerate(model_order):
        subset = combined_df[combined_df['Model'] == model]

        for t in ['RNA', 'RNA Chains', 'Protein Chains']:
            scores = subset[subset['Type'] == t]['TM_Score'].dropna()
            if scores.empty:
                continue

            pos = i + type_offsets[t]
            bplot = ax.boxplot(
                scores,
                positions=[pos],
                widths=box_width,
                patch_artist=True,
                boxprops=dict(facecolor=color_palette[t], edgecolor='black', linewidth=1),
                medianprops=dict(color='black'),
                whiskerprops=dict(color='black'),
                capprops=dict(color='black'),
                flierprops=dict(marker='o', markersize=3, linestyle='none', markerfacecolor=color_palette[t])
            )

            if t not in legend_handles:
                legend_handles[t] = bplot["boxes"][0]

        x_ticks.append(i)
        x_labels.append(model)

    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_labels, rotation=0, ha='center', fontsize=12, fontweight='bold')
    ax.set_ylabel("TM-Score", fontsize=12)
    ax.set_title("Performance on Single RNA", fontsize=14)
    ax.set_ylim(0, 1)
    ax.legend(legend_handles.values(), legend_handles.keys(), bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
    plt.tight_layout()

    output_path = os.path.join(output_dir, output_filename)
    plt.savefig(output_path)
    plt.close()






########################### SINGLE RNA - FIG 1B ###########################

def plot_success_rates(
    input_dir,
    chain_tm_dir,
    output_dir,
    plot_filename="Fig-1B.png",
    summary_filename="Table-S2.csv"
):
    """
    Generate a success rate bar plot and a summary CSV file.

    Args:
        input_dir (str): Directory containing original TM-score CSVs.
        chain_tm_dir (str): Directory containing RNA/protein-specific TM-score CSVs.
        output_dir (str): Directory to save outputs.
        plot_filename (str): Name of the output plot image file.
        summary_filename (str): Name of the output summary CSV file.
    """

    os.makedirs(output_dir, exist_ok=True)

    custom_labels = {
        'tm_af3': 'AF3',
        'tm_chai': 'Chai  ',
        'tm_hf3': 'HF3',
        'tm_drf': 'DRFold',
        'tm_rf2na': 'RF2NA',
        'tm_boltz': 'Boltz-1',
        'tm_rhofold': 'RhoFold+',
        'tm_nufold': 'NuFold',
        'tm_trrosetta': 'trRosettaRNA'
    }

    # === Load original TM-score data ===
    dfs = []
    csv_files = [f for f in os.listdir(input_dir) if f.endswith('.csv')]
    common_ids = None

    for file in csv_files:
        file_path = os.path.join(input_dir, file)
        df = pd.read_csv(file_path)
        base_name = file.replace('.csv', '')
        model_name = custom_labels.get(base_name, base_name)
        df = df[['PDB_ID', 'TM_Score']].drop_duplicates()
        df.rename(columns={'TM_Score': model_name}, inplace=True)
        dfs.append(df)
        if common_ids is None:
            common_ids = set(df['PDB_ID'])
        else:
            common_ids.intersection_update(set(df['PDB_ID']))

    joint_df = None
    for df in dfs:
        df = df[df['PDB_ID'].isin(common_ids)]
        joint_df = df if joint_df is None else pd.merge(joint_df, df, on='PDB_ID', how='inner')

    mean_tmscores = joint_df.drop(columns=['PDB_ID']).mean().sort_values(ascending=False)
    model_order = mean_tmscores.index.tolist()

    success_rates = (joint_df.drop(columns=['PDB_ID']) > 0.5).mean() * 100
    success_df = success_rates.reset_index()
    success_df.columns = ['Model', 'Value']
    success_df['Metric'] = 'Success Rate'

    # === Load RNA/Protein-specific success rates ===
    file_mapping = {
        'af3.csv': 'AF3',
        'boltz.csv': 'Boltz-1',
        'rf2na.csv': 'RF2NA',
        'hf3.csv': 'HF3'
    }

    rna_records = []
    protein_records = []
    for filename, model in file_mapping.items():
        file_path = os.path.join(chain_tm_dir, filename)
        if not os.path.exists(file_path):
            continue
        df = pd.read_csv(file_path)
        df = df[['PDB_ID', 'TM_Score', 'Type']]
        rna_success = (df[df['Type'] == 'rna']['TM_Score'] > 0.5).mean() * 100
        protein_success = (df[df['Type'] == 'protein']['TM_Score'] > 0.5).mean() * 100
        rna_records.append({'Model': model, 'Value': rna_success, 'Metric': 'RNA Success Rate'})
        protein_records.append({'Model': model, 'Value': protein_success, 'Metric': 'Protein Success Rate'})

    all_data = pd.concat([
        success_df,
        pd.DataFrame(rna_records),
        pd.DataFrame(protein_records)
    ], ignore_index=True)

    # === Prepare plotting ===
    metric_order = ['RNA Chains', 'RNA', 'Protein Chain']
    metric_colors = {
        'RNA Chains': '#ee82ee',
        'RNA': 'skyblue',
        'Protein Chain': '#1f77b4'
    }
    metric_rename = {
        'RNA Success Rate': 'RNA Chains',
        'Success Rate': 'RNA',
        'Protein Success Rate': 'Protein Chain'
    }
    all_data['Metric'] = all_data['Metric'].map(metric_rename)

    fig, ax = plt.subplots(figsize=(2.3 * len(model_order), 6))

    bar_width = 3
    bar_spacing = 0.025
    group_spacing = 3.6

    x_locs = []
    x_labels = []
    legend_handles = {}
    current_x = 0

    for model in model_order:
        model_data = all_data[all_data['Model'] == model]
        present_metrics = [metric for metric in metric_order if not model_data[model_data['Metric'] == metric].empty]

        group_start_x = current_x
        for metric in present_metrics:
            value = model_data[model_data['Metric'] == metric]['Value'].values[0]
            bar = ax.bar(
                current_x,
                value,
                width=bar_width,
                color=metric_colors[metric],
                label=metric if metric not in legend_handles else "",
                zorder=3
            )
            if metric not in legend_handles:
                legend_handles[metric] = bar[0]

            ax.text(
                current_x,
                value + 2,
                f'{value:.1f}%',
                ha='center',
                va='bottom',
                fontsize=13
            )

            current_x += bar_width + bar_spacing

        group_width = len(present_metrics) * (bar_width + bar_spacing) - bar_spacing
        x_locs.append(group_start_x + group_width / 2)
        x_labels.append(model)
        current_x += group_spacing

    ax.set_xticks(x_locs)
    ax.set_xticklabels(x_labels, rotation=0, ha='right', fontsize=13, fontweight='bold')
    ax.set_ylim(0, 100)
    ax.set_ylabel('Success Rate [%]', fontsize=15)
    ax.set_title('Single RNA Success Rates', fontsize=17)
    ax.legend(legend_handles.values(), legend_handles.keys(), title='Type', bbox_to_anchor=(1.02, 1), loc='upper left')
    ax.grid(False)

    plt.tight_layout()
    plot_path = os.path.join(output_dir, plot_filename)
    plt.savefig(plot_path)
    plt.close()

    # === Generate and save summary statistics ===
    summary_records = []

    for model in model_order:
        model_scores = joint_df[model]
        n_models = len(model_scores)
        mean_tm = model_scores.mean()
        success = (model_scores > 0.5).mean() * 100
        summary_records.append({
            'Model': model,
            'Type': 'Overall',
            'N Models': n_models,
            'Mean TM-Score': mean_tm,
            'Success %': success
        })

    for filename, model in file_mapping.items():
        file_path = os.path.join(chain_tm_dir, filename)
        if not os.path.exists(file_path):
            continue
        df = pd.read_csv(file_path)
        df = df[['PDB_ID', 'TM_Score', 'Type']]

        for data_type in ['rna', 'protein']:
            sub_df = df[df['Type'] == data_type]
            n_models = len(sub_df)
            mean_tm = sub_df['TM_Score'].mean()
            success = (sub_df['TM_Score'] > 0.5).mean() * 100
            summary_records.append({
                'Model': model,
                'Type': data_type.capitalize(),
                'N Models': n_models,
                'Mean TM-Score': mean_tm,
                'Success %': success
            })

    summary_df = pd.DataFrame(summary_records)
    summary_path = os.path.join(output_dir, summary_filename)
    summary_df.to_csv(summary_path, index=False)



########################### SINGLE RNA - FIG 1C, S1 - SCATTER PLOTS ###########################


def compare_to_boltz(
    input_dir,
    output_dir
):
    """
    Compare models to Boltz-1, generate and save scatter plots.

    Args:
        input_dir (str): Directory containing TM-score CSVs for models.
        output_dir (str): Directory to save output plots.
    """

    os.makedirs(output_dir, exist_ok=True)

    custom_labels = {
        'tm_af3': 'AF3',
        'tm_chai': 'Chai',
        'tm_hf3': 'HF3',
        'tm_drf': 'DRFold',
        'tm_rf2na': 'RF2NA',
        'tm_boltz': 'Boltz-1',
        'tm_rhofold': 'RhoFold+',
        'tm_nufold': 'NuFold',
        'tm_trrosetta': 'trRosettaRNA'
    }

    # Load all CSVs
    csv_files = [f for f in os.listdir(input_dir) if f.endswith('.csv')]
    all_data = {}
    rnatype_info = {}
    length_info = {}

    for csv_file in csv_files:
        label_key = csv_file.split('.')[0]
        if label_key not in custom_labels:
            continue
        label = custom_labels[label_key]
        full_path = os.path.join(input_dir, csv_file)
        df = pd.read_csv(full_path)
        df['Model'] = label

        if label == 'Boltz-1':
            if 'rnatype' not in df.columns or 'length' not in df.columns:
                raise ValueError("Expected 'rnatype' and 'length' columns in Boltz-1 CSV.")
            df['rnatype'] = df['rnatype'].replace({'TetraLoop': 'Loop', 'IntraLoop': 'Loop'})
            rnatype_info = df.set_index('PDB_ID')['rnatype'].to_dict()
            length_info = df.set_index('PDB_ID')['length'].to_dict()

        all_data[label] = df[['PDB_ID', 'TM_Score']].rename(columns={'TM_Score': label})

    if 'Boltz-1' not in all_data:
        raise ValueError("Boltz-1 data not found in the provided directory.")

    boltz_df = all_data['Boltz-1']
    boltz_pdbs = set(boltz_df['PDB_ID'])

    for method, df in all_data.items():
        if method == 'Boltz-1':
            continue

        pdbs = boltz_pdbs & set(df['PDB_ID'])
        if not pdbs:
            print(f"Skipping {method}: no common PDBs found.")
            continue

        merged_df = pd.merge(
            df[df['PDB_ID'].isin(pdbs)],
            boltz_df[boltz_df['PDB_ID'].isin(pdbs)],
            on='PDB_ID'
        )

        merged_df = merged_df.dropna()
        merged_df = merged_df.rename(columns={method: 'Method', 'Boltz-1': 'Boltz'})
        merged_df['rnatype'] = merged_df['PDB_ID'].map(rnatype_info)

        merged_df = merged_df.sort_values(by='Method').reset_index(drop=True)
        merged_df['moving_avg_Method'] = merged_df['Method'].rolling(window=10).mean()
        merged_df['moving_avg_Boltz'] = merged_df['Boltz'].rolling(window=10).mean()

        t_stat, p_value = ttest_rel(merged_df['Method'], merged_df['Boltz'])
        corr = np.corrcoef(merged_df['Method'], merged_df['Boltz'])[0, 1]
        mean_method = merged_df['Method'].mean()
        mean_boltz = merged_df['Boltz'].mean()

        # Plot
        custom_palette = {
            'Loop': '#64bc3e',
            'L-shaped': '#2485dc',
            'Quadruplex': '#ffa500',
            'Complex': '#8d70e6'
        }

        sns.set(style="white")
        plot = sns.jointplot(
            data=merged_df,
            x='Method',
            y='Boltz',
            kind='scatter',
            hue='rnatype',
            palette=custom_palette,
            height=7,
            space=0,
            ratio=6
        )

        maxval = max(merged_df['Method'].max(), merged_df['Boltz'].max())
        maxval = 1 if maxval < 1 else 100 if maxval > 30 else maxval

        sns.lineplot(x=[0, maxval], y=[0, maxval], ax=plot.ax_joint, linestyle='--', color='red')
        sns.lineplot(
            data=merged_df,
            x='moving_avg_Method',
            y='moving_avg_Boltz',
            ax=plot.ax_joint,
            color='purple'
        )

        cutoff = maxval / 10
        for idx, row in merged_df.iterrows():
            if abs(row['Method'] - row['Boltz']) > cutoff:
                plot.ax_joint.text(row['Method'], row['Boltz'], row['PDB_ID'], fontsize=8, alpha=0.7)

        plot.ax_joint.text(maxval * 0.001, maxval * 0.95,
            f'Cc: {corr:.2f}\np: {p_value:.4e}\n<{method}>: {mean_method:.2f}\n<Boltz>: {mean_boltz:.2f}',
            fontsize=8,
        )

        plot.ax_joint.set_xlim(-0.1, 1.15)
        plot.ax_joint.set_ylim(-0.1, 1.15)
        plot.ax_joint.set_xlabel(f'TM-score ({method})', fontsize=12)
        plot.ax_joint.set_ylabel('TM-score (Boltz-1)', fontsize=12)
        plot.ax_joint.legend(loc='lower right', title='Type', fontsize=9, title_fontsize=10)

        plot.figure.suptitle(f'Boltz-1 vs {method} TM-Score', fontsize=14)
        plot.figure.subplots_adjust(top=0.95)

        # Save plot
        if method == 'AF3':
            fig_name = "Fig-1C.png"
        else:
            fig_name = f"Fig-S6-boltz-{method}.png"

        plot_path = os.path.join(output_dir, fig_name)
        plot.figure.savefig(plot_path)
        plt.close(plot.figure)


########################### SINGLE RNA - FIG 1D - CHAIN SCATTER PLOTS ###########################

def plot_tm_chains(chain_tm_dir, output_dir):
    # File paths
    file_paths = {
        'AF3': os.path.join(chain_tm_dir, 'af3.csv'),
        'Boltz-1': os.path.join(chain_tm_dir, 'boltz.csv'),
        'RF2NA': os.path.join(chain_tm_dir, 'rf2na.csv'),
        'HF3': os.path.join(chain_tm_dir, 'hf3.csv')
    }

    # Load data
    all_data = {}
    type_info = {}

    for model, path in file_paths.items():
        df = pd.read_csv(path)
        if 'Type' not in df.columns:
            raise ValueError(f"'Type' column missing in {path}")
        if 'PDB_ID' not in df.columns or 'TM_Score' not in df.columns:
            raise ValueError(f"'PDB_ID' or 'TM_Score' column missing in {path}")

        all_data[model] = df[['PDB_ID', 'TM_Score']].rename(columns={'TM_Score': model})

        if model == 'Boltz-1':
            type_info = df.set_index('PDB_ID')['Type'].to_dict()

    if 'Boltz-1' not in all_data:
        raise ValueError("Boltz-1 data not found.")

    boltz_df = all_data['Boltz-1']
    boltz_pdbs = set(boltz_df['PDB_ID'])

    custom_palette = {
        'RNA Chain': '#ee82ee',
        'Protein Chain': '#1f77b4'
    }

    os.makedirs(output_dir, exist_ok=True)  # Make sure output_dir exists

    # Compare each method with Boltz-1
    for method, df in all_data.items():
        if method == 'Boltz-1':
            continue

        pdbs = boltz_pdbs & set(df['PDB_ID'])
        if not pdbs:
            print(f"Skipping {method}: no common PDBs.")
            continue

        merged_df = pd.merge(
            df[df['PDB_ID'].isin(pdbs)],
            boltz_df[boltz_df['PDB_ID'].isin(pdbs)],
            on='PDB_ID'
        )

        merged_df = merged_df.dropna()
        merged_df = merged_df.rename(columns={method: 'Method', 'Boltz-1': 'Boltz'})
        merged_df['Type'] = merged_df['PDB_ID'].map(type_info)
        merged_df['Type'] = merged_df['Type'].replace({'rna': 'RNA Chain', 'protein': 'Protein Chain'})

        merged_df = merged_df.sort_values(by='Method').reset_index(drop=True)
        merged_df['moving_avg_Method'] = merged_df['Method'].rolling(window=10).mean()
        merged_df['moving_avg_Boltz'] = merged_df['Boltz'].rolling(window=10).mean()

        t_stat, p_value = ttest_rel(merged_df['Method'], merged_df['Boltz'])
        corr = np.corrcoef(merged_df['Method'], merged_df['Boltz'])[0, 1]
        mean_method = merged_df['Method'].mean()
        mean_boltz = merged_df['Boltz'].mean()

        sns.set(style="white")
        plot = sns.jointplot(
            data=merged_df,
            x='Method',
            y='Boltz',
            kind='scatter',
            hue='Type',
            palette=custom_palette,
            height=7,
            space=0,
            ratio=6,
            marginal_kws={"bw_adjust": 0.5}
        )

        maxval = max(merged_df['Method'].max(), merged_df['Boltz'].max())
        maxval = 1 if maxval < 1 else 100 if maxval > 30 else maxval

        sns.lineplot(x=[0, maxval], y=[0, maxval], ax=plot.ax_joint, linestyle='--', color='red')
        sns.lineplot(
            data=merged_df,
            x='moving_avg_Method',
            y='moving_avg_Boltz',
            ax=plot.ax_joint,
            color='purple'
        )

        cutoff = maxval / 10
        for idx, row in merged_df.iterrows():
            if abs(row['Method'] - row['Boltz']) > cutoff:
                if '_chain' in row['PDB_ID']:
                    base, chain = row['PDB_ID'].split('_chain')
                    pdb_clean = f"{base}-{chain}"
                else:
                    pdb_clean = row['PDB_ID']
                plot.ax_joint.text(row['Method'], row['Boltz'], pdb_clean, fontsize=8, alpha=0.7)

        plot.ax_joint.text(maxval * 0.001, maxval * 0.95,
            f'Cc: {corr:.2f}\np: {p_value:.4e}\n<{method}>: {mean_method:.2f}\n<Boltz>: {mean_boltz:.2f}',
            fontsize=8)

        plot.ax_joint.set_xlim(-0.1, 1.15)
        plot.ax_joint.set_ylim(-0.1, 1.15)
        plot.ax_joint.legend(loc='lower right', title='Type', fontsize=9, title_fontsize=10)

        plot.figure.suptitle(f'Boltz-1 vs {method} Chain TM-Score', fontsize=14)
        plot.figure.subplots_adjust(top=0.95)
        plot.ax_joint.set_xlabel(f'TM-score ({method})', fontsize=12)
        plot.ax_joint.set_ylabel('TM-score (Boltz-1)', fontsize=12)

        # ---- Save Figures ----
        if method == 'AF3':
            save_path = os.path.join(output_dir, 'Fig-1D.png')
        else:
            save_path = os.path.join(output_dir, f'Fig-chainTM-boltz-{method}.png')

        plot.figure.savefig(save_path)
        plt.close(plot.figure)


########################### SINGLE RNA - FIG S2B - RNA LENGTH SCATTER PLOTS ###########################

def rna_length_scatter(tm_score_dir, output_dir):
    custom_labels = {
        'tm_af3': 'AF3',
        'tm_chai': 'Chai',
        'tm_hf3': 'HF3',
        'tm_drf': 'DRFold',
        'tm_rf2na': 'RF2NA',
        'tm_boltz': 'Boltz-1',
        'tm_rhofold': 'RhoFold+',
        'tm_nufold': 'NuFold',
        'tm_trrosetta': 'trRosettaRNA'
    }

    os.makedirs(output_dir, exist_ok=True)

    csv_files = [f for f in os.listdir(tm_score_dir) if f.endswith('.csv')]
    all_data = {}
    rnatype_info = {}
    length_info = {}

    for csv_file in csv_files:
        label_key = csv_file.split('.')[0]
        if label_key not in custom_labels:
            continue
        label = custom_labels[label_key]
        full_path = os.path.join(tm_score_dir, csv_file)
        df = pd.read_csv(full_path)
        df['Model'] = label

        if label == 'Boltz-1':
            if 'rnatype' not in df.columns or 'length' not in df.columns:
                raise ValueError("Expected 'rnatype' and 'length' columns in Boltz-1 CSV.")
            df['rnatype'] = df['rnatype'].replace({'TetraLoop': 'Loop', 'IntraLoop': 'Loop'})
            rnatype_info = df.set_index('PDB_ID')['rnatype'].to_dict()
            length_info = df.set_index('PDB_ID')['length'].to_dict()

        all_data[label] = df[['PDB_ID', 'TM_Score']].rename(columns={'TM_Score': label})

    if 'Boltz-1' not in all_data:
        raise ValueError("Boltz data not found.")

    boltz_df = all_data['Boltz-1']
    boltz_pdbs = set(boltz_df['PDB_ID'])

    for method, df in all_data.items():
        if method == 'Boltz-1':
            continue

        pdbs = boltz_pdbs & set(df['PDB_ID'])
        if not pdbs:
            print(f"Skipping {method}: no common PDBs.")
            continue

        merged_df = pd.merge(
            df[df['PDB_ID'].isin(pdbs)],
            boltz_df[boltz_df['PDB_ID'].isin(pdbs)],
            on='PDB_ID'
        )

        merged_df = merged_df.dropna()
        merged_df = merged_df.rename(columns={method: 'Method', 'Boltz-1': 'Boltz'})
        merged_df['length'] = merged_df['PDB_ID'].map(length_info)

        merged_df = merged_df.sort_values(by='Method').reset_index(drop=True)
        merged_df['moving_avg_Method'] = merged_df['Method'].rolling(window=10).mean()
        merged_df['moving_avg_Boltz'] = merged_df['Boltz'].rolling(window=10).mean()

        t_stat, p_value = ttest_rel(merged_df['Method'], merged_df['Boltz'])
        corr = np.corrcoef(merged_df['Method'], merged_df['Boltz'])[0, 1]
        mean_method = merged_df['Method'].mean()
        mean_boltz = merged_df['Boltz'].mean()

        fig, ax = plt.subplots(figsize=(7, 7))

        lengths = merged_df['length']
        scatter = sns.scatterplot(
            data=merged_df,
            x='Method',
            y='Boltz',
            hue=lengths,
            palette='viridis',
            edgecolor='none',
            s=40,
            legend=False,
            ax=ax
        )

        norm = plt.Normalize(lengths.min(), lengths.max())
        sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax)
        cbar.set_label('RNA Length [nt]', fontsize=10)

        maxval = max(merged_df['Method'].max(), merged_df['Boltz'].max())
        maxval = 1 if maxval < 1 else 100 if maxval > 30 else maxval
        ax.plot([0, maxval], [0, maxval], linestyle='--', color='red')

        ax.plot(
            merged_df['moving_avg_Method'],
            merged_df['moving_avg_Boltz'],
            color='purple'
        )

        cutoff = maxval / 10
        for idx, row in merged_df.iterrows():
            if abs(row['Method'] - row['Boltz']) > cutoff:
                ax.text(row['Method'], row['Boltz'], row['PDB_ID'], fontsize=8, alpha=0.7)

        ax.text(maxval * 0.001, maxval * 0.99,
            f'Cc: {corr:.2f}\np: {p_value:.4e}\n<{method}>: {mean_method:.2f}\n<Boltz>: {mean_boltz:.2f}',
            fontsize=8, verticalalignment='top')

        ax.set_xlim(-0.1, 1.15)
        ax.set_ylim(-0.1, 1.15)

        ax.set_xlabel(f'TM-score ({method})', fontsize=12)
        ax.set_ylabel('TM-score (Boltz-1)', fontsize=12)
        ax.set_title(f'Boltz-1 vs {method} TM-Score', fontsize=14)

        plt.tight_layout()

        # ---- Save Figure ----
        if method == 'AF3':
            save_path = os.path.join(output_dir, 'Fig-S7B.png')
        else:
            save_path = os.path.join(output_dir, f'Fig-length-boltz-{method}.png')

        fig.savefig(save_path)
        plt.close(fig)




########################### RNA COMPLEXES - FIG 3A, 3B ###########################

def plot_tm_dockq_complexes(
    input_dir,
    output_dir
):
    """
    Generate TM-Score and DockQ violin plots, and save summary statistics.

    Args:
        input_dir (str): Directory containing input CSV files.
        output_dir (str): Directory to save output plots and summary CSV.
    """

    os.makedirs(output_dir, exist_ok=True)

    custom_labels = {
        'usf1_af3': 'AF3',
        'usf1_hf3': 'HF3',
        'usf1_rf2na': 'RF2NA',
        'usf1_boltz': 'Boltz'
    }

    # Load TM-Score Data
    tm_dict = {}
    dockq_dict = {}

    csv_files = [f for f in os.listdir(input_dir) if f.endswith('.csv')]

    for file in csv_files:
        file_path = os.path.join(input_dir, file)
        df = pd.read_csv(file_path)
        base_name = file.replace('.csv', '')
        model_name = custom_labels.get(base_name, base_name)

        tm_df = df[['PDB_ID', 'TM_Score']].copy()
        tm_df['Model'] = model_name
        tm_df = tm_df.drop_duplicates(subset='PDB_ID')
        tm_dict[base_name] = tm_df

        dockq_df = df[['PDB_ID', 'Total_DockQ_Score']].rename(columns={'Total_DockQ_Score': 'DockQ_Score'}).copy()
        dockq_df['Model'] = model_name
        dockq_dict[base_name] = dockq_df

    # Find common PDB IDs
    common_ids = set(tm_dict['usf1_af3']['PDB_ID'])
    for key in tm_dict.keys():
        common_ids = common_ids.intersection(set(tm_dict[key]['PDB_ID']))

    for key in tm_dict.keys():
        tm_dict[key] = tm_dict[key][tm_dict[key]['PDB_ID'].isin(common_ids)]
        dockq_dict[key] = dockq_dict[key][dockq_dict[key]['PDB_ID'].isin(common_ids)]

    tm_combined = pd.concat(tm_dict.values(), ignore_index=True)
    dockq_combined = pd.concat(dockq_dict.values(), ignore_index=True)

    num_common_ids = len(common_ids)

    # === Plot TM-Score Violin ===
    model_tm_means = tm_combined.groupby('Model')['TM_Score'].mean().sort_values(ascending=False)
    sorted_models_tm = model_tm_means.index.tolist()

    plt.figure(figsize=(13, 7))
    sns.violinplot(
        x='Model',
        y='TM_Score',
        data=tm_combined,
        width=0.8,
        order=sorted_models_tm,
        cut=0,
        scale='width',
        color='mediumseagreen'
    )
    plt.ylabel('USAlign TM-Score', fontsize=18, fontweight='bold')
    plt.xlabel('', fontsize=15, fontweight='bold')
    plt.title('RNA/RNA/Protein Complex Prediction Accuracy', fontsize=16)
    plt.xticks(fontsize=17, fontweight='bold', rotation=0)
    plt.tight_layout()
    tm_plot_path = os.path.join(output_dir, "Fig-3A.png")
    plt.savefig(tm_plot_path)
    plt.close()

    # === Plot DockQ Violin ===
    model_dockq_means = dockq_combined.groupby('Model')['DockQ_Score'].mean().sort_values(ascending=False)
    sorted_models_dockq = model_dockq_means.index.tolist()

    plt.figure(figsize=(14, 7))
    sns.violinplot(
        x='Model',
        y='DockQ_Score',
        data=dockq_combined,
        width=0.8,
        order=sorted_models_dockq,
        cut=0,
        scale='width',
        color='mediumseagreen'
    )
    plt.ylabel('DockQ Score', fontsize=18, fontweight='bold')
    plt.xlabel('')
    plt.title('RNA/RNA/Protein Complex DockQ', fontsize=16)
    plt.xticks(fontsize=17, fontweight='bold', rotation=0)
    plt.tight_layout()
    dockq_plot_path = os.path.join(output_dir, "Fig-3B.png")
    plt.savefig(dockq_plot_path)
    plt.close()

    # === Prepare Summary CSV ===
    summary_records = []

    for model in sorted_models_tm:  # ensure TM-score model ordering
        tm_scores = tm_combined[tm_combined['Model'] == model]
        dockq_scores = dockq_combined[dockq_combined['Model'] == model]

        mean_tm = tm_scores['TM_Score'].mean()
        mean_dockq = dockq_scores['DockQ_Score'].mean()
        count = tm_scores['PDB_ID'].nunique()

        high_quality_ids = dockq_scores[dockq_scores['DockQ_Score'] > 0.23]['PDB_ID'].nunique()
        success_rate = (high_quality_ids / num_common_ids) * 100

        summary_records.append({
            'Method': model,
            'Mean TM_Score': round(mean_tm, 3),
            'Mean DockQ_Score': round(mean_dockq, 3),
            'Count': count,
            'Success_rate_dockq': round(success_rate, 1)
        })

    summary_df = pd.DataFrame(summary_records)
    summary_csv_path = os.path.join(output_dir, "Table-S5.csv")
    summary_df.to_csv(summary_csv_path, index=False)




########################### RNA COMPLEXES - FIG 3D, 3E ###########################

def plot_tm_dockq_scatter(
    input_dir,
    output_dir
):
    """
    Create scatter plots comparing AF3 to other models for TM-Score and DockQ Score.

    Args:
        input_dir (str): Directory containing input CSV files.
        output_dir (str): Directory to save output plots.
    """

    os.makedirs(output_dir, exist_ok=True)

    custom_labels = {
        'usf1_af3': 'AF3',
        'usf1_hf3': 'HF3',
        'usf1_rf2na': 'RF2NA',
        'usf1_boltz': 'Boltz'
    }

    csv_files = [f for f in os.listdir(input_dir) if f.endswith('.csv')]

    dfs = {}
    common_pdb_ids = None

    for csv_file in csv_files:
        label_key = csv_file.split('.')[0]
        if label_key not in custom_labels:
            continue
        label = custom_labels[label_key]
        df = pd.read_csv(os.path.join(input_dir, csv_file))
        df['Model'] = label
        dfs[label] = df

        if common_pdb_ids is None:
            common_pdb_ids = set(df['PDB_ID'])
        else:
            common_pdb_ids &= set(df['PDB_ID'])

    if 'AF3' not in dfs:
        raise ValueError("AF3 model CSV is required for baseline comparisons.")

    if not common_pdb_ids:
        raise ValueError("No common PDB IDs found among models.")

    # Filter to common IDs
    for label in dfs:
        dfs[label] = dfs[label][dfs[label]['PDB_ID'].isin(common_pdb_ids)].drop_duplicates(subset=['PDB_ID'])

    combined_data = pd.concat(dfs.values(), ignore_index=True)

    # ===========================
    # TM-SCORE SCATTER PLOTS
    # ===========================

    tm_pivot = combined_data.pivot(index='PDB_ID', columns='Model', values='TM_Score').reset_index()
    methods = [label for label in custom_labels.values() if label != 'AF3']

    for method in methods:
        if method not in tm_pivot.columns:
            print(f"Skipping {method} for TM-score, not found.")
            continue

        plot_df = tm_pivot.dropna(subset=['AF3', method]).copy()
        plot_df = plot_df.sort_values(by='AF3').reset_index(drop=True)
        plot_df['moving_avg_AF3'] = plot_df['AF3'].rolling(window=10).mean()
        plot_df['moving_avg_' + method] = plot_df[method].rolling(window=10).mean()

        t_stat, p_value = ttest_rel(plot_df['AF3'], plot_df[method])

        sns.set(style="white")
        plot = sns.jointplot(
            data=plot_df,
            x='AF3',
            y=method,
            kind='scatter',
            height=7,
            space=0
        )

        maxval = max(plot_df['AF3'].max(), plot_df[method].max())
        maxval = 1 if maxval < 1 else 100 if maxval > 30 else maxval
        sns.lineplot(x=[0, maxval], y=[0, maxval], ax=plot.ax_joint, linestyle='--', color='red')

        sns.lineplot(
            data=plot_df,
            x='moving_avg_AF3',
            y='moving_avg_' + method,
            ax=plot.ax_joint,
            color='purple'
        )

        cutoff = maxval / 10
        for idx, row in plot_df.iterrows():
            if abs(row['AF3'] - row[method]) > cutoff:
                plot.ax_joint.text(row['AF3'], row[method], row['PDB_ID'], fontsize=8, alpha=0.7)

        corr = np.corrcoef(plot_df['AF3'], plot_df[method])[0, 1]
        mean_af3 = plot_df['AF3'].mean()
        mean_method = plot_df[method].mean()

        plot.ax_joint.text(maxval * 0.05, maxval * 0.85,
            f'Cc: {corr:.2f}\np: {p_value:.4e}\n<AF3>: {mean_af3:.2f}\n<{method}>: {mean_method:.2f}',
            fontsize=10)

        plot.figure.suptitle(f'AF3 vs {method} TM-Score', fontsize=14)
        plot.figure.subplots_adjust(top=0.95)
        plot.ax_joint.set_xlabel('TM-score (AF3)', fontsize=12)
        plot.ax_joint.set_ylabel(f'TM-score ({method})', fontsize=12)

        if method == 'Boltz':
            fig_name = "Fig-3D.png"
        else:
            fig_name = f"Fig-S8-AF3-{method}-tm.png"

        save_path = os.path.join(output_dir, fig_name)
        plot.figure.savefig(save_path)
        plt.close(plot.figure)

    # ===========================
    # DockQ SCATTER PLOTS
    # ===========================

    dfs2 = []
    key_sets = []

    for csv_file in csv_files:
        label_key = csv_file.split('.')[0]
        if label_key not in custom_labels:
            continue
        label = custom_labels[label_key]
        df = pd.read_csv(os.path.join(input_dir, csv_file))
        df['Model'] = label
        df['Type'] = df['Type'].str.strip().str.lower()

        df['Type'] = df['Type'].replace({
            'rna, protein': 'rna-protein',
            'protein, rna': 'rna-protein',
            'rna, rna': 'rna, rna',
            'protein, protein': 'protein, protein'
        })

        df = df[df['Type'] != 'protein, protein']
        df = df[['PDB_ID', 'Type', 'DockQ_Score', 'Model']]
        dfs2.append(df)

        key_sets.append(set(zip(df['PDB_ID'], df['Type'])))

    common_keys = set.intersection(*key_sets)
    if not common_keys:
        raise ValueError("No common (PDB_ID, Type) entries found across models.")

    for i in range(len(dfs2)):
        dfs2[i] = dfs2[i][dfs2[i][['PDB_ID', 'Type']].apply(tuple, axis=1).isin(common_keys)]

    combined_data2 = pd.concat(dfs2, ignore_index=True)

    pivot_df = combined_data2.pivot_table(
        index=['PDB_ID', 'Type'],
        columns='Model',
        values='DockQ_Score'
    ).reset_index()

    type_styles = {
        'rna, rna': {'color': sns.color_palette("tab10")[4], 'marker': '^', 'label': 'RNA–RNA'},
        'rna-protein': {'color': sns.color_palette("tab10")[0], 'marker': 'o', 'label': 'RNA–protein'}
    }

    for model in ['HF3', 'RF2NA', 'Boltz']:
        if model not in pivot_df.columns:
            print(f"Skipping {model} for DockQ, not found.")
            continue

        df_pair = pivot_df.dropna(subset=['AF3', model]).copy()
        df_pair = df_pair.sort_values(by='AF3').reset_index(drop=True)
        df_pair['moving_avg_AF3'] = df_pair['AF3'].rolling(window=10).mean()
        df_pair[f'moving_avg_{model}'] = df_pair[model].rolling(window=10).mean()

        plot = sns.jointplot(
            data=df_pair,
            x='AF3',
            y=model,
            kind='scatter',
            height=7,
            space=0
        )

        for artist in plot.ax_joint.collections:
            artist.remove()

        for t, style in type_styles.items():
            subset = df_pair[df_pair['Type'] == t]
            if not subset.empty:
                plot.ax_joint.scatter(
                    subset['AF3'],
                    subset[model],
                    color=style['color'],
                    marker=style['marker'],
                    label=style['label'],
                    alpha=0.9,
                    s=60,
                    edgecolor='none'
                )

        maxval = max(df_pair['AF3'].max(), df_pair[model].max())
        maxval = 1 if maxval < 1 else 100 if maxval > 30 else maxval
        sns.lineplot(x=[0, maxval], y=[0, maxval], ax=plot.ax_joint, linestyle='--', color='red')
        sns.lineplot(
            data=df_pair,
            x='moving_avg_AF3',
            y=f'moving_avg_{model}',
            ax=plot.ax_joint,
            color='purple'
        )

        cutoff = maxval / 10
        for _, row in df_pair.iterrows():
            if abs(row['AF3'] - row[model]) > cutoff:
                plot.ax_joint.text(row['AF3'], row[model], row['PDB_ID'], fontsize=8, alpha=0.7)

        t_stat, p_value = ttest_rel(df_pair['AF3'], df_pair[model])
        corr = np.corrcoef(df_pair['AF3'], df_pair[model])[0, 1]
        mean_af3 = combined_data2[combined_data2['Model'] == 'AF3']['DockQ_Score'].mean()
        mean_model = combined_data2[combined_data2['Model'] == model]['DockQ_Score'].mean()

        plot.ax_joint.text(maxval * 0.1, maxval * 0.9,
            f'Cc: {corr:.2f}\np: {p_value:.4e}\n<AF3>: {mean_af3:.2f}\n<{model}>: {mean_model:.2f}',
            fontsize=9)

        plot.figure.suptitle(f'AF3 vs {model} DockQ', fontsize=14)
        plot.figure.subplots_adjust(top=0.95)
        plot.ax_joint.set_xlabel('DockQ Score (AF3)', fontsize=12)
        plot.ax_joint.set_ylabel(f'DockQ Score ({model})', fontsize=12)
        plot.ax_joint.legend(loc='lower right', fontsize=9)

        if model == 'Boltz':
            fig_name = "Fig-3E.png"
        else:
            fig_name = f"Fig-S8-AF3-{model}-dockq.png"

        save_path = os.path.join(output_dir, fig_name)
        plot.figure.savefig(save_path)
        plt.close(plot.figure)


########################### RNA COMPLEXES - FIG 3C ###########################


def interface_violinplots(input_dir: str, output_dir: str) -> None:
    """
    Creates and saves a violin plot of DockQ scores grouped by model and interaction type,
    and outputs mean scores as a CSV file.

    Parameters:
    - input_dir: Path to the directory containing the input CSV files.
    - output_dir: Path to the directory where the output figure and table will be saved.
    """
    # File names and model labels
    file_labels = {
        'usf1_af3.csv': 'AF3',
        'usf1_hf3.csv': 'HF3',
        'usf1_rf2na.csv': 'RF2NA',
        'usf1_boltz.csv': 'Boltz'
    }

    df_list = []
    pdb_type_sets = []

    # Load and process each CSV
    for file_name, model_label in file_labels.items():
        file_path = os.path.join(input_dir, file_name)
        df = pd.read_csv(file_path)
        df['Model'] = model_label

        # Normalize Type
        df['Type'] = df['Type'].str.lower().replace({
            'rna, protein': 'RNA/Protein',
            'protein, rna': 'RNA/Protein',
            'rna, rna': 'RNA/RNA'
        })

        # Keep only RNA/Protein and RNA/RNA
        df = df[df['Type'].isin(['RNA/Protein', 'RNA/RNA'])]

        # Track (PDB_ID, Type)
        pdb_type_sets.append(set(zip(df['PDB_ID'], df['Type'])))

        df_list.append(df[['PDB_ID', 'DockQ_Score', 'Type', 'Model']])

    # Find common (PDB_ID, Type) pairs across all models
    common_pairs = set.intersection(*pdb_type_sets)

    # Combine all data and filter by common pairs
    combined_df = pd.concat(df_list, ignore_index=True)
    combined_df['Pair'] = list(zip(combined_df['PDB_ID'], combined_df['Type']))
    combined_df = combined_df[combined_df['Pair'].isin(common_pairs)].drop(columns='Pair')

    # Create combined Model-Type column
    combined_df['Model_Type'] = combined_df['Model'] + '-' + combined_df['Type']

    # Order models by mean DockQ for RNA/Protein
    ordered_models = (
        combined_df[combined_df['Type'] == 'RNA/Protein']
        .groupby('Model')['DockQ_Score']
        .mean()
        .sort_values(ascending=False)
        .index.tolist()
    )

    # Build x-axis order
    model_type_order = []
    for model in ordered_models:
        model_type_order.append(model + '-RNA/Protein')
        model_type_order.append(model + '-RNA/RNA')

    # Plot
    plt.figure(figsize=(14, 6))

    # Define color palette
    palette = {
        model + '-RNA/Protein': sns.color_palette()[0] for model in file_labels.values()
    }
    palette.update({
        model + '-RNA/RNA': sns.color_palette()[4] for model in file_labels.values()
    })

    sns.violinplot(
        x='Model_Type',
        y='DockQ_Score',
        data=combined_df,
        order=model_type_order,
        cut=0,
        palette=palette,
        width=0.6
    )

    # Add vertical separators
    for i in range(1, len(model_type_order), 2):
        plt.axvline(x=i - 0.5 + 1, color='gray', linestyle='--', linewidth=0.8)

    # Formatting
    plt.xticks(
        range(len(model_type_order)),
        model_type_order,
        fontsize=13, rotation=0, ha='center'
    )
    plt.ylabel('DockQ Score', fontsize=14, fontweight='bold')
    plt.xlabel('')
    plt.title('Individual Chain Interface DockQ Scores by Interaction Type', fontsize=15)
    plt.tight_layout()

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Save figure
    output_figure_path = os.path.join(output_dir, 'Fig-3C.png')
    plt.savefig(output_figure_path, dpi=300)
    plt.close()

    # Compute mean DockQ scores
    mean_by_type = (
        combined_df
        .groupby(['Model', 'Type'])['DockQ_Score']
        .mean()
        .reset_index()
        .pivot(index='Model', columns='Type', values='DockQ_Score')
    )

    mean_combined = (
        combined_df
        .groupby('Model')['DockQ_Score']
        .mean()
        .rename('Combined')
    )

    mean_scores = mean_by_type.merge(mean_combined, left_index=True, right_index=True)

    # Save mean scores table
    output_table_path = os.path.join(output_dir, 'Table-S6.csv')
    mean_scores.round(3).to_csv(output_table_path)




########################### SINGLE RNA - TABLE S1 ###########################

def tm_score_table(input_dir: str, output_dir: str):
    os.makedirs(output_dir, exist_ok=True)

    # Label mapping
    custom_labels = {
        'tm_af3': 'AF3',
        'tm_chai': 'Chai',
        'tm_hf3': 'HF3',
        'tm_rf2na': 'RF2NA',
        'tm_boltz': 'Boltz',
        'tm_rhofold': 'RhoFold+',
        'tm_nufold': 'NuFold',
        'tm_trrosetta': 'trRosettaRNA'
    }

    # Collect all dataframes and track common PDB_IDs
    dfs = []
    csv_files = [f for f in os.listdir(input_dir) if f.endswith('.csv')]
    common_ids = None

    for file in csv_files:
        file_path = os.path.join(input_dir, file)
        df = pd.read_csv(file_path)

        base_name = file.replace('.csv', '')
        model_name = custom_labels.get(base_name, base_name)

        df = df[['PDB_ID', 'TM_Score']].drop_duplicates()
        df.rename(columns={'TM_Score': model_name}, inplace=True)

        dfs.append(df)

        if common_ids is None:
            common_ids = set(df['PDB_ID'])
        else:
            common_ids.intersection_update(set(df['PDB_ID']))

    # Merge DataFrames keeping only common PDB_IDs
    joint_df = None
    for df in dfs:
        df = df[df['PDB_ID'].isin(common_ids)]
        if joint_df is None:
            joint_df = df
        else:
            joint_df = pd.merge(joint_df, df, on='PDB_ID', how='inner')

    # Rename and sort
    joint_df = joint_df.rename(columns={'PDB_ID': 'PDB'}).sort_values(by='PDB')

    # Convert score columns to numeric
    for col in joint_df.columns[1:]:
        joint_df[col] = pd.to_numeric(joint_df[col], errors='coerce')

    # Save CSV
    csv_path = os.path.join(output_dir, 'Table-S1.csv')
    joint_df.to_csv(csv_path, index=False)




########################### SINGLE RNA - TABLE S1 ###########################

def tm_and_dockq_tables(input_dir: str, output_dir: str):
    os.makedirs(output_dir, exist_ok=True)

    label_map = {
        'usf1_af3': 'AF3',
        'usf1_hf3': 'HF3',
        'usf1_rf2na': 'RFNA',
        'usf1_boltz': 'Boltz'
    }

    tm_dfs = []
    dockq_dfs = []
    csv_files = [f for f in os.listdir(input_dir) if f.endswith('.csv')]
    common_ids = None

    for file in csv_files:
        file_path = os.path.join(input_dir, file)
        base_name = file.replace('.csv', '')
        model_name = label_map.get(base_name, base_name)

        df = pd.read_csv(file_path)

        # Process DockQ score
        if 'Total_DockQ_Score' in df.columns:
            dockq_df = df[['PDB_ID', 'Total_DockQ_Score']].drop_duplicates()
            dockq_df.rename(columns={'Total_DockQ_Score': model_name}, inplace=True)
            dockq_dfs.append(dockq_df)

        # Process TM score
        if 'TM_Score' in df.columns:
            tm_df = df[['PDB_ID', 'TM_Score']].drop_duplicates()
            tm_df.rename(columns={'TM_Score': model_name}, inplace=True)
            tm_dfs.append(tm_df)

        # Track common PDB_IDs
        if common_ids is None:
            common_ids = set(df['PDB_ID'])
        else:
            common_ids.intersection_update(set(df['PDB_ID']))

    def merge_and_save(dfs, filename):
        joint_df = None
        for df in dfs:
            df = df[df['PDB_ID'].isin(common_ids)]
            if joint_df is None:
                joint_df = df
            else:
                joint_df = pd.merge(joint_df, df, on='PDB_ID', how='inner')

        joint_df = joint_df.rename(columns={'PDB_ID': 'PDB'})
        for col in joint_df.columns[1:]:
            joint_df[col] = pd.to_numeric(joint_df[col], errors='coerce')

        joint_df = joint_df.dropna().sort_values(by='PDB')
        out_path = os.path.join(output_dir, filename)
        joint_df.to_csv(out_path, index=False)

    merge_and_save(tm_dfs, 'Table-S3.csv')
    merge_and_save(dockq_dfs, 'Table-S4.csv')




########################### RANKS - Figure S1 ###########################

def ranks_af3_single(input_dir: str, output_dir: str):
    os.makedirs(output_dir, exist_ok=True)

    tm_path = os.path.join(input_dir, "af3_single.csv")
    conf_path = os.path.join(input_dir, "af3_confidences_single.csv")

    df_tm = pd.read_csv(tm_path)
    df_conf = pd.read_csv(conf_path)

    # Convert scores to numeric
    df_tm['TM_Score'] = pd.to_numeric(df_tm['TM_Score'], errors='coerce')
    df_conf['ranking_score'] = pd.to_numeric(df_conf['ranking_score'], errors='coerce')

    # Merge on common columns
    df_merged = df_tm.merge(df_conf, on=["Predicted_PDB", "Reference_ID"], how="inner")

    # Step 1: Max ranking score per Reference_ID
    df_rank_max = df_merged.groupby("Reference_ID")["ranking_score"].max().reset_index()

    # Step 2: Keep entries with this max ranking
    df_highest_ranked = df_merged.merge(df_rank_max, on=["Reference_ID", "ranking_score"], how="inner")

    # Step 3: Select highest TM score among them
    df_max_tm = df_highest_ranked.loc[
        df_highest_ranked.groupby("Reference_ID")["TM_Score"].idxmax(),
        ["Reference_ID", "TM_Score"]
    ]

    # Order by median TM_Score for boxplot
    median_order = df_tm.groupby('Reference_ID')['TM_Score'].median().sort_values(ascending=False).index

    # Prepare the figure
    plt.figure(figsize=(12, 6))
    sns.boxplot(
        x='Reference_ID', 
        y='TM_Score', 
        data=df_tm, 
        order=median_order, 
        color="steelblue"
    )

    # Align max TM scores with median order
    df_max_tm = df_max_tm.set_index("Reference_ID").loc[median_order].reset_index()

    # Overlay line plot
    plt.plot(df_max_tm["Reference_ID"], df_max_tm["TM_Score"], marker=".", linestyle="-", color="red", label="Top-ranked model")

    # Formatting
    plt.xticks(rotation=90)
    plt.xlabel('PDB ID')
    plt.ylabel('TM-Score')
    plt.title('AF3: TM-Score Distribution for single RNA (n=25)')
    plt.legend()

    # Save the figure
    output_path = os.path.join(output_dir, "Fig-S1-A.png")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()





########################### RANKS - Figure S2 ###########################


def ranks_rf2na_single(input_dir: str, output_dir: str):
    os.makedirs(output_dir, exist_ok=True)

    # File paths
    tm_path = os.path.join(input_dir, "rf2na_single.csv")
    conf_path = os.path.join(input_dir, "rf2na_confidences.csv")

    # Load CSVs
    df_tm = pd.read_csv(tm_path)
    df_conf = pd.read_csv(conf_path)

    # Ensure numeric types
    df_tm['TM_Score'] = pd.to_numeric(df_tm['TM_Score'], errors='coerce')
    df_conf['ranking_score'] = pd.to_numeric(df_conf['ranking_score'], errors='coerce')

    # Merge on Model_ID == ID
    df_merged = df_tm.merge(df_conf, left_on='Model_ID', right_on='ID', how='inner')

    # Get the best-ranked model per PDB_ID
    idx_best = df_merged.groupby('PDB_ID')['ranking_score'].idxmax()
    df_best_ranked = df_merged.loc[idx_best, ['PDB_ID', 'TM_Score']]

    # Order PDB_IDs by median TM score
    median_order = df_tm.groupby('PDB_ID')['TM_Score'].median().sort_values(ascending=False).index

    # Create boxplot
    plt.figure(figsize=(12, 6))
    sns.boxplot(x='PDB_ID', y='TM_Score', data=df_tm, order=median_order, color="steelblue")

    # Align red line to ordered PDBs
    df_best_ranked = df_best_ranked.set_index("PDB_ID").loc[median_order].reset_index()

    # Overlay top-ranked line
    plt.plot(df_best_ranked["PDB_ID"], df_best_ranked["TM_Score"], marker=".", linestyle="-", color="red", label="Top-ranked model")

    # Formatting
    plt.xticks(rotation=90)
    plt.xlabel('PDB ID')
    plt.ylabel('TM-Score')
    plt.title('RF2NA: TM-Score Distribution for single RNA')
    plt.legend()
    plt.tight_layout()

    # Save figure
    output_path = os.path.join(output_dir, "Fig-S2-B.png")
    plt.savefig(output_path, dpi=300)
    plt.close()


########################### RANKS - Figure S1-B ###########################


def ranks_boltz_single(input_dir: str, output_dir: str):
    os.makedirs(output_dir, exist_ok=True)

    # File path
    tm_path = os.path.join(input_dir, "boltz_single.csv")

    # Load and process data
    df_tm = pd.read_csv(tm_path)
    df_tm['TM_Score'] = pd.to_numeric(df_tm['TM_Score'], errors='coerce')

    # Identify top-ranked models (those ending in _0)
    df_top_ranked = df_tm[df_tm['Model_ID'].str.endswith('_0')][['PDB_ID', 'TM_Score']]

    # Order PDB_IDs by median TM_Score
    median_order = df_tm.groupby('PDB_ID')['TM_Score'].median().sort_values(ascending=False).index

    # Filter to match order
    df_top_ranked = df_top_ranked[df_top_ranked['PDB_ID'].isin(median_order)]
    df_top_ranked = df_top_ranked.set_index("PDB_ID").loc[median_order.intersection(df_top_ranked['PDB_ID'])].reset_index()

    # Create plot
    plt.figure(figsize=(12, 6))
    sns.boxplot(x='PDB_ID', y='TM_Score', data=df_tm, order=median_order, color="steelblue")
    plt.plot(df_top_ranked["PDB_ID"], df_top_ranked["TM_Score"], marker=".", linestyle="-", color="red", label="Top-ranked model")

    # Formatting
    plt.xticks(rotation=90)
    plt.xlabel('PDB ID')
    plt.ylabel('TM-Score')
    plt.title('Boltz: TM-Score Distribution for single RNA (n=5)')
    plt.legend()
    plt.tight_layout()

    # Save figure
    output_path = os.path.join(output_dir, "Fig-S1-B.png")
    plt.savefig(output_path, dpi=300)
    plt.close()



########################### RANKS - Figure S2-A ###########################


def ranks_hf3_single(input_dir: str, output_dir: str):
    os.makedirs(output_dir, exist_ok=True)

    # Load file
    tm_path = os.path.join(input_dir, "hf3_single.csv")
    df_tm = pd.read_csv(tm_path)
    df_tm['TM_Score'] = pd.to_numeric(df_tm['TM_Score'], errors='coerce')

    # Select top-ranked models (Model_ID ending in _1)
    df_top_ranked = df_tm[df_tm['Model_ID'].str.endswith('_1')][['PDB_ID', 'TM_Score']]

    # Sort PDBs by median TM_Score
    median_order = df_tm.groupby('PDB_ID')['TM_Score'].median().sort_values(ascending=False).index

    # Filter and align top-ranked data
    df_top_ranked = df_top_ranked[df_top_ranked['PDB_ID'].isin(median_order)]
    df_top_ranked = df_top_ranked.set_index("PDB_ID").loc[median_order.intersection(df_top_ranked['PDB_ID'])].reset_index()

    # Plot boxplot
    plt.figure(figsize=(12, 6))
    sns.boxplot(x='PDB_ID', y='TM_Score', data=df_tm, order=median_order, color="steelblue")
    plt.plot(df_top_ranked["PDB_ID"], df_top_ranked["TM_Score"], marker=".", linestyle="-", color="red", label="Top-ranked model")

    # Formatting
    plt.xticks(rotation=90)
    plt.xlabel('PDB ID')
    plt.ylabel('TM-Score')
    plt.title('HF3: TM-Score Distribution for single RNA (n=5)')
    plt.legend()
    plt.tight_layout()

    # Save figure
    output_path = os.path.join(output_dir, "Fig-S2-A.png")
    plt.savefig(output_path, dpi=300)
    plt.close()



########################### RANKS - Figure S3-A ###########################

def ranks_af3_complexes(input_dir: str, output_dir: str):
    os.makedirs(output_dir, exist_ok=True)

    # File paths
    tm_path = os.path.join(input_dir, "af3_complexes.csv")
    conf_path = os.path.join(input_dir, "af3_confidences_complexes.csv")

    # Load data
    df_tm = pd.read_csv(tm_path)
    df_conf = pd.read_csv(conf_path)

    # Clean and convert data
    df_tm['DockQ_Score'] = pd.to_numeric(df_tm['DockQ_Score'], errors='coerce')
    df_tm['Reference_ID'] = df_tm['Reference_ID'].str.replace(r'\.pdb$', '', regex=True)
    df_tm = df_tm.dropna(subset=['DockQ_Score'])
    df_conf['ranking_score'] = pd.to_numeric(df_conf['ranking_score'], errors='coerce')

    # Merge on predicted/reference ID
    df_merged = df_tm.merge(df_conf, on=["Predicted_PDB", "Reference_ID"], how="inner")

    # Highest-ranked model per Reference_ID
    df_max_rank = df_merged.loc[df_merged.groupby("Reference_ID")["ranking_score"].idxmax(), ["Reference_ID", "DockQ_Score"]]

    # Order for plotting
    median_order = df_tm.groupby('Reference_ID')['DockQ_Score'].median().sort_values(ascending=False).index

    # Oligo RNAs for color highlighting
    oligo_rnas = {
        "7TZR", "8UO6", "9CPI", "8HBA", "8TQX", "8K0P", "8X1V", "8SA3", "8TSV",
        "7XKM", "8EYU", "8SXL", "8HZJ", "8HZL", "8WU5", "8TDZ", "R1283v2", "R1281o"
    }

    # Assign colors
    df_tm['color'] = df_tm['Reference_ID'].apply(lambda x: "mediumorchid" if x in oligo_rnas else "steelblue")
    palette = df_tm.set_index('Reference_ID')['color'].to_dict()

    # Create boxplot
    plt.figure(figsize=(12, 6))
    sns.boxplot(
        x='Reference_ID',
        y='DockQ_Score',
        data=df_tm,
        order=median_order,
        palette=palette
    )

    # Align and plot top-ranked model DockQ scores
    df_max_rank = df_max_rank[df_max_rank["Reference_ID"].isin(median_order)]
    df_max_rank = df_max_rank.set_index("Reference_ID").reindex(median_order.intersection(df_max_rank["Reference_ID"])).reset_index()
    plt.plot(df_max_rank["Reference_ID"], df_max_rank["DockQ_Score"], marker=".", linestyle="-", color="red", label="Top-ranked model")

    # Formatting
    plt.xticks(rotation=90)
    plt.xlabel('Reference ID')
    plt.ylabel('DockQ')
    plt.title('AF3: DockQ Score Distribution for complexes (n=25)')
    plt.legend()
    plt.tight_layout()

    # Save figure
    output_path = os.path.join(output_dir, "Fig-S3-A.png")
    plt.savefig(output_path, dpi=300)
    plt.close()



########################### RANKS - Figure S3-B ###########################

def ranks_boltz_complexes(input_dir: str, output_dir: str):
    os.makedirs(output_dir, exist_ok=True)

    # File path
    tm_path = os.path.join(input_dir, "boltz_complexes.csv")

    # Load data
    df_tm = pd.read_csv(tm_path)
    df_tm['DockQ_Score'] = pd.to_numeric(df_tm['DockQ_Score'], errors='coerce')

    # Extract ranking from Model_ID (e.g., M1239v1_1.pdb → 1)
    df_tm['ranking'] = df_tm['Model_ID'].apply(
        lambda x: int(re.search(r'_(\d+)\.pdb', x).group(1)) if re.search(r'_(\d+)\.pdb', x) else None
    )

    # Get n=1 models for red line
    df_max_rank = df_tm[df_tm['ranking'] == 1][['PDB_ID', 'DockQ_Score']]

    # Order by median DockQ_Score
    median_order = df_tm.groupby('PDB_ID')['DockQ_Score'].median().sort_values(ascending=False).index

    # Filter and align top-ranked models
    df_max_rank = df_max_rank[df_max_rank['PDB_ID'].isin(median_order)]
    df_max_rank = df_max_rank.set_index("PDB_ID").loc[median_order.intersection(df_max_rank['PDB_ID'])].reset_index()

    # Oligo RNA coloring
    oligo_rnas = {
        "7TZR", "8UO6", "9CPI", "8HBA", "8TQX", "8K0P", "8X1V", "8SA3", "8TSV",
        "7XKM", "8EYU", "8SXL", "8HZJ", "8HZL", "8WU5", "8TDZ", "R1283v2", "R1281o"
    }
    df_tm['color'] = df_tm['PDB_ID'].apply(lambda x: "mediumorchid" if x in oligo_rnas else "steelblue")
    palette = df_tm.set_index('PDB_ID')['color'].to_dict()

    # Create plot
    plt.figure(figsize=(12, 6))
    sns.boxplot(
        x='PDB_ID',
        y='DockQ_Score',
        data=df_tm,
        order=median_order,
        palette=palette
    )
    plt.plot(df_max_rank["PDB_ID"], df_max_rank["DockQ_Score"], marker=".", linestyle="-", color="red", label="Top-ranked model")

    # Formatting
    plt.xticks(rotation=90)
    plt.xlabel('Reference ID')
    plt.ylabel('DockQ')
    plt.title('Boltz: DockQ Score Distribution for complexes (n=5)')
    plt.legend()
    plt.tight_layout()

    # Save figure
    output_path = os.path.join(output_dir, "Fig-S3-B.png")
    plt.savefig(output_path, dpi=300)
    plt.close()


########################### RANKS - Figure S4-A ###########################

def ranks_hf3_complexes(input_dir: str, output_dir: str):
    os.makedirs(output_dir, exist_ok=True)

    # File paths
    tm_path = os.path.join(input_dir, "hf3_complexes.csv")
    top_path = os.path.join(input_dir, "hf3_confidences_complexes.csv")

    # Load ranked models
    df_tm = pd.read_csv(tm_path)
    df_tm['DockQ_Score'] = pd.to_numeric(df_tm['DockQ_Score'], errors='coerce')
    df_tm['ranking'] = df_tm['Model_ID'].apply(
        lambda x: int(re.search(r'_(\d+)\.pdb', x).group(1)) if re.search(r'_(\d+)\.pdb', x) else None
    )

    # Load top-ranked scores
    df_toprank = pd.read_csv(top_path)
    df_toprank = df_toprank[['PDB_ID', 'Total_DockQ_Score']].rename(columns={'Total_DockQ_Score': 'DockQ_Score'})
    df_toprank = df_toprank.drop_duplicates(subset='PDB_ID', keep='first')

    # Order by median DockQ_Score
    median_order = df_tm.groupby('PDB_ID')['DockQ_Score'].median().sort_values(ascending=False).index

    # Color oligo RNAs
    oligo_rnas = {
        "7TZR", "8UO6", "9CPI", "8HBA", "8TQX", "8K0P", "8X1V", "8SA3", "8TSV",
        "7XKM", "8EYU", "8SXL", "8HZJ", "8HZL", "8WU5", "8TDZ", "R1283v2", "R1281o"
    }
    df_tm['color'] = df_tm['PDB_ID'].apply(lambda x: "mediumorchid" if x in oligo_rnas else "steelblue")
    palette = df_tm.set_index('PDB_ID')['color'].to_dict()

    # Create boxplot
    plt.figure(figsize=(12, 6))
    sns.boxplot(
        x='PDB_ID',
        y='DockQ_Score',
        data=df_tm,
        order=median_order,
        palette=palette
    )

    # Align and plot top-ranked DockQ scores
    df_toprank = df_toprank.set_index("PDB_ID").reindex(median_order, copy=False).dropna(subset=["DockQ_Score"]).reset_index()
    plt.plot(df_toprank["PDB_ID"], df_toprank["DockQ_Score"], marker=".", linestyle="-", color="red", label="Top-ranked model")

    # Formatting
    plt.xticks(rotation=90)
    plt.xlabel('Reference ID')
    plt.ylabel('DockQ')
    plt.title('HF3: DockQ Score Distribution for complexes (n=5)')
    plt.legend()
    plt.tight_layout()

    # Save figure
    output_path = os.path.join(output_dir, "Fig-S4-A.png")
    plt.savefig(output_path, dpi=300)
    plt.close()




########################### RANKS - Figure S4-B ###########################

def ranks_rf2na_complexes(input_dir: str, output_dir: str):
    os.makedirs(output_dir, exist_ok=True)

    # Input file paths
    tm_path = os.path.join(input_dir, "rf2na_complexes.csv")
    conf_path = os.path.join(input_dir, "rf2na_confidences.csv")

    # Load model scores
    df_tm = pd.read_csv(tm_path)
    df_tm['DockQ_Score'] = pd.to_numeric(df_tm['DockQ_Score'], errors='coerce')

    # Extract numeric rank from Model_ID
    df_tm['ranking'] = df_tm['Model_ID'].apply(
        lambda x: int(re.search(r'_(\d+)\.pdb', x).group(1)) if re.search(r'_(\d+)\.pdb', x) else None
    )

    # Load confidence scores
    df_conf = pd.read_csv(conf_path)
    df_conf['ranking_score'] = pd.to_numeric(df_conf['ranking_score'], errors='coerce')

    # Merge scores and confidences
    df_merged = df_tm.merge(df_conf, left_on='Model_ID', right_on='ID', how='inner')

    # Sort and get top-ranked per PDB_ID (if ties, keep the first one)
    df_merged_sorted = df_merged.sort_values(by=['PDB_ID', 'ranking_score'], ascending=[True, False])
    df_toprank = df_merged_sorted.groupby('PDB_ID', as_index=False).first()[['PDB_ID', 'DockQ_Score']]
    df_toprank = df_toprank.rename(columns={'DockQ_Score': 'Top_DockQ_Score'})

    # Order by median DockQ score
    median_order = df_tm.groupby('PDB_ID')['DockQ_Score'].median().sort_values(ascending=False).index.tolist()

    # Color oligo RNAs
    oligo_rnas = {
        "7TZR", "8UO6", "9CPI", "8HBA", "8TQX", "8K0P", "8X1V", "8SA3", "8TSV",
        "7XKM", "8EYU", "8SXL", "8HZJ", "8HZL", "8WU5", "8TDZ", "R1283v2", "R1281o"
    }
    df_tm['color'] = df_tm['PDB_ID'].apply(lambda x: "mediumorchid" if x in oligo_rnas else "steelblue")
    palette = df_tm.set_index('PDB_ID')['color'].to_dict()

    # Create the boxplot
    plt.figure(figsize=(12, 6))
    sns.boxplot(
        x='PDB_ID',
        y='DockQ_Score',
        data=df_tm,
        order=median_order,
        palette=palette
    )

    # Align and plot top-ranked DockQ scores
    df_toprank = df_toprank.set_index("PDB_ID").reindex(median_order).dropna().reset_index()
    x_positions = list(range(len(df_toprank)))  # use numeric x for line plot
    plt.plot(
        x_positions,
        df_toprank["Top_DockQ_Score"],
        marker=".",
        linestyle="-",
        color="red",
        label="Top-ranked model"
    )

    # Formatting
    plt.xticks(ticks=x_positions, labels=df_toprank["PDB_ID"], rotation=90)
    plt.xlabel('Reference ID')
    plt.ylabel('DockQ')
    plt.title('RF2NA: DockQ Score Distribution for complexes (n=5)')
    plt.legend()
    plt.tight_layout()

    # Save figure
    output_path = os.path.join(output_dir, "Fig-S4-B.png")
    plt.savefig(output_path, dpi=300)
    plt.close()

########################### CONFIDENCE SCORING - Figure 4B, 4D, 4F ###########################

def iptm_vs_dockq(input_dir: str, output_dir: str):
    os.makedirs(output_dir, exist_ok=True)

    # File names and output figure mapping
    files = {
        "usf1_af3.csv": "Fig-4B.png",
        "usf1_boltz.csv": "Fig-4D.png",
        "usf1_hf3.csv": "Fig-4F.png"
    }

    plot_titles = {
        "usf1_af3.csv": "AF3",
        "usf1_boltz.csv": "Boltz-1",
        "usf1_hf3.csv": "HF3"
    }

    # Load and normalize Type labels
    type_path = os.path.join(input_dir, "usf1_boltz.csv")
    type_df = pd.read_csv(type_path)[["PDB_ID", "Type"]]
    type_df["Type"] = type_df["Type"].replace({
        "protein, rna": "rna, protein",
        "protein, protein": "rna, protein"
    })

    for filename, fig_name in files.items():
        df = pd.read_csv(os.path.join(input_dir, filename))
        df = df[["PDB_ID", "iptm", "Total_DockQ_Score"]].dropna()
        df = df.drop_duplicates(subset="PDB_ID", keep="first")
        df = pd.merge(df, type_df, on="PDB_ID")
        df = df.replace([np.inf, -np.inf], np.nan).dropna()

        # Compute global and per-Type correlations
        corr_all, _ = pearsonr(df["iptm"], df["Total_DockQ_Score"])
        l2_loss = np.mean((df["Total_DockQ_Score"] - df["iptm"]) ** 2)

        corr_rna_rna = np.nan
        corr_rna_protein = np.nan

        if not df[df["Type"] == "rna, rna"].empty:
            corr_rna_rna, _ = pearsonr(
                df[df["Type"] == "rna, rna"]["iptm"],
                df[df["Type"] == "rna, rna"]["Total_DockQ_Score"]
            )

        if not df[df["Type"] == "rna, protein"].empty:
            corr_rna_protein, _ = pearsonr(
                df[df["Type"] == "rna, protein"]["iptm"],
                df[df["Type"] == "rna, protein"]["Total_DockQ_Score"]
            )

        sns.set(style="white")
        g = sns.JointGrid(data=df, x="iptm", y="Total_DockQ_Score", height=7)

        palette = {
            "rna, rna": "#ee82ee",
            "rna, protein": "#1f77b4"
        }

        # Main scatterplot
        g.plot_joint(
            sns.scatterplot,
            data=df,
            hue="Type",
            palette=palette,
            alpha=0.8,
            s=60
        )

        # Marginal KDEs
        sns.kdeplot(data=df, x="iptm", ax=g.ax_marg_x, hue="Type", palette=palette, fill=True, common_norm=False, alpha=0.4, legend=False)
        sns.kdeplot(data=df, y="Total_DockQ_Score", ax=g.ax_marg_y, hue="Type", palette=palette, fill=True, common_norm=False, alpha=0.4, legend=False)

        # Diagonal perfect correlation line
        maxval = max(df["iptm"].max(), df["Total_DockQ_Score"].max())
        g.ax_joint.plot([0, maxval], [0, maxval], linestyle="--", color="red")

        # Correlation text block
        corr_text = f"Cc: {corr_all:.2f}\nL2: {l2_loss:.2f}"
        if not np.isnan(corr_rna_rna):
            corr_text += f"\nCc (RNA/RNA): {corr_rna_rna:.2f}"
        if not np.isnan(corr_rna_protein):
            corr_text += f"\nCc (RNA/Protein): {corr_rna_protein:.2f}"

        g.ax_joint.text(
            0.05, 0.95,
            corr_text,
            transform=g.ax_joint.transAxes,
            fontsize=12,
            verticalalignment='top'
        )

        # Labels and titles
        g.set_axis_labels("ipTM", "DockQ", fontsize=13)
        g.figure.suptitle(f"ipTM vs DockQ {plot_titles[filename]}", fontsize=14)
        g.figure.subplots_adjust(top=0.95)
        g.ax_joint.legend(title="", loc="best")

        # Save figure
        output_path = os.path.join(output_dir, fig_name)
        plt.savefig(output_path, dpi=300)
        plt.close()


########################### CONFIDENCE SCORING - Figure 4A, 4C, 4E ###########################

def ptm_vs_tm(input_dir: str, output_dir: str):
    os.makedirs(output_dir, exist_ok=True)

    files = ["tm_af3.csv", "tm_boltz.csv", "tm_hf3.csv"]
    labels = ["AF3", "Boltz", "HF3"]
    output_names = ["Fig-4A.png", "Fig-4C.png", "Fig-4E.png"]

    # Load all data
    all_data = []
    for file, label in zip(files, labels):
        df = pd.read_csv(os.path.join(input_dir, file))
        df = df[["PDB_ID", "ptm", "TM_Score"]].dropna()
        df = df.drop_duplicates(subset="PDB_ID", keep="first")
        df["method"] = label
        all_data.append(df)

    # Find common PDB_IDs across all datasets
    common_ids = set(all_data[0]["PDB_ID"])
    for df in all_data[1:]:
        common_ids &= set(df["PDB_ID"])

    # Filter to common IDs
    filtered_data = [df[df["PDB_ID"].isin(common_ids)].copy() for df in all_data]

    # Generate individual plots
    sns.set(style="white")
    for df, label, out_name in zip(filtered_data, labels, output_names):
        df.rename(columns={"ptm": "pTM", "TM_Score": "TM"}, inplace=True)
        df = df.replace([np.inf, -np.inf], np.nan).dropna()

        g = sns.JointGrid(data=df, x="pTM", y="TM", height=7, space=0)

        # Scatter plot
        g.plot_joint(sns.scatterplot, s=60, alpha=0.8, color="#d382ee")

        # Marginal KDEs
        sns.kdeplot(data=df, x="pTM", ax=g.ax_marg_x, fill=True, alpha=0.4, color="#d382ee")
        sns.kdeplot(data=df, y="TM", ax=g.ax_marg_y, fill=True, alpha=0.4, color="#d382ee")

        # Diagonal reference
        maxval = max(df["pTM"].max(), df["TM"].max())
        g.ax_joint.plot([0, maxval], [0, maxval], linestyle="--", color="red")

        # Stats
        corr, _ = pearsonr(df["pTM"], df["TM"])
        l2_loss = np.mean((df["TM"] - df["pTM"]) ** 2)
        g.ax_joint.text(
            0.05, 0.95,
            f"Cc: {corr:.2f}\nL2: {l2_loss:.2f}",
            transform=g.ax_joint.transAxes,
            fontsize=12,
            verticalalignment='top'
        )

        g.ax_joint.set_xlim(-0.1, 1.15)
        g.ax_joint.set_ylim(-0.1, 1.15)

        g.set_axis_labels("pTM", "TM-score", fontsize=13)
        g.figure.suptitle(f"pTM vs TM-score {label}", fontsize=14)
        g.figure.subplots_adjust(top=0.95)

        # Save figure
        output_path = os.path.join(output_dir, out_name)
        plt.savefig(output_path, dpi=300)
        plt.close()

########################### CONFIDENCE SCORING - Figure 5A, S11-13 ###########################

def training_dependence_violin_single(similarity_dir: str, input_dir: str, output_dir: str):
    os.makedirs(output_dir, exist_ok=True)

    # Input similarity file (Training TM-Scores)
    file1 = os.path.join(similarity_dir, "similarity_tm_scores_single.csv")
    df1 = pd.read_csv(file1)
    df1_max = df1.loc[df1.groupby("ID1")["TM-Score"].idxmax()].rename(columns={"TM-Score": "Training_TM_Score"})

    # Map of input CSVs to display names and output figure names
    files = {
        "tm_af3.csv": ("AF3", "Fig-5A.png"),
        "tm_boltz.csv": ("Boltz", "Fig-S11-A.png"),
        "tm_rf2na.csv": ("RF2NA", "Fig-S11-B.png"),
        "tm_chai.csv": ("Chai", "Fig-S12-A.png"),
        "tm_rhofold.csv": ("RhoFold+", "Fig-S12-B.png"),
        "tm_nufold.csv": ("NuFold", "Fig-S13-A.png"),
        "tm_hf3.csv": ("HF3", "Fig-S13-B.png"),
        "tm_trrosetta.csv": ("trRosettaRNA", "Fig-S13-C.png"),
    }

    for filename, (label, output_name) in files.items():
        file2 = os.path.join(input_dir, filename)
        df2 = pd.read_csv(file2)
        df2_filtered = df2[df2["PDB_ID"].isin(df1_max["ID1"])].drop_duplicates(subset="PDB_ID")

        # Merge training + model TM-scores
        merged_df = df1_max.merge(df2_filtered, left_on="ID1", right_on="PDB_ID", how="left")
        merged_df = merged_df.rename(columns={"TM_Score": "Model_TM_Score"})

        # Drop rows with missing values
        corr_df = merged_df.dropna(subset=["Training_TM_Score", "Model_TM_Score"])
        cc, pval = pearsonr(corr_df["Training_TM_Score"], corr_df["Model_TM_Score"])

        # Bin by similarity to training set
        df_25 = merged_df[merged_df["Training_TM_Score"] < 0.25].copy()
        df_50 = merged_df[(merged_df["Training_TM_Score"] >= 0.25) & (merged_df["Training_TM_Score"] < 0.5)].copy()
        df_75 = merged_df[(merged_df["Training_TM_Score"] >= 0.5) & (merged_df["Training_TM_Score"] < 0.75)].copy()
        df_100 = merged_df[merged_df["Training_TM_Score"] >= 0.75].copy()

        # Assign categories
        df_25["Category"] = "x < 0.25"
        df_50["Category"] = "0.25 < x < 0.50"
        df_75["Category"] = "0.50 < x < 0.75"
        df_100["Category"] = "x > 0.75"

        df_combined = pd.concat([df_25, df_50, df_75, df_100])

        # Plot violin plot
        plt.figure(figsize=(7, 6))
        sns.violinplot(x="Category", y="Model_TM_Score", data=df_combined, inner="quart", cut=0)

        plt.xlabel("TM-Score Category", fontsize=11)
        plt.ylabel("TM-Score", fontsize=11)
        plt.title(f"Training Set Dependence - {label} (Single RNA)", fontsize=12)

        # Annotate Pearson correlation
        plt.text(
            0.05, 0.95,
            f"Cc: {cc:.2f}\n p: {pval:.2e}",
            transform=plt.gca().transAxes,
            fontsize=9,
            verticalalignment="top"
        )

        plt.tight_layout()
        output_path = os.path.join(output_dir, output_name)
        plt.savefig(output_path, dpi=300)
        plt.close()



########################### CONFIDENCE SCORING SCATTER - Figure 5B, S11-13-scatter ###########################

def training_dependence_scatter_single(similarity_dir: str, input_dir: str, output_dir: str):
    os.makedirs(output_dir, exist_ok=True)

    # Load training TM-scores
    tm_path = os.path.join(similarity_dir, "similarity_tm_scores_single.csv")
    tm_df = pd.read_csv(tm_path)
    highest_tm_scores = tm_df.loc[tm_df.groupby("ID1")["TM-Score"].idxmax()].rename(columns={"TM-Score": "Training_TM_Score"})

    # Input files and corresponding method/output names
    file_map = {
        "tm_af3.csv": ("AF3", "Fig-5B.png"),
        "tm_boltz.csv": ("Boltz", "Fig-S11-A-scatter.png"),
        "tm_rf2na.csv": ("RF2NA", "Fig-S11-B-scatter.png"),
        "tm_chai.csv": ("Chai", "Fig-S12-A-scatter.png"),
        "tm_rhofold.csv": ("RhoFold+", "Fig-S12-B-scatter.png"),
        "tm_nufold.csv": ("NuFold", "Fig-S13-A-scatter.png"),
        "tm_hf3.csv": ("HF3", "Fig-S13-B-scatter.png"),
        "tm_trrosetta.csv": ("trRosettaRNA", "Fig-S13-C-scatter.png"),
    }

    for filename, (method_name, output_name) in file_map.items():
        file_path = os.path.join(input_dir, filename)
        usalign_df = pd.read_csv(file_path)

        merged_df = highest_tm_scores.merge(usalign_df, left_on="ID1", right_on="PDB_ID", how="inner")
        merged_df = merged_df.rename(columns={"TM_Score": f"{method_name}_TM_Score"})

        # Compute correlation
        cc, p_value = pearsonr(merged_df["Training_TM_Score"], merged_df[f"{method_name}_TM_Score"])

        # Plot
        plt.figure(figsize=(7, 6))
        sns.regplot(
            x=merged_df["Training_TM_Score"],
            y=merged_df[f"{method_name}_TM_Score"],
            scatter_kws={"s": 50},
            line_kws={"color": "red"}
        )

        plt.xlabel("Training set TM-Score")
        plt.ylabel("TM-Score")
        plt.title(f"Training set TM-Score vs True TM-scores {method_name} (Single RNA)")

        plt.text(
            0.05, 0.9,
            f"Cc: {cc:.2f}\np: {p_value:.2e}",
            transform=plt.gca().transAxes,
            fontsize=9
        )

        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, output_name), dpi=300)
        plt.close()





########################### CONFIDENCE SCORING  - Figure 5C, S14 ###########################

def training_dependence_violin_complexes(similarity_dir: str, input_dir: str, output_dir: str):
    os.makedirs(output_dir, exist_ok=True)

    # Load training similarity scores
    file1 = os.path.join(similarity_dir, "similarity_tm_scores_complexes.csv")
    df1 = pd.read_csv(file1)
    df1_max = df1.loc[df1.groupby("ID1")["TM-Score"].idxmax()].rename(columns={"TM-Score": "Training_TM_Score"})

    # Map input CSVs to method label and output filename
    files = {
        "usf1_af3.csv": ("AF3", "Fig-5C.png"),
        "usf1_boltz.csv": ("Boltz", "Fig-S14-A.png"),
        "usf1_hf3.csv": ("HF3", "Fig-S14-B.png"),
        "usf1_rf2na.csv": ("RF2NA", "Fig-S14-C.png"),
    }

    for filename, (label, output_name) in files.items():
        file2 = os.path.join(input_dir, filename)
        df2 = pd.read_csv(file2)

        df2_filtered = df2[df2["PDB_ID"].isin(df1_max["ID1"])].drop_duplicates(subset="PDB_ID")
        merged_df = df1_max.merge(df2_filtered, left_on="ID1", right_on="PDB_ID", how="left")
        merged_df = merged_df.rename(columns={"TM_Score": f"{label}_TM_Score"})

        corr_df = merged_df.dropna(subset=["Training_TM_Score", f"{label}_TM_Score"])
        cc, pval = pearsonr(corr_df["Training_TM_Score"], corr_df[f"{label}_TM_Score"])

        # Categorize by training similarity
        df_25 = merged_df[merged_df["Training_TM_Score"] < 0.25].copy()
        df_50 = merged_df[(merged_df["Training_TM_Score"] >= 0.25) & (merged_df["Training_TM_Score"] < 0.5)].copy()
        df_75 = merged_df[(merged_df["Training_TM_Score"] >= 0.5) & (merged_df["Training_TM_Score"] < 0.75)].copy()
        df_100 = merged_df[merged_df["Training_TM_Score"] >= 0.75].copy()

        for d, label_ in zip([df_25, df_50, df_75, df_100],
                             ["x < 0.25", "0.25 < x < 0.50", "0.50 < x < 0.75", "x > 0.75"]):
            d["Category"] = label_

        df_combined = pd.concat([df_25, df_50, df_75, df_100])

        # Plot
        plt.figure(figsize=(7, 6))
        sns.violinplot(x="Category", y=f"{label}_TM_Score", data=df_combined, inner="quart", cut=0)
        plt.xlabel("TM-Score Category", fontsize=11)
        plt.ylabel("TM-Score", fontsize=11)
        plt.title(f"Training set dependence - {label} (RNA Complexes)", fontsize=12)

        plt.text(
            0.05, 0.95,
            f"Cc: {cc:.2f}\n p: {pval:.2e}",
            transform=plt.gca().transAxes,
            fontsize=9,
            verticalalignment="top"
        )

        plt.tight_layout()
        output_path = os.path.join(output_dir, output_name)
        plt.savefig(output_path, dpi=300)
        plt.close()

########################### CONFIDENCE SCORING SCATTER - Figure 5C, S14-scatter ###########################

def training_dependence_scatter_complexes(similarity_dir: str, input_dir: str, output_dir: str):
    os.makedirs(output_dir, exist_ok=True)

    # Load similarity TM-scores
    df1 = pd.read_csv(os.path.join(similarity_dir, "similarity_tm_scores_complexes.csv"))
    df1_max = df1.groupby("ID1")["TM-Score"].max().reset_index().rename(columns={"TM-Score": "Training_TM_Score"})

    # Map files to labels and output names
    files = {
        "usf1_af3.csv": ("AF3", "Fig-5D.png"),
        "usf1_boltz.csv": ("Boltz", "Fig-S14-A-scatter.png"),
        "usf1_hf3.csv": ("HF3", "Fig-S14-B-scatter.png"),
        "usf1_rf2na.csv": ("RF2NA", "Fig-S14-C-scatter.png"),
    }

    for filename, (label, output_name) in files.items():
        df2 = pd.read_csv(os.path.join(input_dir, filename))
        df2_filtered = df2[df2["PDB_ID"].isin(df1_max["ID1"])].drop_duplicates(subset="PDB_ID")

        # Merge and rename
        merged_df = df1_max.merge(df2_filtered, left_on="ID1", right_on="PDB_ID", how="inner")
        merged_df = merged_df.rename(columns={"TM_Score": f"{label}_TM_Score"})

        # Compute correlation
        corr_df = merged_df.dropna(subset=["Training_TM_Score", f"{label}_TM_Score"])
        cc, pval = pearsonr(corr_df["Training_TM_Score"], corr_df[f"{label}_TM_Score"])

        # Plot
        plt.figure(figsize=(8, 6))
        sns.regplot(
            x=merged_df["Training_TM_Score"],
            y=merged_df[f"{label}_TM_Score"],
            scatter_kws={"s": 50},
            line_kws={"color": "red"}
        )

        plt.xlabel("Training set TM-Score")
        plt.ylabel("TM-Score")  # From model
        plt.title(f"Training set TM-score vs True TM-scores {label} (RNA Complexes)")

        # Annotate with Cc and p-value
        plt.text(
            0.05, 0.9,
            f"Cc: {cc:.2f}\np: {pval:.2e}",
            transform=plt.gca().transAxes,
            fontsize=9
        )

        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, output_name), dpi=300)
        plt.close()
