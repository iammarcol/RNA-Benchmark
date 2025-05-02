from benchmark_analysis import plot_tm_scores, plot_success_rates, compare_to_boltz, plot_tm_chains, rna_length_scatter, plot_tm_dockq_complexes, plot_tm_dockq_scatter, interface_violinplots, tm_score_table, tm_and_dockq_tables, ranks_af3_single, ranks_rf2na_single, ranks_boltz_single, ranks_hf3_single, ranks_af3_complexes, ranks_hf3_complexes, ranks_boltz_complexes, ranks_rf2na_complexes, iptm_vs_dockq, ptm_vs_tm, training_dependence_violin_single, training_dependence_scatter_single, training_dependence_violin_complexes, training_dependence_scatter_complexes

input_dir= "../OUTPUTS/scores-single"
output_dir= "../OUTPUTS/figures"
chain_tm = "../OUTPUTS/scores-chains"
tm_complexes = "../OUTPUTS/scores-complexes"
ranks_dir = "../OUTPUTS/ranks-csvs"
similarity_dir = "../OUTPUTS/similarity-csv"

# Figure 1
plot_tm_scores(
    input_dir=input_dir,
    chain_tm_dir=chain_tm,
    output_dir=output_dir
)

# Table S5
plot_success_rates(
    input_dir=input_dir,
    chain_tm_dir=chain_tm,
    output_dir=output_dir
)

compare_to_boltz(
    input_dir=input_dir,
    output_dir=output_dir
)

plot_tm_chains(
    chain_tm_dir=chain_tm,
    output_dir=output_dir
)

rna_length_scatter(
    tm_score_dir=input_dir,
    output_dir=output_dir
)

plot_tm_dockq_complexes(
    input_dir=tm_complexes,
    output_dir=output_dir
)

plot_tm_dockq_scatter(
    input_dir=tm_complexes,
    output_dir=output_dir
)

interface_violinplots(
    input_dir=tm_complexes,
    output_dir=output_dir
) 

# Table S1
tm_score_table(
    input_dir=input_dir,
    output_dir=output_dir
)

# Table S3, S4
tm_and_dockq_tables(
    input_dir=tm_complexes,
    output_dir=output_dir
)

# Figure S1-A
ranks_af3_single(
    input_dir=ranks_dir,
    output_dir=output_dir
)

# Figure S2-B
ranks_rf2na_single(
    input_dir=ranks_dir,
    output_dir=output_dir
)

# Figure S1-B
ranks_boltz_single( 
    input_dir=ranks_dir,
    output_dir=output_dir)


# Figure S2-A
ranks_hf3_single( 
    input_dir=ranks_dir,
    output_dir=output_dir)

# Figure S3-A
ranks_af3_complexes(
    input_dir=ranks_dir,
    output_dir=output_dir)

# Figure S3-B
ranks_boltz_complexes(
    input_dir=ranks_dir,
    output_dir=output_dir)

# Figure S4-A
ranks_hf3_complexes(
    input_dir=ranks_dir,
    output_dir=output_dir)


# Figure S4-B
ranks_rf2na_complexes(
    input_dir=ranks_dir,
    output_dir=output_dir
)

# Figure 4B, 4D, 4F
iptm_vs_dockq(
    input_dir=tm_complexes,
    output_dir=output_dir
)

# Figure 4A, 4C, 4D
ptm_vs_tm(
    input_dir=input_dir,
    output_dir=output_dir
)

# Figure 5A, S12-14
training_dependence_violin_single(
    input_dir=input_dir,
    similarity_dir=similarity_dir,
    output_dir=output_dir
)

#  Figure 5B, S12-14-scatter
training_dependence_scatter_single(
    input_dir=input_dir,
    similarity_dir=similarity_dir,
    output_dir=output_dir
)


# Figure 5C, S14
training_dependence_violin_complexes(
    input_dir=tm_complexes,
    similarity_dir=similarity_dir,
    output_dir=output_dir
)

# Figure 5D, S14-scatter
training_dependence_scatter_complexes(
    input_dir=tm_complexes,
    similarity_dir=similarity_dir,
    output_dir=output_dir
)
