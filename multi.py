import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
import io
import base64
import os
import seaborn as sns
import matplotlib.colors as mcolors

# statsmodels for ANOVA and Tukey
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd

###############################################################################
# Collision/Overlap Helpers
###############################################################################
def intersects_any(box, occupied_boxes):
    """
    Returns True if 'box' intersects any bounding box in 'occupied_boxes'.
    Overlap occurs if x-ranges and y-ranges overlap.
    """
    min_x, max_x, min_y, max_y = box
    for (ox_min, ox_max, oy_min, oy_max) in occupied_boxes:
        overlap_x = not (max_x < ox_min or min_x > ox_max)
        overlap_y = not (max_y < oy_min or min_y > oy_max)
        if overlap_x and overlap_y:
            return True
    return False

def draw_significance_indicator_collision_aware(
    ax,
    comparisons,
    bar_positions,
    bar_width,
    bar_height_lookup,
    occupied_boxes,
    max_bar_per_gene,
    base_offset=0.2,
    step=0.15,
    cap_size=0.02,
    star_font_size=8,
    disable_auto_shrink=False
):
    """
    Draw collision-aware significance lines + stars with a bounding box check
    so lines do NOT overlap bars or each other.

    - If disable_auto_shrink=True, we skip the logic that tries to shrink 
      the star font if there's no horizontal space. 
      Otherwise, we do the usual auto-shrink approach.
    """
    # Filter out comparisons with no star marker
    comparisons = [(g, c, e, p, s) for (g, c, e, p, s) in comparisons if s]

    for (gene_index, control_index, experimental_index, p_value, sig_marker) in comparisons:
        x_left  = bar_positions[gene_index] + control_index * bar_width
        x_right = bar_positions[gene_index] + experimental_index * bar_width
        x_min   = min(x_left, x_right)
        x_max   = max(x_left, x_right)
        center_x= (x_left + x_right) / 2

        # Start above the TALLEST bar of this gene
        tallest_bar_top = max_bar_per_gene[gene_index]
        line_y = tallest_bar_top + base_offset

        # Keep moving up if we collide
        while True:
            bottom_y = line_y - (cap_size * 1.5)
            top_y    = line_y + (cap_size * 3)
            new_box  = (x_min, x_max, bottom_y, top_y)
            if intersects_any(new_box, occupied_boxes):
                line_y += step
            else:
                break

        # Horizontal line
        ax.plot([x_left, x_right], [line_y, line_y], lw=1.3, c='black')
        # Caps
        ax.plot([x_left, x_left], [line_y, line_y - cap_size], lw=1.3, c='black')
        ax.plot([x_right, x_right], [line_y, line_y - cap_size], lw=1.3, c='black')

        # Decide star font size
        if disable_auto_shrink:
            effective_font_size = star_font_size
        else:
            needed_width = len(sig_marker) * star_font_size * 0.6
            available_width = (x_max - x_min)
            if needed_width > available_width > 0:
                scale_factor = available_width / needed_width
                effective_font_size = star_font_size * scale_factor
                if effective_font_size < 3:
                    effective_font_size = 3
            else:
                effective_font_size = star_font_size

        star_y = line_y + (cap_size * 0.2)
        ax.text(center_x, star_y, sig_marker,
                ha='center', va='bottom', color='black', fontsize=effective_font_size)

        # Mark bounding box
        occupied_boxes.append((x_min, x_max, bottom_y, top_y))

###############################################################################
# Figure Download Helpers
###############################################################################
def get_figure_bytes(fig, img_format="tiff", dpi=300, transparent=True):
    """
    Return the raw bytes of the figure in user-specified format, dpi, 
    and background color (transparent or white).
    """
    bg_color = 'none' if transparent else 'white'
    buf = io.BytesIO()
    fig.savefig(
        buf,
        format=img_format.lower(),
        dpi=dpi,
        facecolor=bg_color,
        edgecolor=bg_color,
        transparent=transparent
    )
    buf.seek(0)
    return buf.getvalue()

###############################################################################
# Main Streamlit App
###############################################################################
st.title('qPCR Data Analysis using Delta-Delta Ct Method')

# 1) File uploaders
uploaded_files_qpcr = st.file_uploader("Choose qPCR output files", type="txt", accept_multiple_files=True)
uploaded_files_label_map = st.file_uploader("Choose label map files", type="txt", accept_multiple_files=True)

drop_high_cp = st.checkbox("Drop rows with Cp >= 35")
cleaned_data_list = []
matched_files = {}

# 2) Matching qPCR + label map
if uploaded_files_qpcr and uploaded_files_label_map:
    for q_file in uploaded_files_qpcr:
        q_file_name = os.path.basename(q_file.name).replace("_qpcr.txt", "")
        for l_file in uploaded_files_label_map:
            l_file_name = os.path.basename(l_file.name).replace("_platemap.txt", "")
            if q_file_name == l_file_name:
                matched_files[q_file_name] = {'qpcr': q_file, 'label_map': l_file}

    for key, files in matched_files.items():
        try:
            data_qpcr = pd.read_csv(files['qpcr'], sep='\t', skiprows=1)
            columns_to_drop = [0,1,5,6,7]
            data_qpcr.drop(data_qpcr.columns[columns_to_drop], axis=1, inplace=True)

            data_label_map = pd.read_csv(files['label_map'], sep='\t')
            if 'Cell Position' in data_label_map.columns:
                data_label_map.rename(columns={'Cell Position':'Pos'}, inplace=True)

            merged = pd.merge(data_qpcr, data_label_map, on='Pos', how='left').dropna()
            if drop_high_cp:
                merged = merged[merged['Cp']<35]

            cleaned_data_list.append(merged)
            st.write(f"Cleaned Data for File Pair: {key}")
            st.dataframe(merged)

        except Exception as e:
            st.warning(f"Warning for {files['qpcr'].name}/{files['label_map'].name}: {e}")

# 3) Intersection vs. union of Genes
if cleaned_data_list:
    common_genes_set = set(cleaned_data_list[0]['Gene'].unique())
    all_genes_set    = set(cleaned_data_list[0]['Gene'].unique())
    for df in cleaned_data_list[1:]:
        common_genes_set.intersection_update(df['Gene'].unique())
        all_genes_set.update(df['Gene'].unique())

    use_common_genes_only = st.checkbox("Only use common genes across all file pairs", True)
    if use_common_genes_only:
        genes_for_selection = sorted(list(common_genes_set))
    else:
        genes_for_selection = sorted(list(all_genes_set))

    reference_genes = st.multiselect("Select Reference/Housekeeping Genes", genes_for_selection)
    target_genes    = st.multiselect("Select Target Genes", genes_for_selection)

    all_samples_union = sorted(set(s for df in cleaned_data_list for s in df['Sample'].unique()))
    control_samples = st.selectbox("Select Control Sample ID", all_samples_union)
    experimental_samples = st.multiselect("Select Experimental Sample IDs", all_samples_union)

# We'll store T-test/ANOVA in session_state
if "analysis_data" not in st.session_state:
    st.session_state["analysis_data"] = pd.DataFrame()
if "ttest_results" not in st.session_state:
    st.session_state["ttest_results"] = []
if "anova_tukey_results" not in st.session_state:
    st.session_state["anova_tukey_results"] = []

analysis_method = st.radio(
    "Select statistical test method",
    [
        "Automatic (T-test if 2 groups, ANOVA+Tukey if >2 groups)",
        "T-test only",
        "ANOVA+Tukey only"
    ],
    index=0
)

# 4) ΔΔCt Analysis
if st.button('Perform ΔΔCt Analysis') and reference_genes and target_genes:
    combined_df = pd.DataFrame()

    for df in cleaned_data_list:
        ref_df = df[df['Gene'].isin(reference_genes)]
        avg_ref = ref_df.groupby('Sample')['Cp'].mean().reset_index()
        avg_ref.rename(columns={'Cp': 'Average_Ref_Cp'}, inplace=True)

        tar_df = df[df['Gene'].isin(target_genes)]
        merged_df = pd.merge(tar_df, avg_ref, on='Sample')
        merged_df['dCp'] = merged_df['Cp'] - merged_df['Average_Ref_Cp']
        combined_df = pd.concat([combined_df, merged_df])

    # Control dCp
    ctrl_dCp = combined_df[combined_df['Sample']==control_samples].groupby('Gene')['dCp'].mean().reset_index()
    ctrl_dCp.rename(columns={'dCp':'Control_dCp'}, inplace=True)

    combined_df = pd.merge(combined_df, ctrl_dCp, on='Gene')
    combined_df['ddCp'] = combined_df['dCp'] - combined_df['Control_dCp']
    combined_df['Fold_Change'] = 2 ** (-combined_df['ddCp'])

    grouped = combined_df.groupby(['Sample','Gene'])['Fold_Change']
    combined_df['SD']  = grouped.transform('std')
    combined_df['SEM'] = grouped.transform(lambda x: np.std(x, ddof=1)/np.sqrt(len(x)))

    final_df = combined_df[combined_df['Sample'].isin([control_samples]+experimental_samples)]

    # Sort genes as user-specified
    final_df["Gene"] = pd.Categorical(
        final_df["Gene"],
        categories=target_genes,
        ordered=True
    )

    # Sort samples as user-specified
    sample_order = [control_samples] + experimental_samples
    final_df["Sample"] = pd.Categorical(
        final_df["Sample"],
        categories=sample_order,
        ordered=True
    )

    st.session_state["analysis_data"] = final_df
    st.write("Filtered DataFrame with Calculated Values (Control + Experimental):", final_df)

    # Decide T-test vs ANOVA
    n_groups = len(experimental_samples)+1
    do_ttest = False
    do_anova = False
    if analysis_method=="T-test only":
        do_ttest = True
    elif analysis_method=="ANOVA+Tukey only":
        do_anova = True
    else:
        if n_groups==2:
            do_ttest=True
        else:
            do_anova=True

    # Clear old results
    st.session_state["ttest_results"] = []
    st.session_state["anova_tukey_results"] = []

    # T-test
    if do_ttest:
        for gene in target_genes:
            ctrl_vals = final_df[
                (final_df['Sample']==control_samples) & (final_df['Gene']==gene)
            ]['Fold_Change']
            for exp in experimental_samples:
                exp_vals = final_df[
                    (final_df['Sample']==exp) & (final_df['Gene']==gene)
                ]['Fold_Change']
                if not ctrl_vals.empty and not exp_vals.empty:
                    t_stat, p_val = ttest_ind(ctrl_vals, exp_vals, equal_var=False)
                    st.session_state["ttest_results"].append({
                        'Gene': gene,
                        'Control': control_samples,
                        'Experimental': exp,
                        'T-Statistic': t_stat,
                        'p_val': p_val
                    })

    # ANOVA + Tukey
    if do_anova:
        for gene in target_genes:
            subdf = final_df[final_df['Gene']==gene]
            if subdf['Sample'].nunique()>1:
                model = ols('Fold_Change ~ C(Sample)', data=subdf).fit()
                tukey = pairwise_tukeyhsd(
                    endog=subdf['Fold_Change'],
                    groups=subdf['Sample'],
                    alpha=0.05
                )
                for row in tukey.summary().data[1:]:
                    grp1, grp2, meandiff, p_adj, lower, upper, reject = row
                    st.session_state["anova_tukey_results"].append({
                        'Gene': gene,
                        'Group1': grp1,
                        'Group2': grp2,
                        'Mean Diff': meandiff,
                        'p_val': p_adj,
                        'Reject (Significant)': reject
                    })

# 5) Show Statistical Results with Unrounded Values
if st.session_state["ttest_results"]:
    st.write("### T-Test Results")
    ttest_df = pd.DataFrame(st.session_state["ttest_results"])
    # Show unrounded on screen (15 decimals)
    st.dataframe(
        ttest_df.style.format({
            "T-Statistic": "{:.15f}",
            "p_val": "{:.15f}"
        })
    )
    # CSV export (also 15 decimals)
    csv_data_ttest = ttest_df.to_csv(index=False, float_format="%.15f")
    st.download_button(
        label="Download T-Test CSV",
        data=csv_data_ttest.encode("utf-8"),
        file_name="ttest_results.csv",
        mime="text/csv"
    )

if st.session_state["anova_tukey_results"]:
    st.write("### ANOVA + Tukey Results")
    anova_df = pd.DataFrame(st.session_state["anova_tukey_results"])
    # Show unrounded on screen
    st.dataframe(
        anova_df.style.format({
            "Mean Diff": "{:.15f}",
            "p_val": "{:.15f}"
        })
    )
    # CSV export
    csv_data_anova = anova_df.to_csv(index=False, float_format="%.15f")
    st.download_button(
        label="Download ANOVA+Tukey CSV",
        data=csv_data_anova.encode("utf-8"),
        file_name="anova_tukey_results.csv",
        mime="text/csv"
    )

# 6) Graph Customization
if cleaned_data_list:
    all_samples = sorted(set(sample for df in cleaned_data_list for sample in df['Sample'].unique()))

    palette_options = ['husl','Set1','Set2','Set3','Paired','Dark2','Accent']
    selected_palette = st.sidebar.selectbox("Choose a Categorical Color Palette", palette_options, index=0)
    chosen_palette = sns.color_palette(selected_palette, len(all_samples))

    sample_colors = {}
    sample_rename_dict = {}
    for i, s in enumerate(all_samples):
        default_color = mcolors.to_hex(chosen_palette[i])
        chosen_color = st.sidebar.color_picker(f"Color for {s}", default_color)
        sample_colors[s] = chosen_color
        new_label = st.sidebar.text_input(f"Rename legend entry for {s}", value=s)
        sample_rename_dict[s] = new_label

graph_width  = st.sidebar.slider("Graph Width", 5, 20, 10)
graph_height = st.sidebar.slider("Graph Height", 5, 20, 6)
error_bar_type=st.sidebar.selectbox("Error Bar Type", ["SD","SEM"])
bar_width     =st.sidebar.slider("Bar Width",0.1,1.0,0.4)
graph_title   =st.sidebar.text_input("Graph Title","Average Fold Change per Gene by Group")
xlabel        =st.sidebar.text_input("X-axis Label","Gene")
ylabel        =st.sidebar.text_input("Y-axis Label","Average Fold Change")
font_size     =st.sidebar.slider("Font Size (Axis/Labels)",8,20,14)
font_color    =st.sidebar.color_picker("Pick a Font Color",'#000000')
significance_level=st.sidebar.number_input("Significance Level (e.g., 0.05 for p<0.05)",
                                            0.0,1.0,0.05,0.01)
x_labels_diagonal=st.sidebar.checkbox("Rotate X-axis Labels Diagonally", True)
star_font_size = st.sidebar.slider("Star Font Size", 4, 40, 8, step=1)
disable_auto_shrink = st.sidebar.checkbox("Disable auto-shrink for star text?", value=False)
base_offset = st.sidebar.slider("Base Offset Above Bars", 0.0, 2.0, 0.2, step=0.05)
cap_size    = st.sidebar.slider("Significance Cap Size", 0.01, 0.2, 0.02, step=0.01)
collision_step=st.sidebar.slider("Collision Step Size",0.01,1.0,0.15,step=0.01)
hide_ns     = st.sidebar.checkbox("Hide 'ns' Comparisons?", value=True)

download_format= st.sidebar.selectbox("Download Format", ["TIFF","PNG"], index=0)
download_bg   = st.sidebar.selectbox("Background Color", ["Transparent","White"], index=0)
download_dpi  = st.sidebar.selectbox("DPI", [100,300,600,1200], index=1)

# 7) Generate Graph
if 'analysis_data' in st.session_state and not st.session_state['analysis_data'].empty and st.button('Generate Graph'):
    plot_df = st.session_state['analysis_data']
    if 'Fold_Change' not in plot_df.columns:
        st.warning("No 'Fold_Change' column found. Please re-check your analysis.")
    else:
        pivot_df = plot_df.pivot_table(
            index='Gene',  
            columns='Sample',
            values='Fold_Change',
            aggfunc='mean',
            observed=False
        )
        err_pivot_df = plot_df.pivot_table(
            index='Gene',
            columns='Sample',
            values=error_bar_type,
            aggfunc='mean',
            observed=False
        )

        fig, ax = plt.subplots(figsize=(graph_width, graph_height))
        num_groups = len(pivot_df.index)
        num_bars_in_group = len(pivot_df.columns)
        space_between_bars = 0.06
        space_between_groups= bar_width/2

        global bar_positions
        bar_positions = np.arange(num_groups)*(
            num_bars_in_group*bar_width + space_between_groups + space_between_bars
        )

        bar_height_lookup = {}
        occupied_boxes = []

        for gene_index, gene in enumerate(pivot_df.index):
            for sample_index, sample in enumerate(pivot_df.columns):
                x_center = bar_positions[gene_index] + sample_index*bar_width
                mean_val = pivot_df.loc[gene, sample] if not pd.isna(pivot_df.loc[gene, sample]) else 0
                err_val  = err_pivot_df.loc[gene, sample] if not pd.isna(err_pivot_df.loc[gene, sample]) else 0

                legend_label = sample_rename_dict[sample] if gene_index==0 else None

                ax.bar(
                    x_center, 
                    mean_val,
                    bar_width,
                    color=sample_colors[sample],
                    capsize=5,
                    label=legend_label
                )
                ax.errorbar(
                    x_center,
                    mean_val,
                    yerr=err_val,
                    fmt='none',
                    ecolor='black',
                    capsize=5
                )

                bar_top = mean_val+err_val
                bar_height_lookup[(gene_index,sample_index)] = bar_top

                x_min = x_center-(bar_width/2)
                x_max = x_center+(bar_width/2)
                y_min = 0
                y_max = bar_top
                occupied_boxes.append((x_min,x_max,y_min,y_max))

        max_bar_per_gene = {}
        for gene_index in range(num_groups):
            heights_for_gene = [
                bar_height_lookup.get((gene_index,s_idx),0.0)
                for s_idx in range(num_bars_in_group)
            ]
            max_bar_per_gene[gene_index] = max(heights_for_gene) if heights_for_gene else 0

        # collect comparisons
        all_comparisons = []
        # T-test
        if 'ttest_results' in st.session_state:
            ttest_df = pd.DataFrame(st.session_state['ttest_results'])
            if not ttest_df.empty and 'Gene' in ttest_df.columns:
                for _, row in ttest_df.iterrows():
                    g_name= row['Gene']
                    pval  = row.get('p_val', row.get('P-Value',1.0))
                    c_lab = row['Control']
                    e_lab = row['Experimental']

                    if (g_name in pivot_df.index) and (c_lab in pivot_df.columns) and (e_lab in pivot_df.columns):
                        gene_idx= list(pivot_df.index).index(g_name)
                        c_idx   = list(pivot_df.columns).index(c_lab)
                        e_idx   = list(pivot_df.columns).index(e_lab)

                        if pval< significance_level/100: star='****'
                        elif pval<significance_level/10: star='***'
                        elif pval<significance_level:     star='**'
                        else:                              star='ns'

                        if hide_ns and star=='ns':
                            continue

                        all_comparisons.append((gene_idx,c_idx,e_idx,pval,star))

        # ANOVA
        if 'anova_tukey_results' in st.session_state:
            anova_df= pd.DataFrame(st.session_state['anova_tukey_results'])
            if not anova_df.empty and 'Gene' in anova_df.columns:
                for _, row in anova_df.iterrows():
                    g_name= row['Gene']
                    pval  = row.get('p_val', row.get('p-value',1.0))
                    g1    = row['Group1']
                    g2    = row['Group2']

                    if (g_name in pivot_df.index) and (g1 in pivot_df.columns) and (g2 in pivot_df.columns):
                        gene_idx= list(pivot_df.index).index(g_name)
                        idx1    = list(pivot_df.columns).index(g1)
                        idx2    = list(pivot_df.columns).index(g2)

                        if pval< significance_level/100: star='****'
                        elif pval<significance_level/10:  star='***'
                        elif pval<significance_level:      star='**'
                        else:                              star='ns'

                        if hide_ns and star=='ns':
                            continue

                        all_comparisons.append((gene_idx,idx1,idx2,pval,star))

        is_transparent = (download_bg=="Transparent")

        # draw significance lines
        draw_significance_indicator_collision_aware(
            ax=ax,
            comparisons=all_comparisons,
            bar_positions=bar_positions,
            bar_width=bar_width,
            bar_height_lookup=bar_height_lookup,
            occupied_boxes=occupied_boxes,
            max_bar_per_gene=max_bar_per_gene,
            base_offset=base_offset,
            step=collision_step,
            cap_size=cap_size,
            star_font_size=star_font_size,
            disable_auto_shrink=disable_auto_shrink
        )

        ax.set_xticks(
            bar_positions + (num_bars_in_group*bar_width/2) - (bar_width/2)
        )
        ax.set_xticklabels(pivot_df.index, fontsize=font_size)
        ax.set_xlabel(xlabel, fontsize=font_size, color=font_color)
        ax.set_ylabel(ylabel, fontsize=font_size, color=font_color)
        ax.set_title(graph_title, fontsize=font_size+2, color=font_color)

        if x_labels_diagonal:
            plt.xticks(rotation=45, ha='right')

        leg= ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
        if is_transparent and leg is not None:
            leg.get_frame().set_facecolor('none')
            leg.get_frame().set_alpha(0.0)

        sns.despine()
        plt.tight_layout(rect=[0,0,0.85,1])
        st.pyplot(fig)

        # Download figure
        fig_bytes= get_figure_bytes(
            fig,
            img_format= download_format.lower(),
            dpi= download_dpi,
            transparent= is_transparent
        )
        mime_type= "image/tiff" if download_format.lower()=="tiff" else "image/png"
        st.download_button(
            label=f"Download Graph ({download_format.upper()}, DPI={download_dpi})",
            data= fig_bytes,
            file_name= f"plot.{download_format.lower()}",
            mime= mime_type
        )
else:
    st.warning("No 'Fold_Change' column found. Please re-check your analysis.")
