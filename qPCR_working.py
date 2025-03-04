import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
import io
import base64
fig = None
ax = None

# Function to draw significance indicators
def draw_significance_indicator(ax, gene_index, control_index, experimental_index, max_y, bar_width, y_offset=0.1, sig_marker='*', cap_length=0.02):
    # Calculate the x positions for the control and experimental bars
    control_bar_x = bar_positions[gene_index] + control_index * bar_width
    experimental_bar_x = bar_positions[gene_index] + experimental_index * bar_width

    # Calculate the center position for the significance indicator
    center_x = (control_bar_x + experimental_bar_x) / 2

    # Calculate the y position for the significance indicator
    indicator_y = max_y + y_offset

    # Draw the horizontal line of the significance indicator
    ax.plot([center_x - bar_width / 2, center_x + bar_width / 2], [indicator_y, indicator_y], lw=1.0, c='black')

    # Draw the short vertical lines at the ends of the horizontal line
    ax.plot([center_x - bar_width / 2, center_x - bar_width / 2], [indicator_y, indicator_y - cap_length], lw=-1.0, c='black')
    ax.plot([center_x + bar_width / 2, center_x + bar_width / 2], [indicator_y, indicator_y - cap_length], lw=-1.0, c='black')

    # Add the significance marker text
    ax.text(center_x, indicator_y + cap_length - 0.02, sig_marker, ha='center', va='bottom', color='black', fontsize=12)

# Function to download a figure
def get_image_download_link(fig, filename="plot.tiff", text="Download Plot as TIFF"):
    buf = io.BytesIO()
    fig.savefig(buf, format='tiff')
    buf.seek(0)
    img_bytes = buf.read()
    b64 = base64.b64encode(img_bytes).decode()
    href = f'<a href="data:image/tiff;base64,{b64}" download="{filename}">{text}</a>'
    return href

# Function to download a DataFrame as Excel
def get_table_download_link_excel(df, filename="data.xlsx", text="Download analysis as Excel"):
    towrite = io.BytesIO()
    df.to_excel(towrite, index=False)  # Write DataFrame to the buffer
    towrite.seek(0)  # Reset the pointer
    b64 = base64.b64encode(towrite.read()).decode()  # Convert to base64
    href = f'<a href="data:application/octet-stream;base64,{b64}" download="{filename}">{text}</a>'
    return href

# Function to download a DataFrame as TXT
def get_table_download_link_txt(df, filename="data.txt", text="Download analysis as TXT"):
    towrite = io.StringIO()
    df.to_csv(towrite, index=False, sep='\t')  # Write DataFrame to the buffer
    towrite.seek(0)
    b64 = base64.b64encode(towrite.getvalue().encode()).decode()  # Convert to base64
    href = f'<a href="data:file/txt;base64,{b64}" download="{filename}">{text}</a>'
    return href

# Title of the app
st.title('qPCR Data Analysis using Delta-Delta Ct Method')

# File uploaders
uploaded_file_qpcr = st.file_uploader("Choose a qPCR output file", type="txt")
uploaded_file_label_map = st.file_uploader("Choose the label map file", type="txt")
drop_high_cp = st.checkbox("Drop rows with Cp >= 35")

# Initialize empty DataFrame and ttest_results
cleaned_data = pd.DataFrame()
ttest_results = []

# Select box for choosing the sample ID column
sample_id_column = st.selectbox("Choose the column for Sample IDs", ["Name", "Sample"])

# Process files only if both are uploaded
if uploaded_file_qpcr and uploaded_file_label_map:
    try:
        # Read the qPCR file
        data_qpcr = pd.read_csv(uploaded_file_qpcr, sep='\t', skiprows=1)

        # Drop columns by index
        columns_to_drop = [0, 1, 5, 6, 7]
        data_qpcr.drop(data_qpcr.columns[columns_to_drop], axis=1, inplace=True)

        # Read the label map file
        data_label_map = pd.read_csv(uploaded_file_label_map, sep='\t')
        if 'Cell Position' in data_label_map.columns:
            data_label_map.rename(columns={'Cell Position': 'Pos'}, inplace=True)

        # Merge and clean data
        merged_data = pd.merge(data_qpcr, data_label_map, on='Pos', how='left').dropna()
        cleaned_data = merged_data[merged_data['Cp'] < 35] if drop_high_cp else merged_data

        st.write("Cleaned Data:", cleaned_data)

        # Create a dictionary for mapping sample names to colors after defining sample_id_column
        unique_samples = cleaned_data[sample_id_column].unique()
    except Exception as e:
        st.error(f"Error processing files: {e}")
else:
    unique_samples = []

# User selection of colors for each sample ID
sample_colors = {sample: st.sidebar.color_picker(f"Color for {sample}", '#19A09C') for sample in unique_samples}

# Rest of your analysis and plotting code goes here...


# Rest of the analysis logic
if not cleaned_data.empty:
    # User selections for reference genes, target genes, control and experimental samples
    reference_genes = st.multiselect("Select Reference/Housekeeping Genes", options=cleaned_data['Gene'].unique(), default=None)
    target_genes = st.multiselect("Select Target Genes (and order)", options=cleaned_data['Gene'].unique(), default=None)
    control_samples = st.selectbox("Select Control Sample ID", options=cleaned_data[sample_id_column].unique())
    experimental_samples = st.multiselect("Select Experimental Sample IDs", options=cleaned_data[sample_id_column].unique(), default=None)

# Button to perform the analysis
# Button to perform the analysis
if st.button('Perform ΔΔCt Analysis') and reference_genes and target_genes:
    # Step 1: Calculate Average Reference Cp Value per Sample ID
    ref_data = cleaned_data[cleaned_data['Gene'].isin(reference_genes)]
    
    # Calculate the average Cp of the reference genes for each sample
    average_ref_cp = ref_data.groupby(sample_id_column)['Cp'].mean().reset_index()
    average_ref_cp.rename(columns={'Cp': 'Average_Ref_Cp'}, inplace=True)

    # Step 2: Merge the average reference Cp with the target gene data
    target_data = cleaned_data[cleaned_data['Gene'].isin(target_genes)]
    merged_data = pd.merge(target_data, average_ref_cp, on=sample_id_column)

    # Step 3: Calculate dCp for each target gene (Cp - Average_Ref_Cp for the same sample)
    merged_data['dCp'] = merged_data['Cp'] - merged_data['Average_Ref_Cp']

    # Step 4: Calculate the average dCp for the control samples per gene
    control_dCp_avg = merged_data[merged_data[sample_id_column] == control_samples].groupby('Gene')['dCp'].mean().reset_index()
    control_dCp_avg.rename(columns={'dCp': 'Control_Ave_dCp'}, inplace=True)

    # Step 5: Merge the control average dCp into the dataset
    merged_data = pd.merge(merged_data, control_dCp_avg, on='Gene')

    # Step 6: Subtract the control average dCp from all samples to calculate ddCp
    merged_data['ddCp'] = merged_data['dCp'] - merged_data['Control_Ave_dCp']

    # Step 7: Calculate the Fold Change
    merged_data['Fold_Change'] = 2 ** (-merged_data['ddCp'])

    # Step 8: Calculate SD and SEM for Fold Change
    grouped_data = merged_data.groupby([sample_id_column, 'Gene'])
    merged_data['SD'] = grouped_data['Fold_Change'].transform('std')
    merged_data['SEM'] = grouped_data['Fold_Change'].transform(lambda x: np.std(x, ddof=1) / np.sqrt(len(x)))

    # Step 9: Filter the data to only include control and experimental samples
    filtered_data = merged_data[merged_data[sample_id_column].isin([control_samples] + experimental_samples)]

    # Display the DataFrame with filtered calculated ddCp and Fold Change values
    st.write("Filtered DataFrame with Calculated Values (Control and Experimental Samples):", filtered_data)

    # Perform t-tests for each gene and save the results
    ttest_results.clear()
    for gene in target_genes:
       for sample in experimental_samples:
            if sample != control_samples:
                control_fc = merged_data[(merged_data[sample_id_column] == control_samples) & (merged_data['Gene'] == gene)]['Fold_Change']
                experimental_fc = merged_data[(merged_data[sample_id_column] == sample) & (merged_data['Gene'] == gene)]['Fold_Change']
                if not control_fc.empty and not experimental_fc.empty:
                    t_stat, p_val = ttest_ind(control_fc, experimental_fc, equal_var=False)
                    ttest_results.append({'Gene': gene, 'Control': control_samples, 'Experimental': sample, 'T-Statistic': t_stat, 'P-Value': p_val})

    # Save filtered analysis data and t-test results in session state
    st.session_state['analysis_data'] = filtered_data
    st.session_state['ttest_results'] = ttest_results


       

# Displaying t-test results
if 'ttest_results' in st.session_state:
    ttest_df = pd.DataFrame(st.session_state['ttest_results'])
    st.write("T-Test Results:", ttest_df)



# Graph customization options
graph_width = st.sidebar.slider("Graph Width", min_value=5, max_value=20, value=10)
graph_height = st.sidebar.slider("Graph Height", min_value=5, max_value=20, value=6)
error_bar_type = st.sidebar.selectbox("Error Bar Type", ["SD", "SEM"])
graph_title = st.sidebar.text_input("Graph Title", "Average Fold Change per Gene by Group")
xlabel = st.sidebar.text_input("X-axis Label", "Gene")
ylabel = st.sidebar.text_input("Y-axis Label", "Average Fold Change")
font = st.sidebar.selectbox("Font", plt.rcParams['font.family'])
font_size = st.sidebar.slider("Font Size", min_value=8, max_value=20, value=12)
font_color = st.sidebar.color_picker("Pick a Font Color", '#000000')
significance_level = st.sidebar.number_input("Significance Level (e.g., 0.05 for P<0.05)", min_value=0.0, max_value=1.0, value=0.05, step=0.01)

# Define bar_width
bar_width = st.sidebar.slider("Bar Width", min_value=0.1, max_value=1.0, value=0.4)

# Button to generate the graph
if 'analysis_data' in st.session_state and st.button('Generate Graph'):
    plot_df = pd.DataFrame(st.session_state['analysis_data'])

    if 'Fold_Change' in plot_df.columns:
        # Setup the matplotlib figure and axes
        fig, ax = plt.subplots(figsize=(graph_width, graph_height))

        # Create pivot tables for Fold Change and Error (SD or SEM)
        pivot_df = plot_df.pivot_table(index='Gene', columns=sample_id_column, values='Fold_Change', aggfunc='mean')
        error_pivot_df = plot_df.pivot_table(index='Gene', columns=sample_id_column, values=error_bar_type, aggfunc='mean')

        # Define the number of groups and number of bars in each group
        num_groups = len(pivot_df.index)  # Number of genes
        num_bars_in_group = len(pivot_df.columns)  # Number of samples (control and experimental)
        space_between_bars = 0.01  # adjust this value as needed
        # Space between groups
        space_between_groups = bar_width / 2

        # Calculate the positions for each group's bars
        bar_positions = np.arange(num_groups) * (num_bars_in_group * bar_width + space_between_groups)

        # Plot the bars for each sample
        for i, sample in enumerate(pivot_df.columns):
            offsets = bar_positions + i * bar_width
            color = sample_colors.get(sample, '#278C88')  # Default color if not found
            ax.bar(offsets, pivot_df[sample], bar_width, yerr=error_pivot_df[sample], label=sample, color=color, capsize=5)

        # Draw significance indicators
        if 'ttest_results' in st.session_state:
            ttest_df = pd.DataFrame(st.session_state['ttest_results'])
            for _, row in ttest_df.iterrows():
                if row['P-Value'] < significance_level:
                    gene_index = pivot_df.index.get_loc(row['Gene'])
                    control_index = list(pivot_df.columns).index(row['Control'])
                    experimental_index = list(pivot_df.columns).index(row['Experimental'])
                    max_y = max(
                        pivot_df.loc[row['Gene'], row['Control']] + error_pivot_df.loc[row['Gene'], row['Control']],
                        pivot_df.loc[row['Gene'], row['Experimental']] + error_pivot_df.loc[row['Gene'], row['Experimental']]
                    )
                    draw_significance_indicator(ax, gene_index, control_index, experimental_index, max_y, bar_width)

        # Setting graph attributes
        ax.set_xticks(bar_positions + (num_bars_in_group * bar_width / 2) - (bar_width / 2))
        ax.set_xticklabels(pivot_df.index)
        ax.set_xlabel(xlabel, fontsize=font_size, color=font_color)
        ax.set_ylabel(ylabel, fontsize=font_size, color=font_color)
        ax.set_title(graph_title, fontsize=font_size + 2, color=font_color)
        ax.legend()

        # Finalize and display the figure
        st.pyplot(fig)
    else:
        st.error("'Fold_Change' column not found in the DataFrame.")



# Place these lines where you want to display the download links in your Streamlit app
if fig is not None:
    st.markdown(get_image_download_link(fig), unsafe_allow_html=True)
    analysis_data_df = pd.DataFrame(st.session_state['analysis_data'])
    st.markdown(get_table_download_link_excel(analysis_data_df), unsafe_allow_html=True)
    st.markdown(get_table_download_link_txt(analysis_data_df), unsafe_allow_html=True)


