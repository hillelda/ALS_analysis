# Hillel Darshan 8/17/23
import numpy as np
import pandas as pd
import plotly.graph_objects as go

def draw_volcano_plots(df_list, full_add_list, ds_name_list):
    for i, df in enumerate(df_list):
        # draw volcano plots
        df['binary PValue color'] = df['FDR'].apply(lambda x: "red" if x > -np.log10(0.05) else "blue")

        pl = go.Figure(
            go.Scatter(x=df["log2 FC"],
                       y=df["-log 10 FDR"],
                       mode="markers",
                       marker=dict(color=df['binary PValue color'])))
        pl.update_layout(title="ALS " + ds_name_list[i],
                         xaxis_title="Log2 FC",
                         yaxis_title="-Log(10) FDR")

        # pl.show()
        pl.write_html(file=full_add_list[i] + "_plot.html")


def pre_process_als(files_location, data_files, df_list, full_add_list, ds_name_list):
    for cur_file in data_files:
        full_add = files_location + cur_file
        df = pd.read_csv(full_add)
        split_name = cur_file.split('_')
        ds_name_list.append(split_name[0] + " " + split_name[1])

        # refactor logFC from ln to log2
        df["log2 FC"] = np.log2(np.exp(df["logFC"]))
        # refactor FDR to -log10 FDR
        df["-log 10 FDR"] = -np.log10(df["FDR"])

        # Add one-hot columns for significance and sign
        df['sig up als'] = (df['-log 10 FDR'] < 0.05) & (df['logFC'] >= 0)
        df['sig down als'] = (df['-log 10 FDR'] < 0.05) & (df['logFC'] < 0)
        df['unsig up als'] = (df['-log 10 FDR'] >= 0.05) & (df['logFC'] >= 0)
        df['unsig down als'] = (df['-log 10 FDR'] >= 0.05) & (df['logFC'] < 0)

        df.set_index('transcript', inplace=True)
        df_list.append(df)
        full_add_list.append(full_add)

def pre_process_BoNT(files_location_address, cur_file):
    full_add = files_location_address + cur_file
    Bo_NT_df = pd.read_csv(full_add)
    # refactor logFC from ln to log2
    Bo_NT_df["log2 FC"] = np.log2(np.exp(Bo_NT_df["logFC"]))
    # refactor FDR to -log10 FDR
    Bo_NT_df["-log 10 FDR"] = -np.log10(Bo_NT_df["FDR"])

    # Add one-hot columns for significance and sign
    Bo_NT_df['sig up bont'] = (Bo_NT_df['-log 10 FDR'] < 0.05) & (Bo_NT_df['logFC'] >= 0)
    Bo_NT_df['sig down bont'] = (Bo_NT_df['-log 10 FDR'] < 0.05) & (Bo_NT_df['logFC'] < 0)
    Bo_NT_df['unsig up bont'] = (Bo_NT_df['-log 10 FDR'] >= 0.05) & (Bo_NT_df['logFC'] >= 0)
    Bo_NT_df['unsig down bont'] = (Bo_NT_df['-log 10 FDR'] >= 0.05) & (Bo_NT_df['logFC'] < 0)

    Bo_NT_df.drop(Bo_NT_df.columns[0], axis=1, inplace=True)
    Bo_NT_df.set_index('transcript', inplace=True)
    return Bo_NT_df


def cross_sig_dfs(combined_df):
    combined_df['sig up in both'] = (combined_df['sig up als']) & (combined_df['sig up bont'])
    combined_df['sig down in both'] = (combined_df['sig down als']) & (combined_df['sig down bont'])
    combined_df['sig down als & sig up bont'] = (combined_df['sig down als']) & (combined_df['sig up bont'])
    combined_df['sig up als & sig down bont'] = (combined_df['sig up als']) & (combined_df['sig down bont'])

    return combined_df


if __name__ == '__main__':
    working_dir = "../ALS_analysis/working data - ALS/"
    data_files = ["GSE203170_miRNA_SgGens_output.csv", "GSE203170_tRF_SgGens_output.csv", "GSE94888_miRNA_SgGens_output.csv", "GSE94888_tRF_SgGens_output.csv"]
    df_list = []
    full_add_list = []
    ds_name_list = []
    pre_process_als(working_dir, data_files, df_list, full_add_list, ds_name_list)
    draw_volcano_plots(df_list, full_add_list, ds_name_list)

    bont_pd = pre_process_BoNT(working_dir, "DE_tRFs Arik BoNT_A 16.8.23.csv")

    meta_203170 = pd.read_csv(working_dir + 'GSE203170_tRF_meta.csv')
    meta_94888 = pd.read_csv(working_dir + 'GSE94888_tRF_meta.csv')

    combined_df_1 = pd.concat([df_list[1], bont_pd], axis=1, join='outer', sort=False)
    combined_df_2 = pd.concat([df_list[3], bont_pd], axis=1, join='outer', sort=False)

    combined_df_1 = cross_sig_dfs(combined_df_1)
    combined_df_2 = cross_sig_dfs(combined_df_2)

    # add meta data
    column_to_drop = ['isSg']
    combined_df_1.columns = ['trf'] + list(combined_df_1.columns[1:])
    combined_df_1 = combined_df_1.merge(meta_203170, on='trf' , how='left')
    combined_df_1.drop(column_to_drop, axis=1, inplace=True)

    combined_df_2.columns = ['trf'] + list(combined_df_2.columns[1:])
    combined_df_2 = combined_df_1.merge(meta_94888, on='trf', how='left')
    # combined_df_2.drop(column_to_drop, axis=1, inplace=True)

    combined_df_1.to_csv(working_dir + 'GSE203170_ALS_BoNT_analysis_output.csv')
    combined_df_2.to_csv(working_dir + 'GSE94888_ALS_BoNT_analysis_output.csv')

    print("All Done")