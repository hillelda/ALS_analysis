# Hillel Darshan 8/17/23
import numpy as np
import pandas as pd
import plotly.graph_objects as go


if __name__ == '__main__':
    files_location = "../ALS_analysis/working data - ALS/"
    data_files = ["GSE203170_miRNA_SgGens_output.csv", "GSE203170_tRF_SgGens_output.csv", "GSE94888_miRNA_SgGens_output.csv", "GSE94888_tRF_SgGens_output.csv"]
    for cur_file in data_files:
        full_add = files_location + cur_file
        df = pd.read_csv(full_add)
        split_name = cur_file.split('_')
        data_set_name = split_name[0] + " " + split_name[1]

        # refactor logFC from ln to log2
        df["logFC"] = np.log2(np.exp(df["logFC"]))
        # refactor FDR to -log10 FDR
        df["FDR"] = -np.log10(df["FDR"])

        # add column for color significance
        df['binary PValue'] = df['PValue'].apply(lambda x: "red" if x < 0.05 else "blue")

        # draw volcano plots
        pl = go.Figure(
            go.Scatter(x=df["logFC"],
                                  y=df["FDR"],
                                  mode="markers",
                                  marker=dict(color=df['binary PValue'])))
        pl.update_layout(title="ALS " + data_set_name,
                         xaxis_title="Log2 FC",
                         yaxis_title="-Log(10) FDR")

        pl.show()
        pl.write_html(file=full_add + "_plot")



    b = 1 # for debugging