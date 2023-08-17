# Hillel Darshan 8/17/23
import numpy as np
import pandas as pd
import plotly.graph_objects as go



def load_data(filename: str) -> pd.DataFrame:
    """
    Load cvs file to pd frame
    ----------
    filename: str
        Path to house prices dataset

    Returns
    -------
    Design matrix and response vector (Temp)
    """
    df = pd.read_csv(filename)
    return df













if __name__ == '__main__':
    # data_files = ["C:/Users/H/PycharmProjects/ALS_analysis/working data - ALS/GSE203170_tRF_SgGens_output.csv"]
    data_files = ["../ALS_analysis/working data - ALS/GSE203170_tRF_SgGens_output.csv"]
    df = pd.read_csv(data_files[0]) # todo: make it work for all data sets

    # sources_list = df["transcript"].unique()

    # refactor logFC from ln to log2
    df["logFC"] = np.log2(np.exp(df["logFC"]))
    # refactor FDR to -log10 FDR
    df["FDR"] = -np.log10(df["FDR"])

    pl = go.Figure(go.Scatter(x=df["logFC"],
                              y=df["FDR"],
                              mode="markers",
                              marker=dict(color="blue")))
    pl.update_layout(title="ALS " + "SOMTHING",
                     xaxis_title="Log2 ????",
                     yaxis_title="Log(-10) ????")

    pl.show()




    # sources_indexes = np.zeros(len(sources_list))
    # sources_indexes = np.vstack((sources_list, sources_indexes)).T
    # data = df.to_numpy()
    # for source in sources_list:
    #     pass



    b = 1 # for debugging