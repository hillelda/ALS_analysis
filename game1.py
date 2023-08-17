# Hillel Darshan 8/8/23
import numpy as np
import pandas as pd



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
    df = load_data("./small_data_preproccesed.csv")
    sources_list = df["tRNA Source"].unique()
    sources_indexes = np.zeros(len(sources_list))
    sources_indexes = np.vstack((sources_list, sources_indexes)).T
    data = df.to_numpy()
    for source in sources_list:
        pass



    b = 1 # for debugging