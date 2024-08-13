import pandas as pd

def read_color_map(path):
    cmap_df = pd.read_excel(path)
    cmap_df.pop("Color")
    cmap_dict = pd.Series(cmap_df.Hex.values, index=cmap_df.Population).to_dict()
    return cmap_dict