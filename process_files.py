import pandas as pd
import os
import scimap as sm

def clean_file(dir_path, dataset_path, metadata_path, insert_area, sample_col, metadata_info, mkdir=False):

    path_split = dataset_path.split("\\")
    file_name = path_split[-1]

    data_df = pd.read_csv(dataset_path, index_col=False, on_bad_lines="warn")

    # remove indices column
    data_df.drop(data_df.columns[data_df.columns.str.contains('unnamed', case=False)], axis=1, inplace=True)

    # insert area column right after location x and y columns and before parent tumor mask column
    move_col = data_df.pop("Area")
    y_col_idx = data_df.columns.get_loc(insert_area)
    move_col.name = "cellArea"
    data_df.insert(y_col_idx, move_col.name, move_col)

    # edit other column names
    data_df.rename(columns={'ObjectNumber': 'CellID', 'Location_Center_X': 'X_centroid', 'Location_Center_Y': 'Y_centroid', 'class': 'phenotype'}, inplace=True)
    
    # add in patient details from metadata - patient ID, subtype, TMA area
    roi = get_roi(dataset_path)
    data_df = add_metadata_details(data_df, metadata_path, roi, sample_col, metadata_info)

    if data_df.empty:
        return -1

    # create directory in the folder of the source data
    mkdir_path = os.path.join(dir_path, "formatted_files")
    if mkdir == True:
        os.mkdir(mkdir_path, 0o777)

    # assumes file name has 1 underscore, and removes prefix before _
    file_name = file_name.split("_")[-1]

    write_to = os.path.join(mkdir_path, "formatted_" + file_name)
    data_df.to_csv(write_to, index=False)
    return write_to

def get_roi(path):
    # assumes that the sample id containing slide/ROI is in the last part of the file name, after the last underscore
    file = path.split("\\")[-1]
    # includes_roi = file.split("_")[-1]
    # includes_roi = includes_roi.replace(".csv", "")
    includes_roi = file.replace(".csv", "")
    includes_roi = includes_roi.replace("im", "") # im + Sample_ID

    # idx = includes_roi.find("ROI") # this code is for if the sample id in the metadata starts with ROI and the file name does not
    # roi = includes_roi[idx:]
    roi = includes_roi
    return roi

def add_metadata_details(data_df, metadata_path, roi, sample_col, metadata_info):
    '''
    add metadata information to the dataframe being worked with to eventually save in anndata object
    '''
    # meta_df = pd.read_csv(metadata_path)
    meta_df = pd.read_excel(metadata_path)

    # have to match metadata info based on ROI that we are currently on
    selected_row = meta_df.loc[meta_df[sample_col] == roi]
    if selected_row.empty:
        print(f"Missing {roi} in metadata. File couldn't be processed - skipped.")
        return selected_row

    # this assumes that the metadata is exactly aligned with the data csv file
    # insert patient columns at end of sheet
    for meta, col_name in metadata_info.items():
        data_df.loc[:, col_name] = selected_row[meta].values[0]

    return data_df


def create_anndata(fn_list, name):
    # all defaults should be corrected
    adata = sm.pp.mcmicro_to_scimap(feature_table_path=fn_list)
    name = name + ".h5ad"
    adata.write(name)

if __name__ == "__main__":
    
    # parent directory of where the files are located
    # metadata should be a .csv
    # dir_path = r"Z:\Multiplex_IHC_studies\Summer_Interns\2024\DATA\DP20\DP20\Classified_CSV"
    dir_path = r"Z:\Multiplex_IHC_studies\Eric_Berens\HuBrca_TMA_mIHC\DataAnalysis\CSV\Classified_CSV"
    metadata_path = r"Z:\Multiplex_IHC_studies\Eric_Berens\HuBrca_TMA_mIHC\DataAnalysis\HuBrCa_TMA_Metadata.xlsx"
    insert_area = "class" # here, have the name of the column that is right after location x and y; was Parent_TumorMask for HTAN
    sample_col = "Sample_ID" # was ROI for HTAN
    # timepoint, ER, grade, stage, age, Days_From_Tx, Tx, Tissue_ID, Subject_ID (non-unique), Sample_ID (unique)
    metadata_info = {
        "Sample_ID": "Sample_ID",
        "Core_ID": "Core_ID",
        "Tissue_ID": "Tissue_ID",
        "Subject_ID": "Subject_ID",
        "Replicate": "Replicate",
        "ER": "ER",
        "Timepoint": "Timepoint",
        "Area": "sampleArea",
        "Tx": "Tx", 
        "Timepoint": "timepoint", 
        "Grade": "Grade",
        "Stage": "Stage",
        "Age": "Age",
        "Days_From_Tx": "Days_From_Tx",
        "HER2_FISH": "HER2_FISH"
    }

    # metadata_info = {"Area": "sampleArea", "Tx": "Tx", "sampleName": "Tx_Timepoint", "Replicate": "Replicate", "Tumor_ID": "Tumor_ID"}
    # -1 instead of NA for metadata
    mkdir = True
    fn_list = []
    for file in os.listdir(dir_path):
        fn = os.path.join(dir_path, file)
        # to get a representative name for the anndata object
        
        # check that file found is a file
        if os.path.isfile(fn) == False:
            raise ValueError("Found something that is not a file. Ensure everything in the given directory is a file and retry.")

        # clean file and format
        if mkdir == True:
            fn_path = clean_file(dir_path, fn, metadata_path, insert_area, sample_col, metadata_info, True)
            mkdir = False
        else:
            fn_path = clean_file(dir_path, fn, metadata_path, insert_area, sample_col, metadata_info)
        if fn_path != -1:
            fn_list.append(fn_path)
        print("finished", file)
    
    anndata_name = "IM_full_data" # function will add the .h5ad
    create_anndata(fn_list, anndata_name)
