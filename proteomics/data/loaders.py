"""Data read-in: raw Progenesis CSV and PLGS Excel files -> raw DataFrames."""
import pandas as pd
import streamlit as st


def load_proteomics_data(file_input):
    try:
        if hasattr(file_input, 'seek'):
            file_input.seek(0) #every rerun in streamlit the pointer should go to the start of the file
        df = pd.read_csv(file_input, sep=',', encoding='utf-8')
        st.info(f"Data loaded succesfully! Shape: {df.shape[0]} rows, {df.shape[1]} columns.")    

        return df
    except Exception as e:
        st.error(f"Error loading file: {e}")
        raise


def load_plgs_data(file_input):
    #Read every sheet of a PLGS excel workbook. The fraction number is the position of the sheet in the workbook 
    #(first sheet = fraction 1) and combine it into a dataframe
    try:
        if hasattr(file_input, 'seek'):
            file_input.seek(0) #every rerun in streamlit the pointer should go to the start of the file
        sheets = pd.read_excel(file_input, sheet_name=None)  # dict, keeps the sheet order
    except Exception as e:
        st.error(f"Error loading file: {e}")
        raise

    frames = []
    for position, (sheet_name, sheet) in enumerate(sheets.items(), start=1):
        if sheet.empty:
            st.warning(f"Sheet '{sheet_name}' (fraction {position}) is empty and is skipped.")
            continue
        sheet = sheet.copy()
        sheet['Fraction'] = position
        sheet['Sheet'] = sheet_name
        frames.append(sheet)

    if not frames:
        st.error("No data found in any sheet of this workbook.")
        raise ValueError("No data found in any sheet.")

    df = pd.concat(frames, ignore_index=True)
    st.info(f"Data loaded succesfully! {len(frames)} sheet(s) = fractions, {df.shape[0]} rows.")
    return df
