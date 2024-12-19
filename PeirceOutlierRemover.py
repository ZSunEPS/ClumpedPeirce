# If the required packages are not installed, uncomment the following lines to install packages.
# !pip install pandas numpy openpyxl

import pandas as pd
import numpy as np

# Load Peirce's Criterion Table
peircesTable = pd.read_csv("Peirces Table.csv")

# Read data exported from Easotope. Replace "data.csv" with the path of your .csv file.
rawExport = pd.read_csv("Data.csv")

# Create a DataFrame with desired columns and initialize outlier indicators
df = pd.DataFrame({
    "Sample": rawExport["Easotope Name"],
    "d13C": rawExport["d13C VPDB (Final)"],
    "d18O": rawExport["d18O VPDB (Final)"],
    "D47": rawExport["D47 CDES (Final)"],
    "d13C_Outlier": 0,
    "d18O_Outlier": 0,
    "D47_Outlier": 0,
    "Total_Outlier": 0
})

# Drop rows containing missing data
df.dropna(inplace=True)


gdf = df.groupby("Sample") # Group data by "Sample"
vdf = [group for name, group in gdf]  # A vector of subdataframes

# Create a summary table for the evaluated data
summary = pd.DataFrame({
    "Sample": [name for name, group in gdf],
    "n": 0,
    "d13C": 0.0,
    "d13C_SE": 0.0,
    "d18O": 0.0,
    "d18O_SE": 0.0,
    "D47": 0.0,
    "D47_SE": 0.0
})

# Function to perform Peirce's outlier detection
def PeirceOutlierDetect(Obs, Outl, Table, d):
    n = len(Obs)
    if n < 2:  # Not enough data for analysis
        Outl[:] = 0
        return
    
    Obs_mean = round(np.mean(Obs), d)
    Obs_SD = round(np.std(Obs), d + 1)
    Obs_Deviation = np.abs(Obs - Obs_mean)
    LastOutlierDetected = -1
    NewOutlierDetected = 0
    Step = 0
    while NewOutlierDetected > LastOutlierDetected:
        LastOutlierDetected = NewOutlierDetected
        Step += 1
        R = Table.iloc[n - 3, Step -1]
        Obs_MaxDeviation = round(Obs_SD * R, d + 1)
        Outl[:] = Obs_Deviation > Obs_MaxDeviation  # Mark outliers
        NewOutlierDetected = np.sum(Outl)

# Process data per group (or sample)
for k in range(len(vdf)):
    PeirceOutlierDetect(vdf[k]["d13C"].values, vdf[k]["d13C_Outlier"].values, peircesTable, 2)
    PeirceOutlierDetect(vdf[k]["d18O"].values, vdf[k]["d18O_Outlier"].values, peircesTable, 2)
    PeirceOutlierDetect(vdf[k]["D47"].values, vdf[k]["D47_Outlier"].values, peircesTable, 3)

    # Total outlier flag
    vdf[k]["Total_Outlier"] = vdf[k]["d13C_Outlier"] + vdf[k]["d18O_Outlier"] + vdf[k]["D47_Outlier"]
    
    # Filter non-outliers and calculate statistics
    cleanData = vdf[k][vdf[k]["Total_Outlier"] == 0]
    m = len(cleanData)
    summary.loc[k, "n"] = m
    summary.loc[k, "d13C"] = round(np.mean(cleanData["d13C"]), 2)
    summary.loc[k, "d18O"] = round(np.mean(cleanData["d18O"]), 2)
    summary.loc[k, "D47"] = round(np.mean(cleanData["D47"]), 3)
    summary.loc[k, "d13C_SE"] = round(np.std(cleanData["d13C"]) / np.sqrt(m), 2)
    summary.loc[k, "d18O_SE"] = round(np.std(cleanData["d18O"]) / np.sqrt(m), 2)
    summary.loc[k, "D47_SE"] = round(np.std(cleanData["D47"]) / np.sqrt(m), 3)

# The Δ47-Temperature equation is from Anderson et al. (2021)
def Ae21Temperarure(D):
    return np.sqrt(0.0391 * 10**6 / (D - 0.154)) - 273.15

def Ae21TemperarureSE(D, DSE):
    return np.abs(np.sqrt(0.0391) * 10**3 * (-1/2) * (D - 0.154)**(-3/2) * DSE)

summary["Temperature"] = np.round(Ae21Temperarure(summary["D47"]), 0).astype(int)
summary["Temperature_SE"] = np.round(Ae21TemperarureSE(summary["D47"], summary["D47_SE"]), 0).astype(int)

# Convert VPDB to VSMOW
summary["d18Oc (VSMOW)"] = np.round(summary["d18O"] * 1.03091 + 30.91, 2)

# The δ18O fractionation equation is from Kim and O'Neil (1997)
def WaterOxygenIsotopicComposition(d18Ocarb, T):
    return (1000 + d18Ocarb) / (np.exp((18030 / (T + 273.15) - 32.42) / 1000)) - 1000

summary["d18Ow (VSMOW)"] = np.round(WaterOxygenIsotopicComposition(summary["d18Oc (VSMOW)"], summary["Temperature"]), 2)

# Concatenate groups back together
dfFinal = pd.concat([group for group in vdf])

# Write to Excel file
#with pd.ExcelWriter("Report.xlsx") as writer:
#    dfFinal.to_excel(writer, sheet_name="Evaluation", index=False)
#    summary.to_excel(writer, sheet_name="Summary", index=False)
