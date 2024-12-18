using CSV, DataFrames, Statistics, XLSX

peircesTable = CSV.read("Peirces Table.csv", DataFrame) # Load Peirce's Criterion Table

rawExport = CSV.read("Data.csv", DataFrame) # Read data exported from Easotope. Replace "data.csv" with the path of your .csv file.

df = DataFrame(
    Sample = rawExport[:, "Easotope Name"], 
    d13C = rawExport[:, "d13C VPDB (Final)"], 
    d18O = rawExport[:, "d18O VPDB (Final)"], 
    D47 = rawExport[:, "D47 CDES (Final)"], 
    d13C_Outlier = 0, 
    d18O_Outlier = 0, 
    D47_Outlier = 0, 
    Total_Outlier = 0
) # Create a DataFrame with desired columns and initialize outlier indicators
dropmissing!(df) # Drop rows withing Missing data.

gdf = groupby(df, "Sample") # Group data by "Sample".
vdf = collect(gdf) # A vector of subdataframes.

L = length(vdf) # Count samples.
summary = DataFrame(
    Sample = [g[1, "Sample"] for g in vdf], 
    n = 0, 
    d13C = 0.0, 
    d13C_SE = 0.0, 
    d18O = 0.0, 
    d18O_SE = 0.0, 
    D47 = 0.0, 
    D47_SE = 0.0
) # Create a summary for evaluated data.

# Function to perform Peirce's outlier detection
function PeirceOutlierDetect(Obs, Outl, Table::DataFrame, d::Int64)
    n = length(Obs)
    if n < 2  # Not enough data for analysis
        return Outl .= 0
    end
    
    Obs_mean = round(mean(Obs), digits = d)
    Obs_SD = round(std(Obs), digits = d + 1)
    Obs_Deviation = abs.(Obs .- Obs_mean)
     LastOutlierDetected = -1
     NewOutlierDetected = 0
     Step = 0
    while NewOutlierDetected > LastOutlierDetected
        LastOutlierDetected = NewOutlierDetected
        Step += 1
        R = Table[n - 2, Step]
        Obs_MaxDeviation  = round(Obs_SD * R, digits = d + 1)
        Outl .= Obs_Deviation .> Obs_MaxDeviation # Mark outliers
        NewOutlierDetected = sum(Outl)
    end
end

for k in 1 : L # Process data per group (or sample)
    PeirceOutlierDetect(vdf[k][!, "d13C"], vdf[k][!, "d13C_Outlier"], peircesTable, 2)
    PeirceOutlierDetect(vdf[k][!, "d18O"], vdf[k][!, "d18O_Outlier"], peircesTable, 2)
    PeirceOutlierDetect(vdf[k][!, "D47"], vdf[k][!, "D47_Outlier"], peircesTable, 3)

    # Total outlier flag
    vdf[k][!, "Total_Outlier"] = (vdf[k][!, "d13C_Outlier"] .+ vdf[k][!, "d18O_Outlier"] .+ vdf[k][!, "D47_Outlier"])
    
    # Filter non-outliers and calculate statistics
    cleanData = filter(row -> row["Total_Outlier"] == 0, vdf[k])

    m = nrow(cleanData)
    summary[k, "n"] = m
    summary[k, "d13C"] = round(mean(cleanData[!, "d13C"]), digits = 2)
    summary[k, "d18O"] = round(mean(cleanData[!, "d18O"]), digits = 2)
    summary[k, "D47"] = round(mean(cleanData[!, "D47"]), digits = 3)
    summary[k, "d13C_SE"] = round(std(cleanData[!, "d13C"]) / sqrt(m), digits = 2)
    summary[k, "d18O_SE"] = round(std(cleanData[!, "d18O"]) / sqrt(m), digits = 2)
    summary[k, "D47_SE"] = round(std(cleanData[!, "D47"]) / sqrt(m), digits = 3)
end

# The Δ47-Temperature equation is from Anderson et al. (2021).
summary[:, "Temperature"] = sqrt.((0.0391 * 10^6) ./ (summary[!, "D47"] .- 0.154)) .- 273.15
summary[:, "Temperature_SE"] = abs.(sqrt(0.0391) .* 10^3 .* (-1/2) .* (summary[!, "D47"] .- 0.154).^(-3/2) .* summary[!, "D47_SE"])
summary[:, "Temperature"] = Int64.(round.(summary[:, "Temperature"], digits = 0))
summary[:, "Temperature_SE"] = Int64.(round.(summary[:, "Temperature_SE"], digits = 0))

summary[:, "d18Oc (VSMOW)"] = round.(summary[!, "d18O"] .* 1.03091 .+ 30.91, digits = 2)

# The δ18O fractionation equation is from kim and O'Neil (1997).
summary[:, "d18Ow (VSMOW)"] = (1000 .+ summary[:, "d18Oc (VSMOW)"]) ./ (exp.((18030 ./ (summary[:, "Temperature"] .+ 273.15) .- 32.42) ./ 1000)) .- 1000
summary[:, "d18Ow (VSMOW)"] = round.(summary[:, "d18Ow (VSMOW)"], digits = 2)

dfFinal = reduce(vcat, vdf) # Concatenate groups back together

#XLSX.writetable("Report.xlsx", "Evaluation" => dfFinal, "Summary" => summary)
