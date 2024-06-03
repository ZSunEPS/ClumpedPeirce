using CSV, DataFrames, Statistics, XLSX

ptable = CSV.read("Peirces Table.csv", DataFrame)

raw_export = CSV.read("data.csv", DataFrame)

df = DataFrame(
    Sample = raw_export[:, "Easotope Name"], d13C = raw_export[:, "d13C VPDB (Final)"], d18O = raw_export[:, "d18O VPDB (Final)"], D47 = raw_export[:, "D47 CDES (Final)"], 
    d13C_Outlier = 0, d18O_Outlier = 0, D47_Outlier = 0, Total_Outlier = 0
)
dropmissing!(df)

gdf = groupby(df, "Sample")

L = length(gdf)

vdf = [DataFrame(gdf[i]) for i in 1:L]

summary = DataFrame(Sample = [vdf[i][1, "Sample"] for i in 1 : L], n = 0, d13C = 0.0, d13C_SE = 0.0, d18O = 0.0, d18O_SE = 0.0, D47 = 0.0, D47_SE = 0.0)

function PierceOutlierDetect(Obs::Vector{Float64}, Outl::Vector{Int64}, d::Int64)
    n = length(Obs)
    Obs_mean = round(mean(Obs), digits = d)
    Obs_SD = round(std(Obs), digits = d + 1)
    Obs_Deviation = abs.(Obs .- Obs_mean)
    global LastOutlierDetected = -1
    global NewOutlierDetected = 0
    global Step = 0
    while NewOutlierDetected > LastOutlierDetected
        global LastOutlierDetected = NewOutlierDetected
        global Step = Step + 1
        R = ptable[n - 2, Step]
        Obs_MaxDeviation  = round(Obs_SD * R, digits = d + 1)
        for i in 1:n
            if Obs_Deviation[i] > Obs_MaxDeviation
                Outl[i] = 1
            end
        end
        global NewOutlierDetected = sum(Outl)
    end
end

for k in 1 : L
    PierceOutlierDetect(vdf[k][!, "d13C"], vdf[k][!, "d13C_Outlier"], 2)
    PierceOutlierDetect(vdf[k][!, "d18O"], vdf[k][!, "d18O_Outlier"], 2)
    PierceOutlierDetect(vdf[k][!, "D47"], vdf[k][!, "D47_Outlier"], 3)
    vdf[k][!, "Total_Outlier"] = (vdf[k][!, "d13C_Outlier"] .+ vdf[k][!, "d18O_Outlier"] .+ vdf[k][!, "D47_Outlier"])
    temp_v = filter(row -> row["Total_Outlier"] < 1, vdf[k])
    m = nrow(temp_v)
    summary[k, "n"] = m
    summary[k, "d13C"] = round(mean(temp_v[!, "d13C"]), digits = 2)
    summary[k, "d18O"] = round(mean(temp_v[!, "d18O"]), digits = 2)
    summary[k, "D47"] = round(mean(temp_v[!, "D47"]), digits = 3)
    summary[k, "d13C_SE"] = round(std(temp_v[!, "d13C"]) / sqrt(m), digits = 2)
    summary[k, "d18O_SE"] = round(std(temp_v[!, "d18O"]) / sqrt(m), digits = 2)
    summary[k, "D47_SE"] = round(std(temp_v[!, "D47"]) / sqrt(m), digits = 3)
end

summary[:, "Temperature"] = sqrt.((0.0391 * 10^6) ./ (summary[!, "D47"] .- 0.154)) .- 273.15
summary[:, "Temperature_SE"] = abs.(sqrt(0.0391) .* 10^3 .* (-1/2) .* (summary[!, "D47"] .- 0.154).^(-3/2) .* summary[!, "D47_SE"])
summary[:, "Temperature"] = Int64.(round.(summary[:, "Temperature"], digits = 0))
summary[:, "Temperature_SE"] = Int64.(round.(summary[:, "Temperature_SE"], digits = 0))

summary[:, "d18Oc (VSMOW)"] = round.(summary[!, "d18O"] .* 1.03091 .+ 30.91, digits = 2)
summary[:, "d18Ow (VSMOW)"] = (1000 .+ summary[:, "d18Oc (VSMOW)"]) ./ (exp.((18030 ./ (summary[:, "Temperature"] .+ 273.15) .- 32.42) ./ 1000)) .- 1000
summary[:, "d18Ow (VSMOW)"] = round.(summary[:, "d18Ow (VSMOW)"], digits = 2)

df = copy(vdf[1])
for i in 2 : L
    global df = vcat(df, vdf[i])
end

XLSX.writetable("Report.xlsx", "Evaluation" => df, "Summary" => summary)