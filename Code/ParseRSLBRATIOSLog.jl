using DataFrames

function FindAfter(string, substring)

    out = ""

    # Index of the first character after the substring
    idx = findfirst(substring, string)[end] + 1

    for i in string[idx:end]
        if i == string[idx-1]
            break
        else
            out = "$out$i"
        end
    end

    parse(Int64, out)
end

f = open("/home/brad/Downloads/000125000_rsl_reproc/rslbratios.log", "r")

lines = readlines(f)

filtLines = []

for line in lines
    try
        if line[1:27] == "rslbratios: WARNING: (HIGH)"
            append!(filtLines, [line[29:end]])
        end
    catch
    end
end

df = DataFrame((Quadrant=[], Counts=[]))

for line in filtLines[1:2:end]
    push!(df, [FindAfter(line, "quadrant "), FindAfter(line, "with ")])
end

df = subset(df, :Counts => a -> a .> 50)

unique(df.Quadrant)