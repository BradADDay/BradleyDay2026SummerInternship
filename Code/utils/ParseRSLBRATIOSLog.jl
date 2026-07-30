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

function ParseRSBLRATIOS(f; limit=50)

    # Reading the file
    lines = readlines(f)

    filtLines = []

    # Filtering for useful lines
    for line in lines
        try
            if line[1:27] == "rslbratios: WARNING: (HIGH)"
                append!(filtLines, [line[29:end]])
            end
        catch
        end
    end

    # Constructing a dataframe from the data
    df = DataFrame((Quadrant=[], Counts=[]))

    for line in filtLines[1:2:end]
        push!(df, [FindAfter(line, "quadrant "), FindAfter(line, "with ")])
    end

    # Filtering the dataframe for invalid values
    df = subset(df, :Counts => a -> a .> limit)

    # Returning the invalid quadrants
    unique(df.Quadrant), df
end