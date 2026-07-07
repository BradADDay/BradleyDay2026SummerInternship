f = open("/home/brad/Downloads/000125000_rsl_reproc/rslbratios.log", "r")

lines = readlines(f)

keeplines = []

for i in eachindex(lines)
    try
        if lines[i][1:20] == "rslbratios: WARNING:"
            append!(keeplines, i)
        end
    catch
    end
end

lines = lines[keeplines]