using DataFrames, CSV, ProgressBars

include("utils/Deformations.jl")

as = range(0, 0.998, 10)

inpath = "/home/brad/Documents/SummerInternship/tabledataFinal5Merged/"
outpath = "/home/brad/Documents/SummerInternship/tabledataFinal5MergedZeros/"

a = 0.998

for a in as[2:end]
    j=1
    df = DataFrame(CSV.File(joinpath(inpath, "$a.csv")))

    outdf = DataFrame()
    for column in ProgressBar(names(df))
        if column[end-1:end] !== "_1"
            colarray = [parse(Float64, ss) for ss in split(column, ", ")]
            if IsValid(colarray[2], colarray[1], a)
                insertcols!(outdf, j, column => df[!, column])
                j+=1
            end
        end
    end

    CSV.write(joinpath(outpath, "$a.csv"), outdf)
    
end