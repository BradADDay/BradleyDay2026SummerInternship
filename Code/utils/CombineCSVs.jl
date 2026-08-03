using DataFrames, CSV, ProgressBars

include("Deformations.jl")

as = range(0, 0.998, 10)

inpath1 = "/home/brad/Documents/SummerInternship/tabledataFinal5Merged/"
inpath2 = "/home/brad/Documents/SummerInternship/tabledataExtendedH/"
outpath = "/home/brad/Documents/SummerInternship/finaltabledata/"

for a in as
    
    try 
        df1 = DataFrame(CSV.File("$inpath1$a.csv"))
        df2 = DataFrame(CSV.File("$inpath2$a.csv"))
        df3 = hcat(df1, df2)
        println(join([size(df1), size(df2), size(df3)], " "))
    catch
    end
end