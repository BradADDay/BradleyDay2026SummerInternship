using SpectralFitting


DATADIR = "/home/brad/Documents/SummerInternship/data/"

specPath = joinpath(DATADIR, "nu80502304006B01_sr_1000.pha")

data = OGIPDataset(specPath)

plot(data, xlims = (5, 8))