using FITSFiles
using Plots

cm2in(x) = 0.3937008x
cm2px(x) = Int(round(cm2in(x) * 100))

pyplot()

default(titlefont = (12, "serif"), 
    guidefont = (10, "serif"), 
    legendfont = (8, "serif"), 
    tickfont = (8, "serif"), 
    gridalpha=0.,
    minorticks=true,
    dpi=300,
    size=cm2px.((12,9))
)

file = fits("/home/brad/Downloads/000125000_rsl_reproc/rsl000125000_clbr_rise_Hp_spectra.fits")

display(plot(file[3].data.PI, file[3].data.DIFF; xlabel="PI Channel", ylabel="Difference (counts)", title="Rise-Time Screening", label=nothing))

file = fits("/home/brad/Downloads/000125000_rsl_reproc/rsl000125000_clbr_prox_Hp_spectra.fits")

display(plot(file[3].data.PI, file[3].data.DIFF; xlabel="PI Channel", ylabel="Difference (counts)", title="Proximity Screening", label=nothing))

file = fits("/home/brad/Downloads/analysis/xa000125000xtd_p031100010_cl_xpc_clnevt.fits")