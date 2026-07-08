using FITSFiles
using Plots

function GetCoordCheckCommand(file, inpt, instr="RESOLVE")
    RA_NOM = file[1].cards["RA_NOM"]
    DEC_NOM = file[1].cards["DEC_NOM"]
    PA_NOM = file[1].cards["PA_NOM"]

    println("coordpnt input=$inpt outfile=NONE telescop=XRISM instrume=$instr teldeffile=CALDB startsys=DET stopsys=RADEC ra=$RA_NOM dec=$DEC_NOM roll=$PA_NOM ranom=$RA_NOM decnom=$DEC_NOM clobber=yes")
end

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

file = fits("/home/brad/Downloads/analysis/xa000125000rsl_p0px1000_cl2.evt")

GetCoordCheckCommand(file, "3.5,3.5", "RESOLVE")
GetCoordCheckCommand(file, "733,730", "XTEND")

file = fits("/home/brad/Downloads/analysis/xa000125000xtd_src.pi")[2].cards["BACKSCAL"]

file = fits("data/nu80402315002A01_sr_1000.pha")
file = fits("data/xa000125000xtd_src.pi")

file[2].cards["RESPFILE"] = "xa000125000xtd_p031100010_src.rmf"