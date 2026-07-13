using Plots

pyplot()

cm2in(x) = 0.3937008x
cm2px(x) = Int(round(cm2in(x) * 100))

default(titlefont = (12, "serif"), 
    guidefont = (10, "serif"), 
    legendfont = (8, "serif"), 
    tickfont = (8, "serif"), 
    colorbar_titlefont = (10, "serif"),
    gridalpha=0.5,
    minorticks=true,
    dpi=300,
    size=cm2px.((12,9)),
    framestyle=:box,
    minorgridalpha=0.2,
    minorgrid=true
)