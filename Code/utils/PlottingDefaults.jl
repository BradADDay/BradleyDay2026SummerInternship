using Plots

pyplot()

cm2in(x) = 0.3937008x
cm2px(x) = Int(round(cm2in(x) * 100))

default(
    # Fonts
    titlefont = (12, "serif"), 
    guidefont = (10, "serif"), 
    legendfont = (8, "serif"), 
    tickfont = (8, "serif"), 
    colorbar_titlefont = (10, "serif"),
    # Grids
    gridalpha=0.5,
    minorticks=true,
    minorgridalpha=0.2,
    minorgrid=true,
    # Image
    dpi=1200,
    size=cm2px.((12,9)),
    framestyle=:box
)