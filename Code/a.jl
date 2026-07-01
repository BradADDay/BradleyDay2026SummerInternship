using SpectralFitting, XSPECModels, FITSFiles

DATAFILE = "/home/brad/Documents/SummerInternship/Code/model.FITS"

@xspecmodel fits(DATAFILE) struct XS_JohannsenLampPost{T} <: AbstractSpectralModel{T, Additive}
    "Normalisation."
    K::T
    "Photon index."
    a::T
end

function XS_JohannsenLampPost(; K = FitParam(1.0), a = FitParam(1.0))
    XS_JohannsenLampPost{typeof(K)}(K, a)
end

model = XS_JohannsenLampPost()