using WAV

# Append to a named tuple
setindex(nt::NamedTuple, key::AbstractString, value) = merge(nt, (Symbol(key) => value,))

const COMPLETESOUNDPATH = "Code/utils/complete.wav"

# Play a sound upon completion
function CompleteSound(; path=COMPLETESOUNDPATH)
    """A function to play a sound when the program finishes running"""
    y, fs = wavread(path)
    wavplay(y, fs)
end