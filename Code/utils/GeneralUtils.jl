using WAV

# Append to a named tuple
setindex(nt::NamedTuple, key::AbstractString, value) = merge(nt, (Symbol(key) => value,))

# Play a sound upon completion
function CompleteSound()
    """A function to play a sound when the program finishes running"""
    y, fs = wavread("./complete.wav")
    wavplay(y, fs)
end