module NeoKCT
include("KCTLayers.jl")
for n in names(@__MODULE__; all=true)
    if Base.isidentifier(n) && !(n in (Symbol("NeoKCT"), :eval, :include))
        @eval export $n
    end
end
end
