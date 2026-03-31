using Kmers, Kmers.BioSequences

struct AAAlphabet <: Alphabet end

BioSequences.BitsPerSymbol(::AAAlphabet) = BioSequences.BitsPerSymbol{5}()
Base.eltype(::Type{AAAlphabet}) = AminoAcid
bin(x) = Base.bin(x, 64, false)

function BioSequences.encode(::AAAlphabet, aa::AminoAcid)
    reinterpret(UInt8, aa) > reinterpret(UInt8, AA_Gap) && return nothing
    return convert(UInt, reinterpret(UInt8, aa))
end

@inline function BioSequences.decode(::AAAlphabet, x::UInt)
    return reinterpret(AminoAcid, x % UInt8)
end

offset(::Km) where {Km<:Kmer{A, K, N}} where {A<:Alphabet, K, N} = N*64 - BioSequences.bits_per_symbol(A())*K
true_index(k::Km, i::Int) where {Km<:Kmer{A, K, N}} where {A<:Alphabet, K, N} = offset(k)÷BioSequences.bits_per_symbol(A()) + i
carry(c::Int) = c-c÷5
# chunk(k::Km, i::Int) where {Km<:Kmer} = true_index(k, i)÷(12+carry(true_index(k, i)÷12))

function chunk(k::Km, i::Int) where {Km<:Kmer}
    off_idx = true_index(k, i)
    interval = 1:12
    bps = BioSequences.bits_per_symbol(Alphabet(k)) % UInt
    cycle = length(k.data) % bps
    cur_chunk = 0

    while !(off_idx in interval)
        interval = interval[end]+1:interval[end] + 12 + (cycle != 0)
        # println(interval)
        cur_chunk += 1
        cycle = (length(k.data) - cur_chunk) % bps
    end
    # println("\n$interval")
    return cur_chunk
end

@inline function BioSequences.extract_encoded_element(seq::Kmer{<:AAAlphabet}, i::Integer)
    off = offset(seq)
    bps = BioSequences.bits_per_symbol(Alphabet(seq)) % UInt
    start_bit, end_bit = off+bps*(i-1), off+bps*i
    start_chunk, end_chunk = ceil(start_bit/64), ceil(end_bit/64)
    trailing = start_bit % 64
    # println(trailing)
    # println(start_chunk, end_chunk)

    if trailing == 0  # Just looped over (start chunk is off in these cases, there must be a more general implementation)
        return seq.data[start_chunk+1] >> (64-bps)

    elseif start_chunk == end_chunk  # All in one chunk (std case)
        return seq.data[start_chunk] << trailing >> (64-bps)

    else
        overspill = end_bit % 64 # split across 2 chunks
       return (seq.data[start_chunk] << trailing >> (64-bps)) | (seq.data[end_chunk] >> (64-overspill))

    end
end