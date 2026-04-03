"""
    struct PackedArray{T, W <: Unsigned} words::Vector{W} bitmap::BitVector length::Int end -> Return type

A bitpacked array of type T, using words of size W. Keeps a bitmap to remember boundaries.
"""
struct PackedArray{T, W <: Unsigned}
    words::Vector{W}
    bitmap::BitVector
    length::Int
end

bitsizeof(e) = sizeof(e) * 8
PackedArray{T, W}() where {T, W<:Unsigned} = PackedArray{T, W}(W[0], BitVector(falses(bitsizeof(W))), 0)
PackedArray{T}() where {T} = PackedArray{T, UInt64}(UInt64[0], BitVector(falses(64)), 0)

# Computes the packed version of elem on word, returns the pack size by grabbing left most zero
function prepack(elem::T, ::Type{W}=UInt64) where {T, W<:Unsigned}
    !isbits(elem) && throw("Type $T is not bits type and cannot be packed")
    
    # Needs to reinterpret as a UInt, then grab the left most zero to compute pack size.
    # pack = Core.Intrinsics.bitcast(W, elem) 
    pack = W(elem & typemax(W)) # /!\ Will round down to typemax, careful
    pack == 0 && return pack, 1
    return pack, bitsizeof(W) - leading_zeros(pack)
end

# assumes bitmap is always word_size*words long
word_bitmap_slice(word_id::Integer, word_size::DataType=UInt64) = (word_id-1)*bitsizeof(word_size)+1:(word_id-1)*bitsizeof(word_size)+bitsizeof(word_size)

# Last bit position that's set in a word 
function word_last_set(arr::PackedArray{T, W}, word_id::Integer) where {T, W<:Unsigned}
    last_set = findlast(@view arr.bitmap[word_bitmap_slice(word_id, W)])
    return !isnothing(last_set) ? last_set : 0
end

# Version that recomputes a last_set index
function update_word_bitmap!(arr::PackedArray{T, W}, word_id::Integer, bitsize::Integer) where {T, W<:Unsigned}
    arr.bitmap[word_bitmap_slice(word_id, W)[1]+word_last_set(arr, word_id)+bitsize-1] = true
end

# Version when a precomputed last_set index is already known (saves time)
function update_word_bitmap!(arr::PackedArray{T, W}, word_id::Integer, bitsize::Integer, last_set::Integer) where {T, W<:Unsigned}
    arr.bitmap[word_bitmap_slice(word_id, W)[1]+last_set+bitsize-1] = true
end

# function Base.push!(arr::PackedArray{T, W}, word::W) where {T, W<:Unsigned}
#     push!(arr.words, word)
#     append!(arr.bitmap, zeros(bitsizeof(W)))
#     return length(arr.words)
# end

# Adds elem to word at position word_id of array. Adds a new word if not enough space. Returns id of word pushed to
function Base.push!(arr::PackedArray{T1, W}, elem::T2, word_id::Integer) where {T1, T2, W<:Unsigned}
    # Needs to reinterpret elem and compute packed size
    elem, bitsize = prepack(elem, W)

    last_set = word_last_set(arr, word_id)
    # Needs to check if the word has enough space (check boundaries)
    if bitsizeof(W)-last_set >= bitsize
        # println(elem<<(bitsizeof(W)-last_set-bitsize))
        arr.words[word_id] |= (elem<<(bitsizeof(W)-last_set-bitsize))
        update_word_bitmap!(arr, word_id, bitsize, last_set)
        return word_id
    else
        # If it doesn't, push a word, add it to it, update boundaries and return new word_id
        push!(arr.words, W(0))
        append!(arr.bitmap, falses(bitsizeof(W)))
        return push!(arr, elem, length(arr.words))
    end
end

# Push elem to a new word of array. Returns id of word pushed to.
function Base.push!(arr::PackedArray{T, W}, elem) where {T, W<:Unsigned}
    # Without a word_id, we always add a new word with that element (mainly used at initialization)
    elem, _ = prepack(elem, W)
    push!(arr.words, W(0))
    append!(arr.bitmap, falses(bitsizeof(W)))
    return push!(arr, elem, length(arr.words))
end

# Slices bitmap to retrieve value of elem elem_id of word word_id
function get_val(arr::PackedArray{T, W}, word_id::Integer, elem_id::Integer) where {T, W<:Unsigned}
    word = arr.words[word_id]
    bitmap_word = arr.bitmap[word_bitmap_slice(word_id, W)]
    sum(bitmap_word) < elem_id && throw("Out of Bound")
    start = 0; finish = 0; counter = 0; ones = 0
    for bit in bitmap_word
        counter += 1
        ones == elem_id && break
        (start, finish, ones) = bit ? (finish, counter, ones+1) : (start, finish, ones)
    end
    return word << start >> (start + (bitsizeof(W)-finish))
end

# Pushes new empty word to array. Return its id
function new_word!(arr::PackedArray{T, W}) where {T, W<:Unsigned}
    push!(arr.words, W(0))
    append!(arr.bitmap, falses(bitsizeof(W)))
    return length(arr.words)
end

# Assemble the full meaning of word word_id (retrieves all values on word by slicing bitmap)
function assemble_word(arr::PackedArray{T, W}, word_id::Integer) where {T, W<:Unsigned}
    bitmap_word = arr.bitmap[word_bitmap_slice(word_id, W)]
    elem_id = 1
    assembly = W[]
    while sum(bitmap_word) >= elem_id
        push!(assembly, get_val(arr, word_id, elem_id))
        elem_id += 1
    end
    return assembly
end

# Deduplicates words, get permutation array of words, apply to bitmap to update it aswell, returns deduped array and permutations.
# /!\ Might make arr wrong, need to investigate why? Deduped array stays correct, low priority
function permdedup(arr::PackedArray{T, W}) where {T, W<:Unsigned}
    uniq = Vector{W}()              
    indexof = Dict{Pair{W, SubArray}, UInt32}()         
    perm = Vector{UInt32}()
    global_perm = Vector{UInt32}(undef, length(arr.words))

    @inbounds @showprogress "Deduping Count Words..." for i in eachindex(arr.words)
        word = arr.words[i]
        bits = @view arr.bitmap[word_bitmap_slice(i, W)]
        # assembled_word = assemble_word(arr, i)  ## Way too slow (and makes dict scale like shit.... OBVIOUSLY)
        j = get(indexof, word=>bits, 0)
        if j == 0
            push!(uniq, word)
            j = length(uniq)
            indexof[word=>bits] = UInt32(j)
            push!(perm, UInt32(i))
        end
        global_perm[i] = j
    end

    bitmap_slices = word_bitmap_slice.(perm, W)
    deduped_arr = PackedArray{T, W}(uniq, arr.bitmap[vcat(bitmap_slices...)], length(uniq))
    return deduped_arr, perm, global_perm
end