using Random

kmers = psort!(rand(UInt64, 100_000_000))
darr = DeltaArray(kmers)