using Pkg
Pkg.activate(".")

using BenchmarkTools
using SatelliteToolboxGravityModels
using JLD2

@benchmark SatelliteToolboxGravityModels.parse_icgem("/home/dani/.julia/scratchspaces/b3500434-02a6-41be-b4da-749efa17604c/icgem/EGM2008.gfc")

# BenchmarkTools.Trial: 1 sample with 1 evaluation per sample.
#  Single result which took 5.862 s (10.28% GC) to evaluate,
#  with a memory estimate of 1.67 GiB, over 42290806 allocations.

# =================== new ====================
# BenchmarkTools.Trial: 1 sample with 1 evaluation per sample.
#  Single result which took 5.103 s (4.35% GC) to evaluate,
#  with a memory estimate of 1.92 GiB, over 39889494 allocations.

stgm_object = SatelliteToolboxGravityModels.parse_icgem("/home/dani/.julia/scratchspaces/b3500434-02a6-41be-b4da-749efa17604c/icgem/EGM2008.gfc")

N, M = 15, 15
P = zeros((16, 16))
dP = zeros((16, 16))
r = [6378137.0 + 400000.0, 1234.0, -1234.0]

@benchmark SatelliteToolboxGravityModels.GravityModels.gravitational_acceleration(stgm_object, r; max_degree=N, max_order=M, P=P, dP=dP)
# BenchmarkTools.Trial: 10000 samples with 10 evaluations per sample.
#  Range (min … max):  1.833 μs …  3.303 μs  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     1.883 μs              ┊ GC (median):    0.00%
#  Time  (mean ± σ):   1.897 μs ± 75.325 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

#     ▂▄▅█▅▃                                                    
#   ▂▆███████▆▄▃▂▂▂▂▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁ ▂
#   1.83 μs        Histogram: frequency by time        2.33 μs <

#  Memory estimate: 128 bytes, allocs estimate: 3.

# BenchmarkTools.Trial: 10000 samples with 9 evaluations per sample.
#  Range (min … max):  2.222 μs …   6.660 μs  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     2.260 μs               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   2.286 μs ± 137.221 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

#   ▃▇█▇▇▄▃▂▁                                                   ▂
#   ██████████▇▆▆▅▅▅▆▇▇▆▆▅▆▄▄▃▁▄▆▅▅▆▅▆▆▅▅▆▅▆███▇▇▇▅▆▃▄▁▁▁▃▅▁▁▃▄ █
#   2.22 μs      Histogram: log(frequency) by time      2.98 μs <

#  Memory estimate: 128 bytes, allocs estimate: 3.

@benchmark jldsave("foo.jld2"; gm=stgm_object)
# BenchmarkTools.Trial: 1 sample with 1 evaluation per sample.
#  Single result which took 13.246 s (1.65% GC) to evaluate,
#  with a memory estimate of 1.61 GiB, over 48019594 allocations.

# BenchmarkTools.Trial: 52 samples with 1 evaluation per sample.
#  Range (min … max):  93.325 ms … 107.131 ms  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     96.181 ms               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   97.998 ms ±   4.154 ms  ┊ GC (mean ± σ):  0.00% ± 0.00%

#     ▂ █▂    ▂                                   ▅               
#   ▅██▅██▅█▅██▅▅▅▅▁▅▅▁▁▅▅▅▁▅▁▁▁▁▁█▁▁▁▁▁█▁▁▅▁▁▅█▁▁█▁▅▁▁▁▁▁▁▅▅▁▁█ ▁
#   93.3 ms         Histogram: frequency by time          106 ms <

#  Memory estimate: 56.19 KiB, allocs estimate: 1164.


@benchmark gm = load("foo.jld2")
# BenchmarkTools.Trial: 1 sample with 1 evaluation per sample.
#  Single result which took 9.320 s (10.00% GC) to evaluate,
#  with a memory estimate of 3.04 GiB, over 76808926 allocations.

# julia> @benchmark gm = load("foo.jld2")
# BenchmarkTools.Trial: 109 samples with 1 evaluation per sample.
#  Range (min … max):  38.781 ms … 81.458 ms  ┊ GC (min … max):  0.00% … 50.32%
#  Time  (median):     45.394 ms              ┊ GC (median):    12.84%
#  Time  (mean ± σ):   45.895 ms ±  5.091 ms  ┊ GC (mean ± σ):  13.67% ±  6.64%

#            ▂█▄                                                 
#   ▃▃▂▁▃▂▁▁▄███▄▃▁▁▁▁▁▁▁▃▁▁▁▁▁▁▁▃▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▂ ▂
#   38.8 ms         Histogram: frequency by time        75.3 ms <

#  Memory estimate: 73.30 MiB, allocs estimate: 895.