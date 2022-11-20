module DataLoading

using LazyModules

import CSV # lazy import causes world age issues
using Statistics
using DataFrames

include("data_loading/data.jl")
include("data_loading/prior_segmentation.jl")
include("data_loading/cli_wrappers.jl")

end