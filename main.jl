"""
This project allows for data preparation and analysis in the MuSiC (Multi-Sample Single-Cell) pipeline, which has been refactored into Julia.

Author: Amit Halkhoree
Created: Jan 2 2023
Adapted from: LuoHan https://github.com/Prayforhanluo/MuSiC_Python/blob/master/MuSiC.py
"""

include("src/data_loading.jl")
include("src/data_preparation.jl")
include("src/data_analysis.jl")

using .DataLoading
using .DataPreparation
using .DataAnalysis
using DelimitedFiles
using DataFrames
using LinearAlgebra
using Statistics
using NNLS

# Load data
sc_Meta = get_SC_Meta_Info(ARGS[1])
sc_Count = get_SC_Count_Info(ARGS[2])
bulk_Meta = get_bulk_Meta_Info(ARGS[3])
bulk_Count = get_bulk_Count_Info(ARGS[4])

# Define select_ct
select_ct = ["alpha", "beta", "delta", "gamma", "acinar", "ductal"]

# Weight gene by variance of cross-subject
Results_var = music_prop(bulk_Meta, bulk_Count, sc_Meta, sc_Count, ct_cov=false, select_ct=select_ct)

# Weight gene by co-variance of cross-subject
Results_cov = music_prop(bulk_Meta, bulk_Count, sc_Meta, sc_Count, ct_cov=true, select_ct=select_ct)

# Call main function with command-line arguments
main(ARGS)