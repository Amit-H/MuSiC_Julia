# MuSiC implementation in Julia - Under Construction

This project allows for data preparation and analysis in the MuSiC (Multi-Sample Single-Cell) pipeline, which has been refactored into Julia for speed. This project is not complete yet and remains under construction.

## Inspiration
This code was based on the following:
LuoHan's Python implementation of MuSiC: https://github.com/Prayforhanluo/MuSiC_Python/blob/master/MuSiC.py

Bulk tissue cell type deconvolution with multi-subject single-cell expression reference: https://www.nature.com/articles/s41467-018-08023-x#Sec8

## Raison d'etre

The original MuSiC implementation is done in R, and I wanted to see if implementing it in Julia yieleded performance benefits when used on larger datasets.

## Project Structure

- `main.jl`: The main entry point of the application.
- `src/`: This directory contains the main modules of the application.
  - `data_analysis.jl`: Contains functions for data analysis.
  - `data_loading.jl`: Contains functions for loading data.
  - `data_preparation.jl`: Contains functions for preparing data.
- `Project.toml` and `Manifest.toml`: These files are used by Julia's package manager.

## How to Run

To run this project, you need to pass 4 arguments to the `main.jl` script. These arguments are the paths to the files that contain the SC Meta Info, SC Count Info, bulk Meta Info, and bulk Count Info, respectively.

```bash
julia main.jl path_to_sc_meta_info path_to_sc_count_info path_to_bulk_meta_info path_to_bulk_count_info
```

## Dependencies
This project uses the following Julia packages:

- DelimitedFiles
- DataFrames
- LinearAlgebra
- Statistics

The environment can be replicated by calling `instantiate` from Julia's environment manager.
