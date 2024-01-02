struct BulkMatrix
    """
    Calculate relative abundance
    """
    X::Matrix{Float64}
    by_col::Bool

    function BulkMatrix(X, by_col=true)
        new(X, by_col)
    end
end

function DataMatrixConstruction(bulk_matrix::BulkMatrix)
    """
    Abundance computation.
    """
    if any(x -> x < 0, bulk_matrix.X)
        throw(Exception("Negative entry appears"))
    end

    bulk_matrix.X ./= sum(bulk_matrix.X, dims=bulk_matrix.by_col ? 1 : 2)

    return bulk_matrix
end

struct SingleCellMatrix
    Meta::DataFrame
    Count::DataFrame
    nonzero::Bool
    markers::Vector{String}
    select_ct::Vector{String}
    ct_cov::Bool

    function SingleCellMatrix(Meta, Count, nonzero=true, markers=nothing,
                              select_ct=nothing, ct_cov=false)
        new(Meta, Count, nonzero, markers, select_ct, ct_cov)
    end
end

function DataConstruction!(single_cell_matrix::SingleCellMatrix)
    if single_cell_matrix.select_ct !== nothing
        mask = in.(single_cell_matrix.Meta[:, "cellType"], single_cell_matrix.select_ct)
        single_cell_matrix.Meta = single_cell_matrix.Meta[mask, :]
        single_cell_matrix.Count = single_cell_matrix.Count[:, single_cell_matrix.Meta[:, "sampleName"]]
    end

    if single_cell_matrix.nonzero
        nz_gene = sum(eachrow(single_cell_matrix.Count) .!= 0) .!= 0
        single_cell_matrix.Count = single_cell_matrix.Count[nz_gene, :]
    end

    celltypes = single_cell_matrix.select_ct
    samples = Set(single_cell_matrix.Meta[:, "sampleID"])

    println("Creating Relative Abundance Matrix ...")
    single_cell_matrix.theta = Theta(single_cell_matrix, celltypes, samples)

    println("Creating Variance Matrix ...")
    if single_cell_matrix.ct_cov
        single_cell_matrix.sigma_ct = Sigma(single_cell_matrix, celltypes, samples)
    else
        single_cell_matrix.sigma = Sigma(single_cell_matrix, celltypes, samples)
    end

    println("Creating Library Size Matrix ...")
    single_cell_matrix.LS = LibrarySize(single_cell_matrix, celltypes, samples)
    single_cell_matrix.MLS = mean(single_cell_matrix.LS, dims=1)

    println("Creating Design Matrix ...")
    single_cell_matrix.DM = zeros(size(single_cell_matrix.theta))
    theta = Matrix(single_cell_matrix.theta)
    for i in 1:length(celltypes)
        single_cell_matrix.DM[:, i] .= theta[:, i] .* single_cell_matrix.MLS[i]
    end

    println("Markers genes selecting ...")
    if single_cell_matrix.markers === nothing
        Genes = names(single_cell_matrix.Count)
        if !single_cell_matrix.ct_cov
            single_cell_matrix.sigma = DataFrame(single_cell_matrix.sigma, Genes, celltypes)
        else
            single_cell_matrix.sigma_ct = DataFrame(single_cell_matrix.sigma_ct, celltypes)
        end
        single_cell_matrix.DM = DataFrame(single_cell_matrix.DM, Genes, celltypes)
        single_cell_matrix.theta = DataFrame(single_cell_matrix.theta, Genes, celltypes)
    else
        mask = in.(names(single_cell_matrix.Count), single_cell_matrix.markers)
        single_cell_matrix.Count = single_cell_matrix.Count[mask, :]
        Genes = names(single_cell_matrix.Count)
        single_cell_matrix.DM = single_cell_matrix.DM[mask, :]
        single_cell_matrix.theta = single_cell_matrix.theta[mask, :]
        if !single_cell_matrix.ct_cov
            single_cell_matrix.sigma = single_cell_matrix.sigma[mask, :]
            single_cell_matrix.sigma = DataFrame(single_cell_matrix.sigma, Genes, celltypes)
        else
            single_cell_matrix.sigma_ct = single_cell_matrix.sigma_ct[:, mask]
            single_cell_matrix.sigma_ct = DataFrame(single_cell_matrix.sigma_ct, celltypes)
        end
        single_cell_matrix.DM = DataFrame(single_cell_matrix.DM, Genes, celltypes)
        single_cell_matrix.theta = DataFrame(single_cell_matrix.theta, Genes, celltypes)
    end

    println("Data preparation Done!")
end