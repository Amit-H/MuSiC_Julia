function Theta(single_cell_matrix::SingleCellMatrix, celltypes, samples)
    theta = []
    for ct in celltypes
        cross_subject = []
        for sid in samples
            mask_ct = single_cell_matrix.Meta[:, "cellType"] .== ct
            mask_sid = single_cell_matrix.Meta[:, "sampleID"] .== sid
            mask = mask_ct .& mask_sid
            columns = single_cell_matrix.Meta[mask, "sampleName"]
            sub_ct_sid = sum(select(single_cell_matrix.Count, columns), dims=2) / sum(select(single_cell_matrix.Count, columns))
            push!(cross_subject, vec(sub_ct_sid))
        end
        cross_subject = hcat(cross_subject...)
        cross_subject_mean = mean(cross_subject, dims=2)
        push!(theta, vec(cross_subject_mean))
    end
    theta = hcat(theta...)
    return theta
end

function Sigma(single_cell_matrix::SingleCellMatrix, celltypes, samples)
    if single_cell_matrix.ct_cov
        sigma_ct = []
        n_Genes = size(single_cell_matrix.Count, 1)
        sigma = nothing
        for ct in celltypes
            cross_subject = []
            for sid in samples
                mask_ct = single_cell_matrix.Meta[:, "cellType"] .== ct
                mask_sid = single_cell_matrix.Meta[:, "sampleID"] .== sid
                mask = mask_ct .& mask_sid
                columns = single_cell_matrix.Meta[mask, "sampleName"]
                sub_ct_sid = vec(sum(select(single_cell_matrix.Count, columns), dims=2) / sum(select(single_cell_matrix.Count, columns)))
                cross_subject = vcat(cross_subject, sub_ct_sid)
            end
            if sigma === nothing
                sigma = cross_subject
            else
                sigma = vcat(sigma, cross_subject)
            end
        end
        n_subs = length(samples)
        sigma_ct = []
        for i in 1:size(single_cell_matrix.Count, 1)
            index = i .+ n_Genes .* (0:n_subs-1)
            sub_M = sigma[index, :]
            sub_cov = cov(Matrix(sub_M))
            push!(sigma_ct, vec(sub_cov))
        end
        sigma_ct = hcat(sigma_ct...)
        return sigma_ct
    else
        Sigma = []
        for ct in celltypes
            cross_subject = []
            for sid in samples
                mask_ct = single_cell_matrix.Meta[:, "cellType"] .== ct
                mask_sid = single_cell_matrix.Meta[:, "sampleID"] .== sid
                mask = mask_ct .& mask_sid
                columns = single_cell_matrix.Meta[mask, "sampleName"]
                sub_ct_sid = sum(select(single_cell_matrix.Count, columns), dims=2) / sum(select(single_cell_matrix.Count, columns))
                push!(cross_subject, vec(sub_ct_sid))
            end
            cross_subject = hcat(cross_subject...)
            cross_subject_var = var(cross_subject, dims=2, corrected=false)
            push!(Sigma, vec(cross_subject_var))
        end
        Sigma = hcat(Sigma...)
        return Sigma
    end
end

function LibrarySize(single_cell_matrix::SingleCellMatrix, celltypes, samples)
    LS = []
    for ct in celltypes
        cross_subject = []
        for sid in samples
            mask_ct = single_cell_matrix.Meta[:, "cellType"] .== ct
            mask_sid = single_cell_matrix.Meta[:, "sampleID"] .== sid
            mask = mask_ct .& mask_sid
            columns = single_cell_matrix.Meta[mask, "sampleName"]
            sub_ct_sid = sum(select(single_cell_matrix.Count, columns))
            push!(cross_subject, mean(sub_ct_sid))
        end
        push!(LS, cross_subject)
    end
    LS = hcat(LS...)
    LS = DataFrame(LS, names(samples), celltypes)
    return LS
end

struct NNLS
    """
    Simple wrapper for NNLS in Julia.
    """
    X::Matrix{Float64}
    y::Vector{Float64}
    coef::Vector{Float64}
    score::Float64
    r::Vector{Float64}

    function NNLS(X, y)
        new(X, y, zeros(size(X, 2)), 0.0, zeros(size(y)))
        fit!(this)
        resid!(this)
    end
end

function fit!(nnls_model::NNLS)
    """
    Fit the NNLS model.
    """
    nnls_model.coef, nnls_model.score = nnls(nnls_model.X, nnls_model.y)
end

function resid!(nnls_model::NNLS)
    """
    Calculate residuals.
    """
    nnls_model.r = nnls_model.y - nnls_model.X * nnls_model.coef
end

struct ModelResults
    initial_p::DataFrame
    initial_coef::DataFrame
    initial_resid::DataFrame
    weight::DataFrame
    weight_p::DataFrame
    weight_coef::DataFrame
    weight_resid::DataFrame
    R2::Vector{Float64}
    var_p::DataFrame
end

function weight_cal(Sp, Sigma)
    """
    Calculate weight with cross-subject variance for each cell type.
    """
    return sum(Sp.^2 .* Sigma, dims=2)
end

function weight_cal_ct(Sp, Sigma_ct)
    """
    Calculate weight with cross cell type covariance.
    """
    nGenes = size(Sigma_ct, 2)
    n_ct = length(Sp)
    Sp = vec(Sp)
    Sp_2 = Sp .* reshape(Sp, (length(Sp), 1))
    
    weights = []
    
    for i in 1:nGenes
        sig = reshape(Sigma_ct[:, i], (n_ct, n_ct))
        push!(weights, sum(Sp_2 .* sig))
    end
    
    return weights
end

function music_basic(Y, X, S, sigma, iter_max, nu, eps, centered, normalize)
    """
    Weight is estimated with cell type variance.
    """
    if centered
        X = X - mean(X)
        Y = Y - mean(Y)
    end

    if normalize
        X = X / std(X)
        S = S * std(S)
        Y = Y / std(Y)
    else
        Y = Y * 100
    end

    lm_D = NNLS(X, Y)
    r = lm_D.r
    weight_gene = 1.0 / (nu + r.^2 + weight_cal(Sp=lm_D.coef .* S, Sigma=sigma))
    Y_weight = Y .* sqrt.(weight_gene)
    X_weight = X .* sqrt.(weight_gene)
    lm_D_weight = NNLS(X_weight, Y_weight)

    p_weight = lm_D_weight.coef / sum(lm_D_weight.coef)
    r = lm_D_weight.r

    println("    Iteraction NNLS start ...")
    for i in 1:iter_max
        weight_gene = 1.0 / (nu + r.^2 + weight_cal(Sp=lm_D_weight.coef .* S, Sigma=sigma))
        Y_weight = Y .* sqrt.(weight_gene)
        X_weight = X .* sqrt.(weight_gene)
        lm_D_weight = NNLS(X_weight, Y_weight)
        p_weight_new = lm_D_weight.coef / sum(lm_D_weight.coef)
        r_new = lm_D_weight.r

        if sum(abs.(p_weight_new - p_weight)) < eps
            println("    Done, save results.")
            p_weight = p_weight_new
            r = r_new
            R2 = 1.0 - var(Y - X * lm_D_weight.coef) / var(Y)
            var_p = inv(X_weight' * X_weight) * mean(r.^2) / (sum(lm_D_weight.coef)^2)

            return ModelResults(p_weight, lm_D.coef, lm_D.r, weight_gene,
                                p_weight, lm_D_weight.coef, r_new, [R2], var_p)
        end

        p_weight = p_weight_new
        r = r_new
    end

    println("    Done, save results.")
    R2 = 1.0 - var(Y - X * lm_D_weight.coef) / var(Y)
    var_p = inv(X_weight' * X_weight) * mean(r.^2) / (sum(lm_D_weight.coef)^2)

    return ModelResults(p_weight, lm_D.coef, lm_D.r, weight_gene,
                        p_weight, lm_D_weight.coef, r_new, [R2], var_p)
end

function music_basic_ct(Y, X, S, sigma_ct, iter_max, nu, eps, centered, normalize)
    """
    Weight is estimated with cell type covariance.
    """
    if centered
        X = X - mean(X)
        Y = Y - mean(Y)
    end

    if normalize
        X = X / std(X)
        S = S * std(S)
        Y = Y / std(Y)
    else
        Y = Y * 100
    end

    lm_D = NNLS(X, Y)
    weights = weight_cal_ct(Sp=lm_D.coef .* S, Sigma_ct=sigma_ct)
    weight_gene = 1.0 / (nu + lm_D.r.^2 + weights)
    Y_weight = Y .* sqrt.(weight_gene)
    X_weight = X .* sqrt.(weight_gene)

    lm_D_weight = NNLS(X_weight, Y_weight)

    p_weight = lm_D_weight.coef / sum(lm_D_weight.coef)
    r = lm_D_weight.r

    println("    Iteraction NNLS start ...")
    for i in 1:iter_max
        weight_gene = 1.0 / (nu + r.^2 + weight_cal_ct(lm_D_weight.coef .* S, Sigma_ct=sigma_ct))
        Y_weight = Y .* sqrt.(weight_gene)
        X_weight = X .* sqrt.(weight_gene)
        lm_D_weight = NNLS(X_weight, Y_weight)
        p_weight_new = lm_D_weight.coef / sum(lm_D_weight.coef)
        r_new = lm_D_weight.r

        if sum(abs.(p_weight_new - p_weight)) < eps
            println("    Done, save results.")
            p_weight = p_weight_new
            r = r_new
            R2 = 1.0 - var(Y - X * lm_D_weight.coef) / var(Y)
            var_p = inv(X_weight' * X_weight) * mean(r.^2) / (sum(lm_D_weight.coef)^2)

            return ModelResults(p_weight, lm_D.coef, lm_D.r, weight_gene,
                                p_weight, lm_D_weight.coef, r_new, [R2], var_p)
        end

        p_weight = p_weight_new
        r = r_new
    end

    println("    Done, save results.")
    R2 = 1.0 - var(Y - X * lm_D_weight.coef) / var(Y)
    var_p = inv(X_weight' * X_weight) * mean(r.^2) / (sum(lm_D_weight.coef)^2)

    return ModelResults(p_weight, lm_D.coef, lm_D.r, weight_gene,
                        p_weight, lm_D_weight.coef, r_new, [R2], var_p)
end

function music_prop(bulk_Meta, bulk_Count, sc_Meta, sc_Count, markers, clusters, samples,
                    select_ct, ct_cov, iter_max, nu, eps, centered, normalize)
    """
    MuSiC Deconvolution.
    """
    # Data preprocessing
    Nz_genes = mean(bulk_Count, dims=2) .!= 0
    bulk_Count = bulk_Count[Nz_genes, :]
    bulk_gene = bulk_Count.index
    
    if markers == nothing
        sc_markers = bulk_gene
    else
        sc_markers = intersect(markers, bulk_gene)
    end
    
    sc_basis = SingleCellMatrix(sc_Meta, sc_Count, markers=sc_markers, select_ct=select_ct, ct_cov=ct_cov)
    
    cm_gene = sc_basis.DM.index
    if markers == nothing
        if length(cm_gene) < 0.2 * min(length(bulk_gene), length(sc_Count.index))
            throw(ArgumentError("Too few common genes!"))
        end
    else
        if length(cm_gene) < 0.2 * length(Set(markers))
            throw(ArgumentError("Too few common genes!"))
        end
    end
    
    println("Used ", length(cm_gene), " common genes ...")
    
    bulk_Count = bulk_Count[cm_gene, :]
    Yjg = BulkMatrix(convert(Matrix, bulk_Count)).X
    Yjg = DataFrame(Yjg, index = bulk_Count.index, columns = bulk_Count.columns)
    N_bulks = size(Yjg, 2)
    S = sc_basis.MLS
    
    # Deconvolution
    if ct_cov
        Results = ModelResults(DataFrame(), DataFrame(), DataFrame(),
                               DataFrame(), DataFrame(), DataFrame(),
                               DataFrame(), Float64[], DataFrame())
        
        for i in 1:N_bulks
            
            # Remove zero gene in bulk sample i
            Y_i = Yjg[:, i]
            mask = Y_i .!= 0
            Y_i = Y_i[mask]
            D_i = sc_basis.DM[mask, :]
            sigma_ct_i = sc_basis.sigma_ct[:, mask]
            println("Sample ", bulk_Count.columns[i], " has ", length(Y_i), " available genes")
            
            sample_results = music_basic_ct(Y_i, D_i, S,
                                            sigma_ct=sigma_ct_i,
                                            iter_max=iter_max,
                                            nu=nu, eps=eps,
                                            centered=centered,
                                            normalize=normalize)
            
            DataFrame!(sample_results.initial_p, index = bulk_Count.columns, columns = sc_basis.DM.columns)
            DataFrame!(sample_results.initial_coef, index = bulk_Count.columns, columns = sc_basis.DM.columns)
            DataFrame!(sample_results.initial_resid, index = bulk_Count.columns, columns = cm_gene)
            DataFrame!(sample_results.weight, index = cm_gene, columns = bulk_Count.columns)
            DataFrame!(sample_results.weight_p, index = bulk_Count.columns, columns = sc_basis.DM.columns)
            DataFrame!(sample_results.weight_coef, index = bulk_Count.columns, columns = sc_basis.DM.columns)
            DataFrame!(sample_results.weight_resid, index = bulk_Count.columns, columns = cm_gene)
            append!(Results.R2, sample_results.R2[1])
            DataFrame!(sample_results.var_p, index = bulk_Count.columns, columns = sc_basis.DM.columns)
        end

    else
        Results = ModelResults(DataFrame(), DataFrame(), DataFrame(),
                               DataFrame(), DataFrame(), DataFrame(),
                               DataFrame(), Float64[], DataFrame())
        
        for i in 1:N_bulks
            
            # Remove zero gene in bulk sample i
            Y_i = Yjg[:, i]
            mask = Y_i .!= 0
            Y_i = Y_i[mask]
            D_i = sc_basis.DM[mask, :]
            sigma_i = sc_basis.sigma[mask, :]
            println("Sample ", bulk_Count.columns[i], " has ", length(Y_i), " available genes")
            sample_results = music_basic(Y_i, D_i, S,
                                         sigma=sigma_i,
                                         iter_max=iter_max,
                                         nu=nu, eps=eps,
                                         centered=centered,
                                         normalize=normalize)
            
            DataFrame!(sample_results.initial_p, index = bulk_Count.columns, columns = sc_basis.DM.columns)
            DataFrame!(sample_results.initial_coef, index = bulk_Count.columns, columns = sc_basis.DM.columns)
            DataFrame!(sample_results.initial_resid, index = bulk_Count.columns, columns = cm_gene)
            DataFrame!(sample_results.weight, index = cm_gene, columns = bulk_Count.columns)
            DataFrame!(sample_results.weight_p, index = bulk_Count.columns, columns = sc_basis.DM.columns)
            DataFrame!(sample_results.weight_coef, index = bulk_Count.columns, columns = sc_basis.DM.columns)
            DataFrame!(sample_results.weight_resid, index = bulk_Count.columns, columns = cm_gene)
            append!(Results.R2, sample_results.R2[1])
            DataFrame!(sample_results.var_p, index = bulk_Count.columns, columns = sc_basis.DM.columns)
        end
    end
    
    return Results
end
