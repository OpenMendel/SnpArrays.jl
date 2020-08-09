"""
    kinship_pruning(grm::Matrix; method=:bottom_up, cutoff=0.125, min_samples=0)

Return a BitVector indicating which samples remains after kinship pruning. 
Valid values for `method` are `:bottom_up`, `:top_down`, `:gcta`, and `:plink`.
When `min_samples` is defined, pruning stops if that number of samples remain.

- `:bottom_up` runs bottom-up greedy method, iteratively eliminating neighbors of sample with the lowest degree.
- `:top_down` runs top-down greedy method, iteratively eliminating sample with the highest degree.
- `:gcta` runs the method implemented in GCTA.
- `:plink` runs the method implemented in PLINK.

In general, `:bottom_up` greedy method performs the best (keeping the hightest number of subjects) and has a guaranteed approximamtion ratio.
With `min_samples` defined, `:top_down` greedy method keeps the least number of off-diagonal values.
"""
function kinship_pruning(grm::AbstractMatrix; method=:bottom_up, cutoff=0.125, min_samples=0)
    if method == :bottom_up
        r = _kinship_pruning_bottomup(grm; cutoff=cutoff, min_samples=min_samples)
    elseif method == :top_down
        r = _kinship_pruning_topdown(grm; cutoff=cutoff, min_samples=min_samples)
    elseif method == :gcta
        r = _kinship_pruning_gcta(grm; cutoff=cutoff)
    elseif method == :plink
        r = _kinship_pruning_plink(grm; cutoff=cutoff, min_samples=min_samples)
    else
        @error "Invalid value of argument method. Valid values are: :bottom_up, :top_down, :gcta, and :plink."
    end
    r
end

function _find_neighbors!(neighbors::Vector{Integer}, i::Integer, grm::AbstractMatrix, r::BitVector, degree::Integer; cutoff=0.125)
    m = size(grm, 1)
    n_neighbors = 0
    @inbounds for j in 1:m
        if r[j] && i != j && grm[i, j] > cutoff
            n_neighbors += 1
            neighbors[n_neighbors] = j
        end
    end
    @assert n_neighbors == degree
    neighbors
end

function _count_degrees(grm::AbstractMatrix; cutoff=0.125)
    m = size(grm, 1)
    degrees = zeros(Int, m)
    @inbounds for i in 1:m
        for j in 1:m
            if j == i; continue; end
            if grm[i, j] > cutoff
                degrees[i] += 1
            end
        end
    end
    degrees
end

function _eliminate_node!(i::Integer, grm::AbstractMatrix, r::BitVector, degrees::Vector{Integer}; cutoff=0.125)
    m = size(grm, 1)
    r[i] = false
    degrees[i] = 0
    @inbounds for j in 1:m
        if r[j] && grm[i, j] > cutoff
            degrees[j] -= 1
        end
    end
    nothing 
end

function _kinship_pruning_gcta(grm::AbstractMatrix; cutoff=0.125)
    m = size(grm, 1)
    degrees = _count_degrees(grm; cutoff=cutoff)
    
    r = trues(m)
    @inbounds for i in 1:m, j in (i+1):m
        if j == i; continue; end
        if grm[i, j] > cutoff
            if degrees[i] > degrees[j] 
                r[i] = false
            else
                r[j] = false
            end
        end
    end 
    r
end

function _kinship_pruning_plink(grm::AbstractMatrix; cutoff=0.125, min_samples=0)
    m = size(grm, 1)
    degrees = _count_degrees(grm; cutoff=cutoff) 
    r = trues(m) 
    nbs = [0]
    while sum(degrees) > 0 && count(r) > min_samples
        if any(degrees .== 1)
            to_add = findfirst(x -> x == 1, degrees)
            to_eliminate = _find_neighbors!(nbs, to_add, grm, r, degrees[to_add]; cutoff=cutoff)
            _eliminate_node!(to_eliminate[1], grm, r, degrees; cutoff=cutoff)
        else
           to_eliminate = findfirst(x -> x == maximum(degrees), degrees)
            _eliminate_node!(to_eliminate, grm, r, degrees; cutoff=cutoff) 
        end
    end
    if sum(degrees) != 0
        @warn "not all edges were removed: $(sum(degrees) ÷ 2) edges remain"
    end        
    r
end

function _kinship_pruning_topdown(grm::AbstractMatrix; cutoff=0.125, min_samples=0)
    m = size(grm, 1)
    degrees = _count_degrees(grm; cutoff=cutoff)
    r = trues(m)
    while sum(degrees) > 0 && count(r) > min_samples
        to_eliminate = findfirst(x -> x == maximum(degrees), degrees)
        _eliminate_node!(to_eliminate, grm, r, degrees; cutoff=cutoff)
    end
    if sum(degrees) != 0
        @warn "not all edges were removed: $(sum(degrees) ÷ 2) edges remain"
    end
    r
end

function _kinship_pruning_bottomup(grm::AbstractMatrix; cutoff=0.125, min_samples=0)
    m = size(grm, 1)
    degrees = _count_degrees(grm; cutoff=cutoff)
    nbs = Array{Int}(undef, m)
    r = trues(m)
    while sum(degrees) > 0 && count(r) > min_samples
        minpositive = minimum(degrees[degrees .> 0])
        to_add = findfirst(x -> x == minpositive, degrees)
        to_eliminate = _find_neighbors!(nbs, to_add, grm, r, degrees[to_add]; cutoff=cutoff) 
        @inbounds for i in 1:degrees[to_add]
            _eliminate_node!(to_eliminate[i], grm, r, degrees; cutoff=cutoff)
        end
    end
    if sum(degrees) != 0
        @warn "not all edges were removed: $(sum(degrees) ÷ 2) edges remain"
    end
    r
end
