using LinearAlgebra, Graphs, Tables, Random, Statistics

# Define the SCM struct
"""
    struct SCM
    variables::Vector{<:AbstractString}
    coefficients::Vector{<:Vector{<:AbstractFloat}}
    residuals::Vector{<:Vector{<:AbstractFloat}}
    dag::DiGraph

A struct representing a Structural Causal Model (SCM).

# Fields
- `variables::Vector{<:AbstractString}`: A list of variable names.
- `coefficients::Vector{<:Vector{<:AbstractFloat}}`: A list of coefficient vectors for each variable.
- `residuals::Vector{<:Vector{<:AbstractFloat}}`: A list of residuals for each variable.
- `dag::DiGraph`: The directed graph representing the structure of the SCM.
"""
struct SCM
    variables::Vector{String}
    coefficients::Vector{Vector{Float64}}
    residuals::Vector{Vector{Float64}}
    dag::DiGraph
end

function ols_compute(X, y)
    X = hcat(ones(size(X, 1)), X)
    coef = X \ y
    yhat = X * coef
    resids = y - yhat
    return coef, resids
end

# Function to estimate equations and return an SCM struct
"""
    estimate_equations(t, est_g::DiGraph)::SCM

Estimate linear equations from the given table `t` based on the structure of the directed graph `est_g`.

# Arguments
- `t`: A table containing the data for estimation (supports any Tables.jl-compatible format).
- `est_g::DiGraph`: A directed graph representing the structure of the SCM.

# Returns
- `SCM`: A struct containing the estimated variables, their corresponding coefficients, residuals, and the DAG.
"""
function estimate_equations(t, est_g::DiGraph)::SCM
    Tables.istable(t) || throw(ArgumentError("Argument supports just Tables.jl types"))
    
    columns = Tables.columns(t)
    schema = Tables.schema(t)
    variables = propertynames(schema.names)
    
    # Check if it is a DAG
    if is_cyclic(est_g)
        throw(ArgumentError("The provided graph is cyclic -> est_g::DiGraph should be a DAG."))
    end

    adj_list = collect(edges(est_g))

    var_names = String[]
    coefficients = Vector{Vector{Float64}}()
    residuals = Vector{Vector{Float64}}()
    nodes = variables

    for node in nodes
        node_index = findfirst(==(node), nodes)
        preds = [nodes[e.src] for e in adj_list if e.dst == node_index]

        if !isempty(preds)
            X = hcat([columns[pred] for pred in preds]...)
            y = columns[node]

            coef, resid = ols_compute(X, y)

            if isa(coef, Vector)
                push!(var_names, string(node))
                push!(coefficients, coef)
                push!(residuals, resid)
            else
                println("Warning: Coefficients not stored for node $node. Expected vector, got $coef")
            end
        else
            y = columns[node]
            intercept = mean(y)
            resid = y .- intercept
            push!(var_names, string(node))
            push!(coefficients, [intercept])
            push!(residuals, resid)
        end
    end

    return SCM(var_names, coefficients, residuals, est_g)
end

# Function to generate data from the SCM
"""
    generate_data(scm::SCM, N::Int)::NamedTuple

Generate data from the given SCM.

# Arguments
- `scm::SCM`: The structural causal model.
- `N::Int`: The number of data points to generate.

# Returns
- `NamedTuple`: A NamedTuple containing the generated data.
"""
function generate_data(scm::SCM, N::Int)::NamedTuple
    columns = Dict{Symbol, Vector{Float64}}()
    
    sorted_indices = topological_sort_by_dfs(scm.dag)
    sorted_variables = [scm.variables[i] for i in sorted_indices]
    variable_index_map = Dict(variable => index for (index, variable) in enumerate(scm.variables))

    for node in sorted_variables
        idx = variable_index_map[node]
        coef = scm.coefficients[idx]
        residual_std = std(scm.residuals[idx])
        
        if length(coef) == 1
            columns[Symbol(node)] = coef[1] .+ residual_std * randn(N)
        else
            preds = [Symbol(scm.variables[i]) for i in inneighbors(scm.dag, idx)]
            if isempty(preds)
                columns[Symbol(node)] = coef[1] .+ residual_std * randn(N)
            else
                X = hcat(ones(N), [columns[pred] for pred in preds]...)
                columns[Symbol(node)] = X * coef .+ residual_std * randn(N)
            end
        end
    end
    
    return NamedTuple(columns)
end