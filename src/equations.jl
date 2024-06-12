using DataFrames, LinearAlgebra, Graphs, Tables, Random, Statistics

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
    Tables.istable(t) || throw(ArgumentError("Argument does not support Tables.jl"))
    
    df = DataFrame(t)
    
    # Ensure all variable names are valid Julia symbols
    if !(all(occursin(r"^[A-Za-z_][A-Za-z0-9_]*$", string(name)) for name in names(df)))
        rename!(df, Symbol(name) => Symbol("var$(name)") for name in names(df))
    end

    # Check if it is a DAG
    if is_cyclic(est_g)
        throw(ArgumentError("The provided graph is cyclic -> est_g::DiGraph should be a DAG."))
    end

    adj_list = collect(edges(est_g))

    function ols(X, y)
        X = hcat(ones(size(X, 1)), X) 
        coef = X \ y
        yhat = X * coef
        resids = y - yhat
        return coef, resids
    end

    variables = String[]
    coefficients = Vector{Vector{Float64}}()
    residuals = Vector{Vector{Float64}}()
    nodes = names(df)

    for node in nodes
        println("node: ", node , " variables: ", variables)
        node_index = findfirst(==(node), nodes)
        
        preds = [nodes[e.src] for e in adj_list if e.dst == node_index]
        
        if !isempty(preds)
            X = hcat([df[!, pred] for pred in preds]...)
            y = df[!, node]
            
            coef, resid = ols(X, y)
            
            if isa(coef, Vector)
                push!(variables, string(node))
                push!(coefficients, coef)
                push!(residuals, resid)
            else
                println("Warning: Coefficients not stored for node $node. Expected vector, got $coef")
            end
        else
            y = df[!, node]
            intercept = mean(y)  
            resid = y .- intercept
            push!(variables, string(node))
            push!(coefficients, [intercept])
            push!(residuals, resid)
        end
    end

    return SCM(variables, coefficients, residuals, est_g)
end

# Function to generate data from the SCM
"""
    generate_data(scm::SCM, N::Int)::DataFrame

Generate data from the given SCM.

# Arguments
- `scm::SCM`: The structural causal model.
- `N::Int`: The number of data points to generate.

# Returns
- `DataFrame`: A DataFrame containing the generated data.
"""
function generate_data(scm::SCM, N::Int)::DataFrame
    df = DataFrame()
    
    sorted_indices = topological_sort_by_dfs(scm.dag)
    
    sorted_variables = [scm.variables[i] for i in sorted_indices]
    
    variable_index_map = Dict(variable => index for (index, variable) in enumerate(scm.variables))

    for node in sorted_variables
        idx = variable_index_map[node]
        coef = scm.coefficients[idx]
        residual_std = std(scm.residuals[idx])
        
        if length(coef) == 1
            df[!, Symbol(node)] = coef[1] .+ residual_std * randn(N)  
        else
            preds = [Symbol(scm.variables[i]) for i in inneighbors(scm.dag, idx)]
            if isempty(preds)
                df[!, Symbol(node)] = coef[1] .+ residual_std * randn(N)  
            else
                X = hcat(ones(N), [df[!, pred] for pred in preds]...)  
                df[!, Symbol(node)] = X * coef .+ residual_std * randn(N)
            end
        end
    end
    
    return df
end
