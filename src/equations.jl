using DataFrames, LinearAlgebra, Graphs

# Define the SCM struct
"""
    struct SCM
    variables::Vector{<:AbstractString}
    coefficients::Vector{<:Vector{<:AbstractFloat}}

A struct representing a Structural Causal Model (SCM).

# Fields
- `variables::Vector{<:AbstractString}`: A list of variable names.
- `coefficients::Vector{<:Vector{<:AbstractFloat}}`: A list of coefficient vectors for each variable.
"""
struct SCM 
    variables::Vector{<:AbstractString}
    coefficients::Vector{<:Vector{<:AbstractFloat}}
end

# Function to estimate equations and return an SCM struct
"""
    estimate_equations(df::DataFrame, est_g::DiGraph)::SCM

Estimate linear equations from the given DataFrame `df` based on the structure of the directed graph `est_g`.

# Arguments
- `df::DataFrame`: A DataFrame containing the data for estimation.
- `est_g::DiGraph`: A directed graph representing the structure of the SCM.

# Returns
- `SCM`: A struct containing the estimated variables and their corresponding coefficients.
"""
function estimate_equations(df::DataFrame, est_g::DiGraph)::SCM
    
    # Ensure all variable names are valid Julia symbols
    if !(all(occursin(r"^[A-Za-z_][A-Za-z0-9_]*$", string(name)) for name in names(df)))
        rename!(df, Symbol(name) => Symbol("var$(name)") for name in names(df))
    end

    # Ensure it is a DAG
    est_g = pdag2dag!(est_g) 

    adj_list = collect(edges(est_g))

    # Ordinary least squares function
    function ols(X, y)
        return inv(X' * X) * X' * y
    end

    variables = String[]
    coefficients = Vector{Vector{AbstractFloat}}()

    nodes = names(df)

    for node in nodes
        node_index = findfirst(==(node), nodes)
        
        preds = [e.src for e in adj_list if e.dst == node_index]
        
        if !isempty(preds)
            X = hcat([df[!, pred] for pred in preds]...)
            y = df[!, node]
            
            coef = ols(X, y)
            
            if isa(coef, Vector)
                push!(variables, string(node))
                push!(coefficients, coef)
            else
                println("Warning: Coefficients not stored for node $node. Expected vector, got $coef")
            end
        else
            push!(variables, string(node))
            push!(coefficients, AbstractFloat[])
        end
    end

    return SCM(variables, coefficients)
end

