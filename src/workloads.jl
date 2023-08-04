using PrecompileTools

@setup_workload begin
    # define simple linear model with added noise
    df___ = let N = 200 # number of data points
        x = randn(N)
        v = x + randn(N)*0.25
        w = x + randn(N)*0.25
        z = v + w + randn(N)*0.25
        s = z + randn(N)*0.25
        (x=x, v=v, w=w, z=z, s=s)
    end
    @compile_workload begin
        ges(df___; penalty=1.0, parallel=false)
        pcalg(df___, 0.01, gausscitest; stable=false)
    end

    dag__ = digraph([1 => 3, 3 => 6, 2 => 5, 5 => 8, 6 => 7, 7 => 8, 1 => 4, 2 => 4, 4 => 6, 4 => 8])

    @compile_workload begin
        dsep(dag__, 6, 8, 1)
        collect(list_covariate_adjustment(dag__, 6, 8, Int[], setdiff(Set(1:8), [1, 2])))
    end
end