function forward_step(dag)
    # N is number of variables
    #score has to be initialised

    # record the local score calculated each time. Thus when we transition to the second phase,
    # many of the operators can be scored without an explicit call the the scoring function
    # record_local_score{trial}{j} record the local scores when Xj as a parent 
    record_local_score = [[] for i in 1:n]

    score_new = score
    count1 = 0
    update1 = []
    G_step1 = []
    score_record1 = []
    graph_record1 = []
    while true
        count1 = count1 + 1
        score = score_new
        append!(score_record1, score)
        append!(graph_record1, cpdag)
        min_chscore = 1e7
        min_desc = []
        for s in ve(cpdag)
            for d in ve(cpdag)
                # find a pair (s,d) that is not adjacent in the current graph , and restrict the number of parents
                if !has_edge(g, s, d) && !has_edge(g, d, s) && length(inneighbors(g, d) <= maxP)
                    Td = neighbors(d)
                    Ts = union(inneighbors(g,s),outneighbors(g,s))
                    NTs = setdiff(ve(cpdag), Ts) 
                    T0 = intersect(Td, NTs)
                    # for any subset of T0
                    sub = powerset(T0)
                    S = zeros(length(sub))
                    # S indicate whether we need to check sub{k}.
                    # 0: check both conditions.
                    # 1: only check the first condition
                    # 2: check nothing and is not valid.
                    for k in 1:length(sub)
                        if (S[k] < 2) #S indicates whether we need to check subset(k)
                            V1 = 

                    # find the neighbours of Xj that are not adjacent to Xi



isadjacent(dg, v, w) = has_edge(dg, v, w) || has_edge(dg, w, v)

function insert_validity_test1(G, s, d, T)
    # V=Insert_validity_test1(G, X, Y, T,1); % do validity test for the operator Insert; V=1 means valid, V=0 mean invalid;
    # here G is CPDAG
    V = 0

    # condition 1
    Td = neighbors(g, s) # neighbors of Xj
    Ts = union(inneighbors(g,s),outneighbors(g,s))
    NA = np.intersect1d(Tj, Ti)  # find the neighbours of Xj and are adjacent to Xi
    V = check_clique(G, list(np.union1d(NA, T).astype(int)))  # check whether it is a clique
    return 
end



end

function backward_step(g)
end

function turning_step(g)
end

"""
    apply_ges_rules(g, dg)


`g` is an empty graph.
`dg` is an directed graph containing edges of known direction and 
both `v=>w` and `w=>v `if the direction of `Edge(v,w)`` is unknown.

Returns the CPDAG as DiGraph. 
"""  
function apply_ges_rules(g, dg, X, score; VERBOSE = false)
    score = 

function gesalg(t, score, maxP)
    @assert Tables.istable(t)
    c = Tables.columns(t)
    sch = Tables.schema(t)
    n = length(sch.names)



    if (X.shape[0] < X.shape[1]):
        warnings.warn("The number of features is much larger than the sample size!")