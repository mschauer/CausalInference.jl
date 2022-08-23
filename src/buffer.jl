function pdag2dag(g, dg)
    """
    Covert a PDAG to its corresponding DAG
    Parameters
    ----------
    G : Partially Direct Acyclic Graph
    Returns
    -------
    Gd : Direct Acyclic Graph
    """
    nodes = vertices(g)
    # first create a DAG that contains all the directed edges in PDAG
    cg = copy(g)
    Gp = deepcopy(G)
    inde = zeros(nv(Gp))  # index whether the ith node has been removed. 1:removed; 0: not
    while 0 in inde
        for i in 1:nv(Gp)
            if inde[i] == 0
                sign = 0
                if (length(np.intersect1d(np.where(Gp.graph[:, i] == 1)[0],
                                       np.where(inde == 0)[0])) == 0):  # Xi has no out-going edges
                    sign = sign + 1
                    Nx = np.intersect1d(
                        np.intersect1d(np.where(Gp.graph[:, i] == -1)[0], np.where(Gp.graph[i, :] == -1)[0]),
                        np.where(inde == 0)[0])  # find the neighbors of Xi in P
                    Ax = np.intersect1d(np.union1d(np.where(Gp.graph[i, :] == 1)[0], np.where(Gp.graph[:, i] == 1)[0]),
                                        np.where(inde == 0)[0])  # find the adjacent of Xi in P
                    Ax = np.union1d(Ax, Nx)
                    if len(Nx) > 0:
                        if check2(Gp, Nx, Ax):  # according to the original paper
                            sign = sign + 1
                    else:
                        sign = sign + 1
                if sign == 2:
                    # for each undirected edge Y-X in PDAG, insert a directed edge Y->X in G
                    for index in np.intersect1d(np.where(Gp.graph[:, i] == -1)[0], np.where(Gp.graph[i, :] == -1)[0]):
                        Gd.add_edge(Edge(nodes[index], nodes[i], Endpoint.TAIL, Endpoint.ARROW))
                    inde[i] = 1

    return Gd



function check2(g, dg, Nx, Ax):
    s = 1
    for i in 1:length(Nx)
        j = np.delete(Ax, np.where(Ax == Nx[i])[0])
        g
        if len(np.where(G.graph[Nx[i], j] == 0)[0]) != 0:
            s = 0
            break
    return s