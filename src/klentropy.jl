using SpecialFunctions, NearestNeighbors, Distances, Distributions, Random

"""
    n_ball(n::Number)
Computes the volume of a n-dimensional unit sphere.
"""
function n_ball(n::Number)
    return π^(n/2.) / gamma(n/2. + 1.)
end


"""
    kl_entropy(data::Array{Float64, 2}; k=5)
Compute the nearest-neighbor estimate of the differential entropy of data.

data is a 2d array, with every column representing one data point. 
For further information, see

"A class of Rényi information estimators for multidimensional densities"
Nikolai Leonenko, Luc Pronzato, and Vippal Savani
The Annals of Statistics, 2008
https://projecteuclid.org/euclid.aos/1223908088

keyword arguments:
k=5: number of nearest neighbors
"""
function kl_entropy(data::Array{Float64, 2}; k=5)
    d, N = size(data)
    kdtree = KDTree(data)
    _, dist = knn(kdtree, data, k+1, true)
    H = log(N) - digamma(k) + log(n_ball(d)) + d/N * sum(map(l->log(l[end]), dist))
end


"""
    kl_mutual_information(x, y; k=5, bias_correction=true)
compute the nearest-neighbor 'KGS' estimate of the mutual information between x and y.

x and y are 2d arrays, with every column representing one data point. 
For further information, see

"Estimating Mutual Information"
Alexander Kraskov, Harald Stoegbauer, and Peter Grassberger
Physical Review E
https://arxiv.org/pdf/cond-mat/0305641.pdf

"Demystifying Fixed k-Nearest Neighbor Information Estimators"
Weihao Gao, Sewoong Oh, Pramod Viswanath
EEE International Symposium on Information Theory - Proceedings
https://arxiv.org/pdf/1604.03006.pdf

keyword arguments:
k=5: number of nearest neighbors
bias_correction=true: flag to apply Gao's bias correction
"""
function kl_mutual_information(x, y; k=5, bias_correction=true)
    dist = bias_correction ? Euclidean() : Chebyshev()
   
    d_x, N = size(x)
    d_y, _ = size(y)

    xy = vcat(x,y)

    kdtree_x = KDTree(x, dist)
    kdtree_y = KDTree(y, dist)
    kdtree_xy = KDTree(xy, dist)

    I = 0.
    
    for i in 1:N
        _, dist = knn(kdtree_xy, xy[:,i], k+1, true)
        n_x = length(inrange(kdtree_x, x[:,i], dist[end], false))-1
        n_y = length(inrange(kdtree_y, y[:,i], dist[end], false))-1
        I += digamma(n_x+1) + digamma(n_y+1)
    end

    I = -I/N + log(N) + digamma(k)

    if bias_correction
        I += log( n_ball(d_x) * n_ball(d_y)/n_ball(d_x+d_y) )
    end

    return I
end

"""
    kl_renyi(data::Array{Float64, 2}, q; k=5)
Compute the nearest-neighbor estimate of the Renyi-alpha entropy of data.

data is a 2d array, with every column representing one data point. 
For further information, see

"A class of Rényi information estimators for multidimensional densities"
Nikolai Leonenko, Luc Pronzato, and Vippal Savani
The Annals of Statistics, 2008
https://projecteuclid.org/euclid.aos/1223908088

keyword arguments:
k=5: number of nearest neighbors
"""
function kl_renyi(data, q, k=5)
    d, N = size(data)
    Vd = π^(d/2.) / gamma(d/2. + 1.)
    Ck = ( gamma(k)/gamma(k+1-q) )^(1/(1-q))
    kdtree = KDTree(data)
    _, dist = knn(kdtree, data, k+1, true)
    Iq = 1/N * sum(map(l->( (N-1)*Ck*Vd*l[end]^d )^(1-q), dist))

    return log(Iq)/(1-q)
end


"""
    kl_cond_mi(x, y, z; k=5, bias_correction=true)
compute the nearest-neighbor 'KGS' estimate of the conditional mutual information between x and y given z.

x, y, and z are 2d arrays, with every column representing one data point. 
keyword arguments:
k=5: number of nearest neighbors
bias_correction=true: flag to apply Gao's bias correction
"""
function kl_cond_mi(x, y, z; k=5, bias_correction=true)
    dist = bias_correction ? Euclidean() : Chebyshev()    
    
    xz = vcat(x,z)
    yz = vcat(y,z)
    xyz = vcat(x,y,z)

    d_x, N = size(x)
    d_y, _ = size(y)
    d_z, _ = size(z)
    
    kdtree_z = KDTree(z, dist)
    kdtree_xz = KDTree(xz, dist)
    kdtree_yz = KDTree(yz, dist)
    kdtree_xyz = KDTree(xyz, dist)
    
    CMI = 0.

    for i in 1:N
        _, dist = knn(kdtree_xyz, xyz[:,i], k+1, true)
        n_xz = length(inrange(kdtree_xz, xz[:,i], dist[end], false))-1
        n_yz = length(inrange(kdtree_yz, yz[:,i], dist[end], false))-1
        n_z = length(inrange(kdtree_z, z[:,i], dist[end], false))-1
        CMI = CMI + digamma(n_xz) + digamma(n_yz) - digamma(n_z)
    end

    CMI = digamma(k) - CMI/N

    if bias_correction
        CMI += log( (n_ball(d_x + d_z) * n_ball(d_y+d_z))/(n_ball(d_x+d_y+d_z)*n_ball(d_z)))
    end
    
    return CMI
end


"""
    kl_perm_mi_test(x, y; k=5, B=100, bias_correction=true)
compute permutation test of independence of x and y.

keyword arguments:
k=5: number of nearest neighbors to use for mutual information estimate
B=100: number of permutations
bias_correction=true: flag to apply Gao's bias correction
"""
function kl_perm_mi_test(x, y; k=5, B=100, bias_correction=true)
    MI = kl_mutual_information(x, y, k=k, bias_correction=bias_correction)
    samples = Float64[]

    for i in 1:B
        push!(samples, kl_mutual_information(x, y[:, shuffle(1:end)], k=k, bias_correction=bias_correction))
    end
    
    p = length(filter(d->MI<d, samples))/B
    return p
end


"""
    kl_perm_cond_mi_test(x, y, z; k=5, B=100, kp=5, bias_correction=true)
compute permutation test of conditional independence of x and y given z.

For further information, see:
"Conditional independence testing based on a nearest-neighbor estimator of conditional mutual information"
Jakob Runge
Proceedings of the 21st International Conference on Artificial Intelligence and Statistics (AISTATS) 2018, Lanzarote, Spain.
http://proceedings.mlr.press/v84/runge18a/runge18a.pdf

keyword arguments:
k=5: number of nearest neighbors to use for mutual information estimate
B=100: number of permutations
bias_correction=true: flag to apply Gao's bias correction

"""
function kl_perm_cond_mi_test(x, y, z; k=5, B=100, kp=5, bias_correction=true)
    d, N = size(z)
    CMI = kl_cond_mi(x,y,z,k=k,bias_correction=bias_correction)
    samples = Float64[]
    
    kdtree_z = KDTree(z, Chebyshev())
    z_knn, _ = knn(kdtree_z, z, kp+1)
    
    for b in 1:B
        # create permutation for independence test
        U = Int64[]
        P = collect(1:N)
        j = 0
        Ns = map(shuffle, z_knn)
        
        for i in shuffle(collect(1:N))
            j = Ns[i][1]
            m = 1
            while (j ∈ U) && (m < kp)
                m += 1
                j = Ns[i][m]
            end
            P[i]=j
            push!(U,j)
        end
        push!(samples, kl_cond_mi(x[:,P],y,z,k=k,bias_correction=bias_correction))
    end

    p = length(filter(d->CMI<d, samples))/B 
    
    return p
end
