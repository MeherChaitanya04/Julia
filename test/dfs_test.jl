using MatrixNetworks

type Biconnected_components_output
    bcc_edge_list::Vector{Int64}  #list of bcc edges
    components::Vector{Int64} #Indices of bcc corresponding to the edges in List
    number::Int64 # Total number of Biconnected components
    A::MatrixNetwork #MatrixNetwork
end

I = [1,2,2,5,5,1,5,3,3,6,6,5,4,5]
J = [2,1,5,2,1,5,3,5,6,3,5,6,5,4]
V = [1,1,1,1,1,1,1,1,1,1,1,1,1,1]

A = sparse(I,J,V)

function BiCC_util(A::MatrixNetwork, vertex::Int64, bcc_edges::Vector{Int64}, low::Vector{Int64}, dt::Vector{Int64}, pred::Vector{Int64}, components::Vector{Int64}, cs::Vector{Int64}, cn::Int64, t::Int64, css::Int64, bcc_edges_count::Int64)
    t = t+1
    dt[vertex] = t
    low[vertex] = t
    children = 0
    for i=A.rp[vertex]:A.rp[vertex+1]-1
        adj = A.ci[i];
        if dt[adj] == 0
            children += 1
            pred[adj] = vertex
            css = css+1
            cs[2*css-1] = vertex
            cs[2*css] = adj
            (cn,t,css,bcc_edges_count) = BiCC_util(A, adj, bcc_edges, low, dt, pred, components, cs, cn, t, css, bcc_edges_count) 
            if low[vertex] > low[adj]
                low[vertex] = low[adj]
            end
            if ((dt[vertex] == 1 && children > 1) || (dt[vertex] > 1 && low[adj] >= dt[vertex]))
                while css > 0
                    bcc_edges_count += 1
                    bcc_edges[2*bcc_edges_count-1] = cs[2*css-1]
                    bcc_edges[2*bcc_edges_count] = cs[2*css]
                    components[bcc_edges_count] = cn
                    if (cs[2*css-1] == vertex && cs[2*css] == adj)
                        css -= 1
                        break
                    end
                    css -= 1
                end
                cn += 1
            end
        elseif (adj != pred[vertex] && dt[adj] < low[vertex])
            low[vertex] = dt[adj]
            css = css+1
            cs[2*css-1] = vertex
            cs[2*css] = adj
        end
    end
    return (cn,t,css,bcc_edges_count) 
end

function BiCC(A::MatrixNetwork)
    bcc_edges=zeros(Int64,2*length(A.ci))
    bcc_edges_count=0
    cn=1
    low=zeros(Int64,A.n)
    dt=zeros(Int64,A.n)
    pred=zeros(Int64,A.n)
    components=-1 * ones(Int64,length(A.ci))
    t=0
    cs=zeros(Int64,2*length(A.ci))
    css=0
    for i=1:A.n
        if dt[i] == 0
            (cn,t,css,bcc_edges_count) = BiCC_util(A, i, bcc_edges, low, dt, pred, components, cs, cn, t, css, bcc_edges_count)
        end
    end 
    while css > 0
        bcc_edges_count += 1
        bcc_edges[2*bcc_edges_count-1] = cs[2*css-1]
        bcc_edges[2*bcc_edges_count] = cs[2*css]
	components[bcc_edges_count] = cn
        css -= 1
    end
    print(bcc_edges)
    print(components)
end

function dfs_test()
    print(A);
B = MatrixNetwork(A);
print(B) 
   print(BiCC(B))
end

dfs_test();
