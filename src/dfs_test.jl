using MatrixNetworks

type Biconnected_components_output
    bcc_edge_list::Vector{Int64}  #list of bcc edges
    components::Vector{Int64} #Indices of bcc corresponding to the edges in List
    number::Int64 # Total number of Biconnected components
    A::MatrixNetwork #MatrixNetwork
end

I = [1,2,2,3,3,4,4,5,5,6,6,7,7,1]
J = [2,1,3,2,4,3,5,4,6,5,7,6,1,7]
V = [1,1,1,1,1,1,1,1,1,1,1,1,1,1]

A = sparse(I,J,V)
# TODO : Reduce the number of function parameters
function BiCC_util(A::MatrixNetwork, vertex::Int64, bcc_edges::Vector{Int64}, low::Vector{Int64}, dt::Vector{Int64}, pred::Vector{Int64}, components::Vector{Int64}, cs::Vector{Int64}, cn::Int64, t::Int64, css::Int64, bcc_edges_count::Int64)
    # Initialize discovery time and low value
    t = t+1
    dt[vertex] = t
    low[vertex] = t
    children = 0
    for i=A.rp[vertex]:A.rp[vertex+1]-1
        adj = A.ci[i];
        # If adj is not visited yet, then recur for it
        if dt[adj] == 0
            children += 1
            pred[adj] = vertex
	    # store the edge in stack
            css = css+1
            cs[2*css-1] = vertex
            cs[2*css] = adj
            (cn,t,css,bcc_edges_count) = BiCC_util(A, adj, bcc_edges, low, dt, pred, components, cs, cn, t, css, bcc_edges_count) 
            # Check if the subtree rooted with 'vertex' has a connection to one of the ancestors of 'u'
	    if low[vertex] > low[adj]
                low[vertex] = low[adj]
            end
	    # If vertex is an articulation point, the pop all the edges till we encounter the edge vertex <---> adj  
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
	# Check for the back edges and update the low pointers 
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
    # If the vertex is not visited- call bcc utility on this vertex
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
