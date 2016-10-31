using MatrixNetworks

type Biconnected_components_output
    bcc_edge_list::Vector{Int64}  #list of bcc edges
    components::Vector{Int64} #Indices of bcc corresponding to the edges in List
    number::Int64 # Total number of Biconnected components
    A::MatrixNetwork #MatrixNetwork
end

I = [1,2,3,2,3,1,3,4,4,5,5,3]
J = [2,1,2,3,1,3,4,3,5,4,3,5]
V = [1,1,1,1,1,1,1,1,1,1,1,1]

A = sparse(I,J,V)

function BiCC(A::MatrixNetwork)
    n=A.n
    rp=A.rp
    ci=A.ci
    bcc_edges=zeros(Int64,2*length(ci))
    bcc_edges_count=0
    cn=1
    low=zeros(Int64,n)
    dt=zeros(Int64,n)
    pred=zeros(Int64,n)
    components=-1 * ones(Int64,length(ci))
    t=0
    cs=zeros(Int64,2*length(ci))
    css=0
    rs=zeros(Int64,2*n)
    rss=0 #recursion stack (v,ri)
    children = zeros(Int64, n)
    #start dfs at 1
    for sv=1:n
        v=sv
        if dt[v]>0 
            continue
        end
        rss=rss+1
        rs[2*rss-1]=v
        rs[2*rss]=rp[v] #add v,rp[v] to stack
        low[v]=t
        dt[v]=t;
        t=t+1;
        
        while rss>0
            v=rs[2*rss-1]
            ri=rs[2*rss]
            rss=rss-1
            while ri<rp[v+1]
                w=ci[ri]
                ri=ri+1
                if dt[w] == 0
                    low[w]=w; dt[w]=t; t=t+1;
                    css=css+1    
                    cs[2*css-1]=v
                    cs[2*css]=w
                    children[v]=children[v]+1
                    pred[w]=v
                    rss=rss+1
                    rs[2*rss-1]=v
                    rs[2*rss]=ri
                    v=w
                    ri=rp[w]
                    continue
                elseif(w != pred[v] && dt[w] < low[v])
 		    if low[v] > dt[w]
			low[v] = dt[w]
		    end
                    css=css+1    
                    cs[2*css-1]=v
                    cs[2*css]=w
                end
            end
             for ri=rp[v]:rp[v+1]-1
                w=ci[ri]
                if components[w] == -1
                    if low[v] > low[w]
                        low[v]=low[w]
                    end
                end
                if (( (dt[v] == 1 ) && ( children[v] > 1 ) ) || ( dt[v] > 1 && low[w] >= dt[v]))
                    while css > 0 
                        bcc_edges_count = bcc_edges_count + 1
                        bcc_edges[2*bcc_edges_count-1] = cs[2*css-1]
                        bcc_edges[2*bcc_edges_count] = cs[2*css]
                        components[bcc_edges_count] = cn
			if (cs[2*css-1] == v && cs[2*css] ==w)
				css=css-1
				break
			end
                        css = css-1
                    end
			
                    cn=cn+1
                end
            end    
                  
          end
	end
	print(components)
          return bcc_edges
end

function dfs_test()
    print(A);
B = MatrixNetwork(A);
#print("\n");print("Matrix B is ")
#print(B) 
print("\n")
 print(BiCC(B))
end

dfs_test();
