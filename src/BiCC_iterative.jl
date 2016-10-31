
using MatrixNetworks


"""

Biconnected Component Decomposition
----------------------------------

Any connected graph decomposes into a tree of biconnected components called the block-cut tree of the graph. 
The blocks are attached to each other at shared vertices called cut vertices or articulation points.

This implementation is based on the algorithm provided by Tarjan in "Depth-First Search and Linear Graph Algorithms".  

TODO: 1) See how to output the biconnected components
      2) Compare its performance with C

"""
type Biconnected_components_output
    bcc_edge_list::Vector{Int64}  #list of bcc edges
    components::Vector{Int64} #Indices of bcc corresponding to the edges in List
    number::Int64 # Total number of Biconnected components
    A::MatrixNetwork #MatrixNetwork
end

#I = [1,2,2,3,3,4,4,5,5,1]
#J = [2,1,3,2,4,3,5,4,1,5]
#V = [1,1,1,1,1,1,1,1,1,1]

#I = [1,2,2,5,5,1,5,3,3,6,6,5,4,5]
#J = [2,1,5,2,1,5,3,5,6,3,5,6,5,4]
#V = [1,1,1,1,1,1,1,1,1,1,1,1,1,1]

#I = [1,2,2,3,3,1,4,5,5,6,4,6]
#J = [2,1,3,2,1,3,5,4,6,5,6,4]
#V = [1,1,1,1,1,1,1,1,1,1,1,1]

I = [1,2,2,3,3,4,1,4,3,5,5,6,6,7,7,8,6,8,5,9,9,10,5,10,1,3,2,5]
J = [2,1,3,2,4,3,4,1,5,3,6,5,7,6,8,7,8,6,9,5,10,9,10,5,3,1,5,2]
V = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]


A = sparse(I,J,V)



function BiCC(A::MatrixNetwork)
    n=A.n
    rp=A.rp
    ci=A.ci
    bcc_edges=zeros(Int64,2*length(ci))
    bcc_edges_count=0
    cn=1
    low=zeros(Int64,n)
    dt=-1*ones(Int64,n)
    pred=zeros(Int64,n)
    components=-1 * ones(Int64,length(ci))
    t=0
    cs=zeros(Int64,2*length(ci))
    css=0
    rs=zeros(Int64,4*n)
    rss=0 #recursion stack (v,ri)
    root_children= 0
    #start dfs at 1
    for sv=1:n
        v=sv
        if dt[v]>0 
            continue
        end
	root_children = 0
        rss=rss+1
        rs[3*rss-2]=v
	rs[3*rss-1]=v
        rs[3*rss]=rp[v] #add v,rp[v] to stack
        low[v]=t
        dt[v]=t;
        t=t+1;
        while rss>0
		grandparent = rs[3*rss-2]
		parent = rs[3*rss-1]
		children_index = rs[3*rss]
		child_index = children_index
		if (child_index != rp[parent+1])
			child = ci[child_index]
			children_index = children_index+1
			rs[3*rss] = children_index
			if (grandparent == child)		
				continue		# skip tree edges
			end
			if (dt[child] >= 0)
				if (dt[child] <= dt[parent])
					if (dt[child] < low[parent]) #Update Low
						low[parent] = dt[child]
					end
					css+=1
					cs[2*css-1] = parent
					cs[2*css] = child
				end
			else
				low[child] = dt[child] = t
				t+=1
				rss+=1
        			rs[3*rss-2]=parent
				rs[3*rss-1]=child
        			rs[3*rss]=rp[child]
				css+=1
				cs[2*css-1] = parent
				cs[2*css] = child
			end
		else
			# do StopIteration
			rss-=1
			print("here")
			if (rss > 1)
				if (low[parent] >= dt[grandparent])  # Grandparent is an Articulation point
					while true
						u = cs[2*css-1]
						v = cs[2*css]
						css-=1
						bcc_edges_count+=1
						bcc_edges[2*bcc_edges_count-1]=u
						bcc_edges[2*bcc_edges_count]=v
						components[bcc_edges_count]=cn   #Assign #Bcc to the corresponding edges
						if (u == grandparent && v == parent)
							break
						end
					end
					cn+=1
				end
				if (low[parent] < low[grandparent])
					low[grandparent] = low[parent]
				end
			elseif (rss == 1)
				root_children +=1
				while true
					u = cs[2*css-1]
					v = cs[2*css]
					css-=1
					bcc_edges_count+=1
					bcc_edges[2*bcc_edges_count-1]=u
					bcc_edges[2*bcc_edges_count]=v
					components[bcc_edges_count]=cn
					if (u == grandparent && v == parent)
						break
					end
				end
				cn+=1
			end
		end
	end
  end	
print(bcc_edges)
print("\n")
print(components)
   
end
 



function dfs_test()
		print(A);
		B = MatrixNetwork(A);
		print(B)
		print(BiCC(B))
end

dfs_test();


