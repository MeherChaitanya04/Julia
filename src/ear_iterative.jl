using MatrixNetwork

"""
DFS
---
    compute depth first search distances and returns the distance (d), the discover (dt),
    the finish time(ft), and the predecessor array (pred) in the tuple (d,dt,ft, pred).\n
    pred[i] = 0 if vertex i is in a component not reachable from u and i != u.

Functions
---------
- (d,dt,ft,pred) = dfs(A::MatrixNetwork,u::Int64,full::Int64,target::Int64)
- (d,dt,ft,pred) = dfs(A::MatrixNetwork,u::Int64)
- (d,dt,ft,pred) = dfs{T}(A::SparseMatrixCSC{T,Int64},u::Int64,full::Int64,target::Int64)
- (d,dt,ft,pred) = dfs{T}(A::SparseMatrixCSC{T,Int64},u::Int64)
- (d,dt,ft,pred) = dfs(ei::Vector{Int64},ej::Vector{Int64},u::Int64,full::Int64,target::Int64)
- (d,dt,ft,pred) = dfs(ei::Vector{Int64},ej::Vector{Int64},u::Int64)

Example
-------
~~~
A = load_matrix_network("dfs_example")
(d,dt,ft,pred)  = dfs(A,1)
~~~
"""

function ear(A::MatrixNetwork,u::Int64,full::Int64,target::Int64)

    (rp,ci) = (A.rp,A.ci)
    n=length(rp)-1
    d=-1*ones(Int64,n)
    dt=-1*ones(Int64,n)
    ft=-1*ones(Int64,n)
    pred=zeros(Int64,n)
    rs=zeros(Int64,2*n)
    rss=0 # recursion stack holds two nums (v,ri)
    be = zeros(Int64, (length(A.ci) - A.n +1), 3) #back edges array   
    be_index=0; 
    ears = zeros(Int64,2 * length(A.ci))
    ear_number = zeros(Int64,length(A.ci))
    ear_index = 0
    ear_visited = zeros(Int64, A.n)
    #start dfs at u
    t=0
    targethit=0
    for i=1:n
        if i==1
            v=u
        else
            v=mod(u+i-1,n)+1
            if d[v]>0
                continue
            end
        end
        d[v]=0
        dt[v]=t
        t=t+1
        ri=rp[v]
        rss=rss+1
        rs[2*rss-1]=v
        rs[2*rss]=ri # add v to the stack
        
        while rss>0
            v=rs[2*rss-1]
            ri=rs[2*rss]
            rss=rss-1 # pop v from the stack
            if v==target || targethit == 1
                ri=rp[v+1]
                targethit=1 # end the algorithm if v is the target
            end
            while ri<rp[v+1]
                w=ci[ri]
                ri=ri+1
                if d[w]<0
                    d[w]=d[v]+1
                    pred[w]=v
                    rss=rss+1
                    rs[2*rss-1]=v
                    rs[2*rss]=ri # add v to the stack
                    v=w
                    ri=rp[w]
                    dt[v]=t
                    t=t+1
                    continue #discover a new vertex!
                elseif (pred[v]!=w)
		    be_index+=1
                    be[be_index][1]=dt[w]; be[be_index][2]=w; be[be_index][3]=v;
		end
			
            end
            ft[v]=t
            t=t+1 # finish with v
        end
        if full == 0
            break
        end
    end
    firstid = be[:,1]
    sortids = sortperm(firstid);
    cn=1
    for i=1:length(sortids)
        startVertex = be[sortids[i]][2]
        ear_visited[startVertex] = 1
        nextVertex = be[sortids[i]][3]
        ear_index += 1
        ears[2*ear_index-1] = startVertex
        ears[2*ear_index] = nextVertex
        ear_number[ear_index] = cn
        while (ear_visited[nextVertex] != 1)
            ear_visited[nextVertex] = 1
            startVertex = nextVertex
            nextVertex = pred[startVertex]    
            ear_index += 1
            ears[2*ear_index-1] = startVertex
            ears[2*ear_index] = nextVertex
            ear_number[ear_index] = cn
        end 
        cn += 1
    end
    #be_sorted = be[sortids,:]
    print(ears)
    print(ear_number)
    return (d,dt,ft,pred)
end

############################
### Additional functions ###
############################

function dfs(A::MatrixNetwork,u::Int64)
    return dfs(A,u,0,0);
end

## CSC sparse matrices:
function dfs{T}(A::SparseMatrixCSC{T,Int64},u::Int64)
    return dfs(MatrixNetwork(A),u)
end

function dfs{T}(A::SparseMatrixCSC{T,Int64},u::Int64,full::Int64,target::Int64)
    return dfs(MatrixNetwork(A),u,full,target)
end

## Triplet Format:
function dfs(ei::Vector{Int64},ej::Vector{Int64},u::Int64)
    return dfs(MatrixNetwork(ei,ej),u)
end

function dfs(ei::Vector{Int64},ej::Vector{Int64},u::Int64,full::Int64,target::Int64)
    return dfs(MatrixNetwork(ei,ej),u,full,target)
end

## CSR sparse matrices:
function dfs{T}(rp::Vector{Int64},ci::Vector{Int64},vals::Vector{T},n::Int64,u::Int64)
    return dfs(MatrixNetwork(n,rp,ci,vals),u)
end

function dfs{T}(rp::Vector{Int64},ci::Vector{Int64},vals::Vector{T},n::Int64,u::Int64,full::Int64,target::Int64)
    return dfs(MatrixNetwork(n,rp,ci,vals),u,full,target)
end
