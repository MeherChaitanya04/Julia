using MatrixNetworks

I = [1,2,3,4,5,4]
J = [2,3,1,5,6,6]
V = [1,1,1,1,1,1]

A = sparse(I,J,V)
function dfs(A::MatrixNetwork,u::Int64,full::Int64,target::Int64)

    (rp,ci) = (A.rp,A.ci)
    n=length(rp)-1
    d=-1*ones(Int64,n)
    dt=-1*ones(Int64,n)
    ft=-1*ones(Int64,n)
    pred=zeros(Int64,n)
    rs=zeros(Int64,2*n)
    rss=0 # recursion stack holds two nums (v,ri)
    
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
                end
            end
            ft[v]=t
            t=t+1 # finish with v
        end
        if full == 0
            break
        end
    end
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


function dfs_test()
#print(A);
B = MatrixNetwork(A);
#print(B) 
dfs(B, 1,1,0)
end

dfs_test();
