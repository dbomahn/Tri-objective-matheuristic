# Feasibility-pump and shared solution-processing routines.
# Extracted from phd/codes/FPGPRfunctions.jl; see docs/source_audit.md.

function flip(x_h,j,e)
    if x_h[e[j]]==1
        x_h[e[j]] = 0
    else
        x_h[e[j]] = 1
    end
    return x_h
end
function flipoper(Tabu,x_t,x_r)
    e = sortperm(abs.(x_t-x_r),rev=true)
    xi = []
    x_h = copy(x_r)
    j = 1
    M=length(x_t) #
    while j<=M && xi==[]
        x_h = flip(x_h,j,e)
        if x_h ∉ Tabu
            xi=x_h
        else
            j+=1
        end
    end
    if xi==[]
        while j<=M
            x_h=copy(x_r)
            Num = Int64(rand(ceil(length(x_r)/2):length(x_r)-1))
            R = sample(1:M,Num, replace=false)
            for i in R
                x_h = flip(x_h,r)
                if x_h ∉ Tabu
                    xi = x_h
                end
            end
            j+=1
        end
    end
    return xi
end
function FBcheck(xx,n)
    for k=1:n
        JuMP.fix(mx[k],xx[k]; force=true)
    end
    optimize!(mip)
    if termination_status(mip) == MOI.OPTIMAL
        return true
    else
        return false
    end
end
function fbsearch(x_r) #solveLP
    idx0 = findall(k->k==0, x_r)
    idx1 = findall(k->k==1, x_r)
    @objective( dist, Min, sum(x[i] for i in idx0) + sum(1-(x[j]) for j in idx1) )
    optimize!(dist)
    if termination_status(dist) == MOI.OPTIMAL
        return JuMP.value.(x)
    else
        return 0;
    end
end
function getobjval(x,C)
    return [dot(x,C[1,:]),dot(x,C[2,:]),dot(x,C[3,:])]
end
function dominated(x,P)
    st = false
    for k=1:length(P)
        if all( x .>= P[k])# && any(x > P[k])
            st=true; break
        else
            continue
        end
    end

    return st
end
function domFilter(sol,obj)
    copysol = Dict(); copyobj = Dict();
    for i=1:length(obj)
    	copysol[i] = sol[i]
        copyobj[i] = obj[i]
    end
    for i=1:length(obj)-1
        for j=i+1:length(obj)
            if all(obj[i] .>= obj[j]) == true #dominated by PF[j]
                copyobj[i]=nothing; copysol[i]=nothing;break
            elseif all(obj[j] .>= obj[i]) == true
                copyobj[j]=nothing; copysol[j]=nothing;
            end
        end
    end
    finalsol = filter!(a->a!=nothing, collect(values(copysol)))
    finalobj = filter!(a->a!=nothing, collect(values(copyobj)))
    return finalsol,finalobj
end
function Postpro(initX,initY,dvar,LB,newsol)
    P = union(initX,dvar[end-newsol+1:end]); Pobj = union(initY,LB[end-newsol+1:end])
    #Filter dominated solutions
    copysol = Dict(); copyobj = Dict();
    for i=1:length(Pobj)
        copysol[i] = P[i]; copyobj[i] = Pobj[i]
    end
    for i=1:length(Pobj)-1
        for j=i+1:length(Pobj)
            if all(Pobj[i] .>= Pobj[j]) == true #dominated by PF[j]
                copyobj[i]=nothing; copysol[i]=nothing; break
            elseif all(Pobj[j] .>= Pobj[i]) == true
                copyobj[j]=nothing; copysol[j]=nothing;
            end
        end
    end
    finalsol = filter!(a->a!=nothing, collect(values(copysol)))
    finalobj = filter!(a->a!=nothing, collect(values(copyobj)))
    return finalsol,finalobj
end
#############################    Feasibility Pump   ############################
function fractionalFP(dvar,n,C,TL)
    X = copy(dvar); IGPair=[]; Tabu = []; newsol=0; t0=time(); #Y = copy(LB);
    while time()-t0 < TL
        I,G = sample(1:length(X), 2, replace=false)
        x1 = X[I]; x2 = X[G];
        λ = round(rand(Float64, 1)[1]; digits=1)
        x_t = x1*λ + x2*(1-λ); SearchDone = false; iter=0; Max_iter = 10
        while [I,G]∉IGPair && iter<Max_iter && SearchDone == false
            x_r = round.(Int,x_t)
            if ( (FBcheck(x_r,n) == true) && (x_r∉X) )
                push!(X,x_r); #Y = [Y; getobjval(x_r,C)'];
                newsol+=1; SearchDone = true
            else
                if x_r ∈ Tabu
                    x_r = flipoper(Tabu,x_t,x_r)
                    if x_r==[]
                        SearchDone = true
                    else
                        if ( (FBcheck(x_r,n) == true) && (x_r∉X) )
                            push!(X,x_r); #Y = [Y; getobjval(x_r,C)'];
                            newsol+=1; SearchDone = true;
                        end
                    end
                end
                if time()-t0 >= TL
                    break
                end
                if SearchDone == false
                    push!(Tabu,x_r)
                    x_t = fbsearch(x_r)
                    if x_t == 0 #when there's no new feasible lp sol
                        SearchDone = true
                    end
                end
            end
			iter+=1
        end
        push!(IGPair,[I,G])
    end
    return X,newsol
end
####################### Running Path Relinking  ##########################
