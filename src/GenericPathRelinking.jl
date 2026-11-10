# Generic path-relinking routines.
# Extracted from phd/codes/FPGPRfunctions.jl; see docs/source_audit.md.

function createNB(SI,C,dif,exploredSI)
    neibour = []; neiobj = [];
    for i in dif
        cpSI = copy(SI)
        if cpSI[i] == round(cpSI[i]) #if vari is int value
            if cpSI[i] == 1
                cpSI[i] = 0
            else
                cpSI[i] = 1
            end
            push!(neibour, cpSI);
            push!(neiobj, getobjval(cpSI,C))
        else #if vari is fractional
            cpSI[i] = 1; push!(neibour, cpSI); push!( neiobj, getobjval(cpSI,C) )
            cpSI = copy(SI)
            cpSI[i] = 0; push!(neibour, cpSI); push!( neiobj, getobjval(cpSI,C) )
        end
    end
    idx = findall(i-> neibour[i] in exploredSI, 1:length(neibour))
    neibour = setdiff(neibour,exploredSI)
    neiobj = neiobj[setdiff(1:end, idx),:]
    return neibour,neiobj
end
function nextSI(neibour,neiobj,C,SI)
    SIobj = getobjval(SI,C)
    for i=1:length(neiobj)
        if length(neiobj) == 1  #if there is one candiate sol
            return neibour[1]#,neiobj[1]
        else length(neiobj) > 1 # if there are multiple candiates, check the improved ratio
            ratiotb = zeros(length(neiobj),3)
            for i=1:length(neiobj)
                ratiotb[i,:] = neiobj[i]./SIobj
            end
            ranktb = zeros(length(neiobj),3)
            for i=1:3
                ranktb[:,i] = tiedrank(ratiotb[:,i])
            end
            ranksum = [sum(ranktb[i,:]) for i=1:length(neiobj)]
            mostimp = findall(x-> x == maximum(ranksum), ranksum)
            k = rand(mostimp)[1]
            return neibour[k]#, neiobj[k]
        end
    end
end

function GPR(candX,candY,C,n,TL)
    IGPair=[]; exploredSI = []; newsol=0; t0=time();
    while time()-t0 < TL
	    I,G = sample(1:length(candX), 2, replace=false)
        SI = candX[I]; SG = candX[G]; infeas=0; #Maxiter = length(candX)*10;  #if the algorithm finds infeasible solutions in a row, move to the next iteration
        while all.(SI != SG) && [I,G]∉IGPair && iter<20 && (time()-t0<TL)
            dif = findall(i-> SI[i]!=SG[i], 1:n)
            neibour,neiobj = createNB(SI,C,dif,exploredSI)
            if ( (length(neiobj)==0) || (time()-t0 >= TL) )
                break
            else
                for l=1:length(neiobj)
                    if ( (FBcheck(neibour[l],n)==true) && neibour[l]∉candX ) #dominance check (dominated(neiobj[l],LB)==false)
                        push!(candX, neibour[l]); push!(candY, neiobj[l]);
                        newsol+=1;
                    end
                end
            end
            SI = nextSI(neibour,neiobj,C,SI)
            push!(exploredSI,SI);
            iter+=1
        end
        push!(IGPair,[I,G])
    end
    return candX,candY,newsol
end


function FPGPR(dvar,LB,n,C,TL)
    candX = copy(dvar); candY = copy(LB); IGPair=[];  newsol=0;
	gMax = 20; giter=0;	exploredSI = []; Tabu = []; t0=time();
    while time()-t0 < TL #&& newsol<=length(candY)*0.5
        I,G = sample(1:length(candY), 2, replace=false) #BASIC VER
		SI = candX[I]; SG = candX[G];
        while all.(SI != SG) && [I,G]∉IGPair && giter<gMax
			dif = findall(i-> SI[i]!=SG[i], 1:n)
			neibour,neiobj = createNB(SI,C,dif,exploredSI)
			if length(neibour)==0
				break
			else
				for l=1:length(neibour)
					SearchDone = false; fiter=0; fMax = 10;
					while fiter<fMax && SearchDone == false && time()-t0 <TL
						x_r = round.(Int,neibour[l])
						fx = getobjval(x_r,C)
			            if ( (FBcheck(x_r,n) == true) && (fx ∉ candY) ) #checking feasibility and dominance   #z[1]<=fmax[1]-1 && z[2]<=fmax[2]-1 && z[3]<=fmax[3]-1)
			                push!(candX,x_r); push!(candY,fx)#Y = [Y; getobjval(x_r,C)'];
			                newsol+=1; SearchDone = true
			            else
			                if x_r ∈ Tabu
			                    x_r = flipoper(Tabu,neibour[l],x_r)
			                    if x_r==[]
			                        SearchDone = true
			                    else
			                        fx = getobjval(x_r,C)
			                        if ( (FBcheck(x_r,n) == true) && (fx ∉ candY) )  #z[1]<=fmax[1]-1 && z[2]<=fmax[2]-1 && z[3]<=fmax[3]-1 )
			                            push!(candX,x_r); push!(candY,fx)#Y = [Y; getobjval(x_r,C)'];
			                            newsol+=1; SearchDone = true;# push!(sucess,[I,G])
			                        end
			                    end
			                end
			                if SearchDone == false
			                    push!(Tabu,x_r)
			                    x_t = fbsearch(x_r)
			                    if x_t == 0 #when there's no new feasible lp sol
			                        SearchDone = true
			                    end
			                end
			            end
						fiter+=1
        			end
				end
			end
			SI = nextSI(neibour,neiobj,C,SI)
			push!(exploredSI,SI);
			giter+=1
		end
		push!(IGPair,[I,G])
    end

	while time()-t0 < TL
		I,G = sample(1:length(candY), 2, replace=false) #BASIC VER
		SI = candX[I]; SG = candX[G]; giter=0

        while all.(SI != SG) && [I,G]∉IGPair && giter<gMax
			dif = findall(i-> SI[i]!=SG[i], 1:n)
			neibour,neiobj = createNB(SI,C,dif,exploredSI)
			if length(neibour)==0
				break
			else
				for l=1:length(neibour)
					if FBcheck(neibour[l],n) && neibour[l]∉ candX
						push!(candX, neibour[l]); candY = [candY; neiobj[l]']; #push!(candY, neiobj[l]);
						newsol+=1; gn+=1;
					end
				end
			end
			SI = nextSI(neibour,neiobj,C,SI)
			push!(exploredSI,SI);
			giter+=1
		end
		push!(IGPair,[I,G])
	end

    return candX,candY,newsol
end
