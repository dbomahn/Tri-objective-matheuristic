using DataStructures,DelimitedFiles,DataFrames,JuMP,CPLEX,GLPK,LinearAlgebra,CPUTime,CSV,Statistics,StatsBase,MathProgBase,MathOptInterface,MathOptFormat,JLD2,SparseArrays
const MPB = MathProgBase; const MOF = MathOptFormat; const MOI = MathOptInterface;
mutable struct Data
    mpsfile::String; vlpfile::String
    m::Int; n::Int; C::Array{}; B::Array{}; RHS::Dict{}; signs::Array{};
    function Data(mpsfile::String, vlpfile::String)
        mpsmodel = CPLEX.CplexMathProgModel(); MPB.loadproblem!(mpsmodel,mpsfile)
        B = MPB.getconstrmatrix(mpsmodel); m,n=size(B)
        lb = MPB.getconstrLB(mpsmodel);ub = MPB.getconstrUB(mpsmodel);
        RHS = Dict()
        for i=1:m
            if ub[i]==Inf
                RHS[i] = lb[i]
            else
                RHS[i] = ub[i]
            end
        end
        signs = []
        for i=1:m
            if ub[i] == Inf
                push!(signs,"l")
            elseif lb[i] == -Inf
                push!(signs,"u")
            else
                push!(signs, "s")
            end
        end
        vlp = readdlm(vlpfile); os = findall(i->vlp[i,1]=="o", 1:length(vlp[:,1]))
        inio = os[1]; endo = os[end]; Pmtx = vlp[inio:endo,2:4]; C = zeros(3,n)
        for k=1:length(Pmtx[:,1]) #(len,id) enumerate
            x = Pmtx[k,:][1]; y = Pmtx[k,:][2];
            C[x,y] = Pmtx[k,:][3]
        end
        new(mpsfile,vlpfile,m,n,C,B,RHS,signs)
    end
end
mutable struct Valu
    x::String; y::String; dvar::Array{}; LB::Array{}; LBmtx::Array{};
    function Valu(x,y)
        JLD2.@load x dv;
        dv0 = Array(dv);
        dv1 = round.(dv0; digits=4);
        objs = round.(digits=4, readdlm(y));
        ind = findall(i-> 0 in objs[i,:]  , 1:size(objs)[1]);
        dv2 = dv1[setdiff(1:end, ind), :];
        LBmtx = objs[setdiff(1:end, ind), 2:end];
        dvar = [Vector(dv2[i, :]) for i = 1:size(LBmtx)[1]];
        LB = [Vector(LBmtx[i, :]) for i = 1:size(LBmtx)[1]];
        new(x,y,dvar,LB,LBmtx)
    end
end
