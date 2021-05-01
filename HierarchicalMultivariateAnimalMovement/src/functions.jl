function find_const(x::Vector{Int16}, k::Int16)

    # x = inits_zeta[1]
    # x = Int16.([1; 1; 1; 2; 2; 2; 3; 1; 2; 1; 1; 2; 1])
    # k = 1
    #

    times   = findall(x->x == k, x)
    n       = Int16(0)
    Vstart  = zeros(Int16,length(times))
    Vend    = zeros(Int16,length(times))

    Vstart[1] = times[1]
    n         += one(n)

    for i in 1:(length(times)-1)

        if times[i+1] != (times[i]+1)
#            global n
            Vend[n]   = times[i]
            n         += one(n)
            Vstart[n] = times[i+1]
        end
    end
    Vend[n] = times[length(times)]

    ret         = Matrix{Int16}(undef, n,2)
    ret[:,1]    = deepcopy(Vstart[1:n])
    ret[:,2]    = deepcopy(Vend[1:n])

    return ret

end

#find_const(x,1)

function check_unique(C::Matrix{Int16})
   un  = true
   seen = Set{SubArray{Int16,1,Array{Int16,2},Tuple{Int64,Base.Slice{Base.OneTo{Int64}}},true}}()
   i = Int16(1)
   while (un == true)&& (i<=size(C,1))
       #print(i, "", un,"\n")
       @inbounds x = view(C,i,:)
       if !in(x, seen)
           push!(seen, x)
       else
           un = false
       end
       i += one(i)
   end
   un

end

# function check_unique(C::Matrix{Int16}, k::Int16, kmax::Int16,Indk::Vector{Int16})
#
#     un      = true
#     i       = Int16(1)
#     view_K  = view(C,k,:)
#
#     while (un == true)&& (i<=kmax)
#
#         if view_K == view(C,i,:)
#             if i != k
#                 un = false
#             end
#         end
#         i += one(i)
#     end
#     un
#     #sum([C[k,:] == C[iii,:] for iii in 1:kmax])==1
#     #sum([Aindex[k,:] == Aindex[iii,:] for iii in 1:kmax])==1
# end

function check_unique(C::Matrix{Int16}, k::Int16, kmax::Int16)

    un      = true
    i       = Int16(1)
    view_K  = view(C,k,:)

    while (un == true)&& (i<=kmax)

        if view_K == view(C,i,:)
            if i != k
                un = false
            end
        end
        i += one(i)
    end
    un
    #sum([C[k,:] == C[iii,:] for iii in 1:kmax])==1
    #sum([Aindex[k,:] == Aindex[iii,:] for iii in 1:kmax])==1
end

function compute_AppCondMeanAndVariance_MvNormal(IndexA::Vector{T},CovInv::Matrix{Float64} ) where{T<:Integer}

    nvar        = size(CovInv,1)
    IndexB      = T[1:nvar;]
    IndexB      = deleteat!(IndexB,IndexA)

    CovInv_AgivenB = CovInv[IndexA,IndexA]
    CovAgivenB     = inv(CovInv_AgivenB)

    CovABInvCovB      = -CovAgivenB*view(CovInv,IndexA,IndexB)

    return (CovABInvCovB, CovAgivenB ,CovInv_AgivenB )


end


function compute_CondMeanAndVariance_MvNormal(IndexA::Vector{T},Mean::Vector{Float64}, CovInv::Matrix{Float64}, Obs::Vector{Float64} ) where {T<:Integer}

    nvar        = size(CovInv,1)
    IndexB      = T[1:nvar;]
    IndexB      = deleteat!(IndexB,IndexA)

    CovABInvCovB, CovAgivenB ,CovInv_AgivenB = compute_AppCondMeanAndVariance_MvNormal(IndexA,CovInv )
    MuAgivenB      = Mean[IndexA] + CovABInvCovB*(view(Obs,IndexB)- view(Mean,IndexB) )

    return ( MuAgivenB, CovAgivenB ,CovInv_AgivenB )


end


function logsumexp(X::Vector{Float64})
    #print("INIZIO \n",X,"\n")
    alpha = -Inf
    r = 0.0
    for x = X
        if x <= alpha
            r += exp(x - alpha)
        else
            r *= exp(alpha - x)
            r += 1.0
            alpha = x
        end
    end
    #print(,r,"FINE \n")
    return log(r) + alpha
end
function sample_discretevar(rng::MersenneTwister,ProbNonNorm::Vector{Float64})
#X = ProbNonNorm
    App = exp.(ProbNonNorm .-logsumexp(ProbNonNorm))
    #return findall(rand(rng,Multinomial(1,App/sum(App) )) .==1)[1];
    #app2  = findall(rand(rng,Multinomial(1,App/sum(App) )) .==1)[1];
    return findall(cumsum(App).>=rand(rng,Uniform(0.0,sum(App))))[1]
    #return 1

end
#
# int Class_Utils::sample_DiscreteVar(double *logprob_NonNormalized, int K)
# {
#     int k;
#     double sum = 0.0;;
#     double prob[K];
#     Class_Utils::log_sum_exp(&sum, K, logprob_NonNormalized);
#
#     for(k=0;k<K;k++)
#     {
#         prob[k] = exp(logprob_NonNormalized[k]-sum);
#     }
#
#     double u = runif(0.0,1.0);
#     //REprintf("u %f ",u);
#     sum = 0.0;
#     k = -1;
#     do{
#         k++;
#         sum += prob[k];
#         //REprintf("%f ",sum );
#     }while(sum<u);
#
#     //REprintf("K=%i\n",k);
#     return(k);
# }
