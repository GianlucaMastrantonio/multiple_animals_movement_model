
#
# InformationalCriteria(ModelOutput::Dict)
#
#     PosteriorSamples    = ModelOutput["PosteriorSamples"]
#     Likelihood          = ModelOutput["Likelihood"]
#     Clusterization      = ModelOutput["Clusterization"]
#     compute_InformationalCriteria(PosteriorSamples,Likelihood,Clusterization)
# end
#
# function compute_InformationalCriteria(MCMCout::Dict,Likelihood::Likelihood_OUmodel, Clusterization::Clusterization_HMM)
#
#     PosteriorSamples    = ModelOUT_TOT["PosteriorSamples"]
#     # Likelihood          = ModelOUT_TOT["Likelihood"]
#     # Clusterization      =ModelOUT_TOT["Clusterization"]
#
#     #
#     kmax        = Likelihood.kmax
#     nc          = Likelihood.data.ncol
#     nanim       = Likelihood.data.nanimals
#     nt          = Likelihood.data.nt
#
#     nsamp_mcmc  = size(PosteriorSamples["mu0"],1)
#
#     PostLikelihood = zeros(Float64,nsamp_mcmc)
#     AppMatrixL     = zeros(Float64,nt,kmax)
#     for i_mcmc in 1:nsamp_mcmc
#
#             for iobs in 2:nt
#
#                 k   = PosteriorSamples["zeta"][i_mcmc,iobs-i]
#                 yt  =
#
#             end
#     end
#     LogPdfMatrix = zeros(Float64,nt,kmax)
#     @inbounds for i  in 1:(nt-1)
#         for k in 1:kmax
#              LogPdfMatrix[i+1,k] = BayesianAnimalMovementModels.compute_logpdf(Likelihood, Int16(k), Int16(i+1))
#         end
#     end
#
#     for k1 in 1:kmax
#         for k2 in 1:kmax
#
#         end
#     end
#
# end
