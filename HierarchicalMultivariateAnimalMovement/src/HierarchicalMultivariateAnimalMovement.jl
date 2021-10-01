module HierarchicalMultivariateAnimalMovement

##### Packages
using Distributions, Random
using LinearAlgebra, PDMats
using Impute
using SpecialFunctions
#import: Distributions:
##### Include
include(joinpath("checks.jl"))
include(joinpath("functions.jl"))
include(joinpath("settings/missing.jl"))

include(joinpath("settings/vectorparameters.jl"))
include(joinpath("settings/posdefmats.jl"))
include(joinpath("settings/zeta.jl"))
include(joinpath("settings/pi.jl"))
#
# include(joinpath("settings/mu0.jl"))
# include(joinpath("settings/psi.jl"))
# include(joinpath("settings/sigma.jl"))
# include(joinpath("settings/pi.jl"))

include(joinpath("settings/data.jl"))
include(joinpath("settings/options.jl"))
include(joinpath("settings/utils.jl"))


include(joinpath("mcmc/rhodp.jl"))
include(joinpath("mcmc/gamma.jl"))
include(joinpath("mcmc/ak.jl"))
include(joinpath("mcmc/beta.jl"))
include(joinpath("mcmc/matrixm.jl"))
include(joinpath("mcmc/mu0.jl"))
include(joinpath("mcmc/muC.jl"))
include(joinpath("mcmc/sigma.jl"))
include(joinpath("mcmc/pi.jl"))
include(joinpath("mcmc/psi.jl"))
include(joinpath("mcmc/rho.jl"))

include(joinpath("mcmc/mu0_lv2.jl"))
include(joinpath("mcmc/muC_lv2.jl"))
include(joinpath("mcmc/sigma_lv2.jl"))
include(joinpath("mcmc/psi_lv2.jl"))
include(joinpath("mcmc/rho_lv2.jl"))
include(joinpath("mcmc/prob_lv2.jl"))


include(joinpath("mcmc/data.jl"))
include(joinpath("mcmc/zeta.jl"))
include(joinpath("mcmc/missing.jl"))
include(joinpath("mcmcalgorithms.jl"))
include(joinpath("informationalcriteria.jl"))
##### Functions
export
    OptionsLikelihood,
    OptionsClusterization,
    OptionsMCMC,
    MCMCalgorithm,
    MCMCalgorithm_joint,
    MCMCalgorithm_singleHDP,
    InformationalCriteria,
    MCMCutils,
    OptionHierarchicalParameters


end # module BayesianAnimalMovementModels
