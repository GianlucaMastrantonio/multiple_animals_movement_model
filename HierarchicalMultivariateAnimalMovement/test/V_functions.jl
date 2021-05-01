
IndexA   = Int32[2;4]
#IndexA   = Int32[1;2;3;4;5;6]
CovInv = rand(InverseWishart(20, Matrix{Float64}(I,6,6)  ) )
Mean    = [1.0;2.0;3.0;-0.4;0.1;-10.9]
Obs     = rand(6)

CondMean, CondVar, CondVarInv = BayesianAnimalMovementModels.compute_CondMeanAndVariance_MvNormal(IndexA,Mean, CovInv,Obs)

@testset "Conditional Mean and Variance" begin
    IndexB = [1;3;5;6]
    Cov = inv(CovInv)
    CondMean2 = Mean[IndexA]+Cov[IndexA,IndexB]*inv(Cov[IndexB,IndexB])*(Obs[IndexB]-Mean[IndexB]  )
    CondVar2 = Cov[IndexA,IndexA]-Cov[IndexA,IndexB]*inv(Cov[IndexB,IndexB])*Cov[IndexB,IndexA]
    @test CondMean ≈ CondMean2
    @test CondVar  ≈ CondVar2
end
