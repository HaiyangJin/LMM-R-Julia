
## Setup
using Random, Printf
using DataFrames, DataFrameMacros, Chain, Arrow, AlgebraOfGraphics, CairoMakie
using StatsModels, MixedModels, MixedModelsMakie
AOG = AlgebraOfGraphics;

# Load data
df_cf = DataFrame(Arrow.Table("data/example_sd.arrow"))

contr = Dict(
    :Participant => Grouping(),
    :Condition => SeqDiffCoding(levels = ["CFS", "monocular"]),
    :Congruency => SeqDiffCoding(levels = ["incongruent", "congruent"]),
    :Alignment => SeqDiffCoding(levels = ["misaligned", "aligned"]),
    :answer => SeqDiffCoding(levels = ["different", "same"]),
);

# Make varaibles for reducing random effects
f_tmp = ModelFrame(@formula(resp ~ Condition * Congruency * Alignment * answer),
            df_cf, contrasts = contr);
# coefnames(f_tmp)

df_contra = @chain f_tmp begin
    modelmatrix
    DataFrame(coefnames(f_tmp))
    rename(
        Symbol("(Intercept)") => :intercept,
        Symbol("Condition: monocular") => :Cond_main,
        Symbol("Congruency: congruent") => :Cong_main,
        Symbol("Alignment: aligned") => :Ali_main,
        Symbol("answer: same") => :Ans_main,
        Symbol("Condition: monocular & Congruency: congruent") => :Cond_Cong,
        Symbol("Condition: monocular & Alignment: aligned") => :Cond_Ali,
        Symbol("Congruency: congruent & Alignment: aligned") => :Cong_Ali,
        Symbol("Condition: monocular & answer: same") => :Cond_Ans,
        Symbol("Congruency: congruent & answer: same") => :Cong_Ans,
        Symbol("Alignment: aligned & answer: same") => :Ali_Ans,
        Symbol("Condition: monocular & Congruency: congruent & Alignment: aligned") => :Cond_Cong_Ali,
        Symbol("Condition: monocular & Congruency: congruent & answer: same") => :Cond_Cong_Ans,
        Symbol("Condition: monocular & Alignment: aligned & answer: same") => :Cond_Ali_Ans,
        Symbol("Congruency: congruent & Alignment: aligned & answer: same") => :Cong_Ali_Ans,
        Symbol("Condition: monocular & Congruency: congruent & Alignment: aligned & answer: same") => :Cond_Cong_Ali_Ans
    )
    # @transform(:coefname = 
	# 	:coefname == "(Intercept)" ? "intercept" :
	# 	:coefname == "Condition: monocular" ? "Cond_main" :
	# 	:coefname == "Congruency: congruent" ? "Cong_main" :
	# 	:coefname == "Alignment: aligned" ? "Ali_main" :
    #   :coefname == "answer: same" ? "Ans_main" :
	# 	:coefname == "Condition: monocular & Congruency: congruent" ? "Cond_Cong" :
	# 	:coefname == "Condition: monocular & Alignment: aligned" ? "Cond_Ali" :
	# 	:coefname == "Congruency: congruent & Alignment: aligned" ? "Cong_Ali" :
    #   :coefname == "Condition: monocular & answer: same" ? " Cond_Ans" :
    #   :coefname == "Congruency: congruent & answer: same" ? "Cong_Ans" :
    #   :coefname == "Alignment: aligned & answer: same" ? "Ali_Ans" :
	# 	:coefname == "Condition: monocular & Congruency: congruent & Alignment: aligned" ? "Cond_Cong_Ali" :
    #   :coefname == "Condition: monocular & Congruency: congruent & answer: same" ? "Cond_Cong_Ans" :
    #   :coefname == "Condition: monocular & Alignment: aligned & answer: same" ? "Cond_Ali_Ans" :
    #   :coefname == "Congruency: congruent & Alignment: aligned & answer: same" ? "Cong_Ali_Ans" :
    #   :coefname == "Condition: monocular & Congruency: congruent & Alignment: aligned & answer: same" ? "Cond_Cong_Ali_Ans" :
	# "unknown")
    # DataFrame([:intercept,
    #     :Cond_main, :Cong_main, :Ali_main, :Ans_main, 
    #     :Cond_Cong, :Cond_Ali, :Cond_Ans, :Cong_Ali, :Cong_Ans, :Ali_Ans,
    #     :Cond_Cong_Ali, :Cond_Cong_Ans, :Cond_Ali_Ans, :Cong_Ali_Ans,
    #     :Cond_Cong_Ali_Ans])
    select(Not(:intercept))
end

# compare contrast in Julia and R
df_r_con = select(df_cf, r"_")  # contrast from R
sum(names(df_r_con) .== names(df_contra))  # same column names
df_r_con == df_contra # values are the same

# concatenate 
df_lmm = @chain df_cf begin
    select!(Not(r"_"))
    hcat(df_contra)
end

## Behavioral responses
# with signal detection models -- GLMM with bernoulli family and probit link

# (an equivalent) ZCP model
f_resp_zcp_ = @formula(resp ~ Condition * Congruency * Alignment * answer +
        zerocorr(Condition * Congruency * Alignment * answer | Participant));
glmm_resp_zcp_ = fit(
    MixedModel, 
    f_resp_zcp_, 
    df_lmm, 
    Bernoulli(), 
    ProbitLink(), 
    contrasts = contr)

issingular(glmm_resp_zcp_)

# ZCP model
f_resp_zcp = @formula(resp ~ Condition * Congruency * Alignment * answer +
		zerocorr(Cond_main + Cong_main + Ali_main + Ans_main +
			Cond_Cong + Cond_Ali + Cond_Ans + Cong_Ali + Cong_Ans + Ali_Ans +
			Cond_Cong_Ali + Cond_Cong_Ans + Cond_Ali_Ans + Cong_Ali_Ans +
			Cond_Cong_Ali_Ans | Participant));
	
glmm_resp_zcp = fit(
	MixedModel, 
	f_resp_zcp, 
	df_lmm,
	Bernoulli(),
	ProbitLink(),
    contrasts = contr)

issingular(glmm_resp_zcp)

glmm_resp_zcp.rePCA
VarCorr(glmm_resp_zcp)

# Reduced model
f_resp_rdc = @formula(resp ~ Condition * Congruency * Alignment * answer +
		zerocorr(Cond_main + Cong_main + Ali_main + Ans_main +
			Cond_Ali + Cong_Ans + Ali_Ans + Cond_Ans + # Cond_Cong + Cong_Ali + 
			Cong_Ali_Ans + Cond_Cong_Ali  # + Cond_Ali_Ans + Cond_Cong_Ans + 
			| Participant)); # Cond_Cong_Ali_Ans
	
glmm_resp_rdc = fit(
    MixedModel, 
    f_resp_rdc, 
    df_lmm,
    Bernoulli(),
    ProbitLink(),
    contrasts = contr)

issingular(glmm_resp_rdc)

# Extended model
f_resp_etd = @formula(resp ~ Condition * Congruency * Alignment * answer +
		(Cond_main + Cong_main + Ali_main + Ans_main +
			Cond_Ali + Cong_Ans + Ali_Ans + Cond_Ans + # Cond_Cong + Cong_Ali + 
			Cong_Ali_Ans + Cond_Cong_Ali  # + Cond_Ali_Ans + Cond_Cong_Ans + 
			| Participant)); # Cond_Cong_Ali_Ans
	
glmm_resp_etd = fit(
    MixedModel, 
    f_resp_etd, 
    df_lmm,
    Bernoulli(),
    ProbitLink(),
    contrasts = contr)

issingular(glmm_resp_etd)

glmm_resp_etd.rePCA
VarCorr(glmm_resp_etd)

# Extended model 2
f_resp_etd2 = @formula(resp ~ Condition * Congruency * Alignment * answer +
		(Ans_main + # # Cong_main + Ali_main + Cond_main + 
			Cond_Ans + # Cond_Cong + Cong_Ali + # Cong_Ans + Ali_Ans + Cond_Ali + 
			Cong_Ali_Ans + Cond_Cong_Ali  # + Cond_Ali_Ans + Cond_Cong_Ans + 
			| Participant)); # Cond_Cong_Ali_Ans
	
glmm_resp_etd2 = fit(
    MixedModel, 
    f_resp_etd2, 
    df_lmm,
    Bernoulli(),
    ProbitLink(),
    contrasts = contr)

issingular(glmm_resp_etd2)

glmm_resp_etd2.rePCA
VarCorr(glmm_resp_etd2)

# Extended model 3
f_resp_etd3 = @formula(resp ~ Condition * Congruency * Alignment * answer +
		(0 + Ans_main + # # Cong_main + Ali_main + Cond_main + 
			Cond_Ans + # Cond_Cong + Cong_Ali + # Cong_Ans + Ali_Ans + Cond_Ali + 
			Cong_Ali_Ans   # + Cond_Ali_Ans + Cond_Cong_Ans + + Cond_Cong_Ali
			| Participant)); # Cond_Cong_Ali_Ans
	
glmm_resp_etd3 = fit(
    MixedModel, 
    f_resp_etd3, 
    df_lmm,
    Bernoulli(),
    ProbitLink(),
    contrasts = contr)

issingular(glmm_resp_etd3)

# Optimal model
MixedModels.likelihoodratiotest(glmm_resp_etd3, glmm_resp_rdc)

mods = [glmm_resp_etd3, glmm_resp_rdc];
DataFrame(
    model = [:glmm_resp_etd3, :glmm_resp_rdc],
    npar = dof.(mods),
    deviance = deviance.(mods),
    AIC = aic.(mods),
    BIC = bic.(mods),
    AICc = aicc.(mods),
)

# glmm_resp_rdc is used as the optimal model
glmm_resp_opt = glmm_resp_rdc 

# bootstrap (save multiple files)
# this loops took more than 20 hours
# for i in 8:10
#     bs_resp = parametricbootstrap(MersenneTwister(42+i), 400, glmm_resp_opt)
#     Arrow.write(@sprintf("output_jl/bs_resp_%d.arrow", i), DataFrame(bs_resp.allpars))
# end

# the rest analysis for responses is done in lmm_CF_jl.Rmd

# design = Dict(
#     :answer => ["same", "different"],
# 		:Alignment => ["aligned", "misaligned"],
# 		:Congruency => ["congruent", "incongruent"],
# 		:Condition => ["monocular", "CFS"],
# )
# preds = effects(design, glmm_resp_opt)


#################################
## Response times
df_rt = @subset(df_lmm, :isCorrect == 1)

# (an equivalent) ZCP model
f_rt_zcp_ = @formula(logRT ~ Condition * Congruency * Alignment +
		zerocorr(Condition * Congruency * Alignment | Participant));
	
lmm_rt_zcp_ = fit(
    MixedModel, 
    f_rt_zcp_, 
    df_rt,
    contrasts = contr)

issingular(lmm_rt_zcp_)

# ZCP model
f_rt_zcp = @formula(logRT ~ Condition * Congruency * Alignment +
		zerocorr(Cond_main + Cong_main + Ali_main + 
			Cond_Cong + Cond_Ali + Cong_Ali + 
			Cond_Cong_Ali | Participant));
	
lmm_rt_zcp = fit(
    MixedModel, 
    f_rt_zcp, 
    df_rt,
    contrasts = contr)

issingular(lmm_rt_zcp)

lmm_rt_zcp.rePCA
VarCorr(lmm_rt_zcp)

# Reduced model
f_rt_rdc = @formula(logRT ~ Condition * Congruency * Alignment +
		zerocorr(Cond_main + # Cong_main + Ali_main + 
        Cond_Ali # + Cond_Cong + Cong_Ali + 
			| Participant)); # Cond_Cong_Ali
	
lmm_rt_rdc = fit(
    MixedModel, 
    f_rt_rdc, 
    df_rt,
    contrasts = contr)

issingular(lmm_rt_rdc)

# Extended model
f_rt_etd = @formula(logRT ~ Condition * Congruency * Alignment +
        (Cond_main + # Cong_main + Ali_main + 
        Cond_Ali # + Cond_Cong + Cong_Ali + 
            | Participant)); # Cond_Cong_Ali

lmm_rt_etd = fit(
    MixedModel, 
    f_rt_etd, 
    df_rt,
    contrasts = contr)

issingular(lmm_rt_etd)

lmm_rt_etd.rePCA
VarCorr(lmm_rt_etd)

# Extended model 2
f_rt_etd2 = @formula(logRT ~ Condition * Congruency * Alignment +
        (Cond_main # Cong_main + Ali_main + 
             # + Cond_Cong + Cond_Ali + Cong_Ali
            | Participant)); # Cond_Cong_Ali

lmm_rt_etd2 = fit(
    MixedModel, 
    f_rt_etd2, 
    df_rt,
    contrasts = contr)

issingular(lmm_rt_etd2)

# Optimal model
MixedModels.likelihoodratiotest(lmm_rt_etd2, lmm_rt_rdc)

mods = [lmm_rt_etd, lmm_rt_rdc];
DataFrame(
    model = [:lmm_rt_etd, :lmm_rt_rdc],
    npar = dof.(mods),
    deviance = deviance.(mods),
    AIC = aic.(mods),
    BIC = bic.(mods),
    AICc = aicc.(mods),
)

# the reduced model is used as the optimal model
lmm_rt_opt = lmm_rt_rdc

shrinkageplot(lmm_rt_opt)

# bootstrap
bs_rt = parametricbootstrap(MersenneTwister(42), 4000, lmm_rt_opt);
Arrow.write("output_jl/bs_rt.arrow", DataFrame(bs_rt.allpars))

bs_int = DataFrame(shortestcovint(bs_rt))

# the rest analysis for RT is done in lmm_CF_jl.Rmd

# design_rt = Dict(
# 		:Alignment => ["aligned", "misaligned"],
# 		:Congruency => ["congruent", "incongruent"],
# 		:Condition => ["monocular", "CFS"],
# )

# preds = effects(design_rt, lmm_rt_opt)
