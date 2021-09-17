
## Setup
using DataFrames, DataFrameMacros, Chain, Arrow, AlgebraOfGraphics, CairoMakie
using StatsModels, MixedModels, MixedModelsMakie, Random
AOG = AlgebraOfGraphics;

# Load data
df_cf = DataFrame(Arrow.Table("data/example_sd.arrow"))
contr = Dict(
    :Participant => Grouping(),
    :Condition => SeqDiffCoding(levels = ["CFS", "monocular"]),
    :Congruency => SeqDiffCoding(levels = ["incongruent", "congruent"]),
    :Aligment => SeqDiffCoding(levels = ["misaligned", "aligned"]),
    :answer => SeqDiffCoding(levels = ["different", "same"]),
);

# Make varaibles for reducing random effects
df_contra = @chain df_cf begin
    ModelFrame(@formula(resp ~ Condition * Congruency * Alignment * answer),
               _,
               contrasts = contr,)
    modelmatrix
    DataFrame([:intercept,
        :cond_main, :cong_main, :ali_main, :ans_main, 
        :cond_cong, :cond_ali, :cond_ans, :cong_ali, :cong_ans, :ali_ans,
        :cond_cong_ali, :cond_cong_ans, :cond_ali_ans, :cong_ali_ans,
        :cond_cong_ali_ans])
    hcat(df_cf, _)
    @transform(:logRT = log(:RT))
end

## Behavioral responses
# with signal detection models -- GLMM with bernoulli family and probit link

# (an equivalent) ZCP model
f_resp_zcp_ = @formula(resp ~ Condition * Congruency * Alignment * answer +
        zerocorr(Condition * Congruency * Alignment * answer | Participant));

glmm_resp_zcp_ = fit(
    MixedModel, 
    f_resp_zcp_, 
    df_contra, 
    Bernoulli(), 
    ProbitLink(), 
    contrasts = contr)

issingular(glmm_resp_zcp_)

# ZCP model
f_resp_zcp = @formula(resp ~ Condition * Congruency * Alignment * answer +
		zerocorr(cond_main + cong_main + ali_main + ans_main +
			cond_cong + cond_ali + cond_ans + cong_ali + cong_ans + ali_ans +
			cond_cong_ali + cond_cong_ans + cond_ali_ans + cong_ali_ans +
			cond_cong_ali_ans | Participant));
	
glmm_resp_zcp = fit(
	MixedModel, 
	f_resp_zcp, 
	df_contra,
	Bernoulli(),
	ProbitLink(),
    contrasts = contr)

issingular(glmm_resp_zcp)

glmm_resp_zcp.rePCA
VarCorr(glmm_resp_zcp)

# Reduced model
f_resp_rdc = @formula(resp ~ Condition * Congruency * Alignment * answer +
		zerocorr(cond_main + cong_main + ali_main + ans_main +
			cond_ali + cong_ali + cong_ans + # ali_ans + cond_ans + cond_cong + 
			cong_ali_ans # + cond_ali_ans + cond_cong_ans + cond_cong_ali + 
			| Participant)); # cond_cong_ali_ans
	
glmm_resp_rdc = fit(
    MixedModel, 
    f_resp_rdc, 
    df_contra,
    Bernoulli(),
    ProbitLink(),
    contrasts = contr)

issingular(glmm_resp_rdc)

# Extended model
f_resp_etd = @formula(resp ~ Condition * Congruency * Alignment * answer +
        (cond_main + cong_main + ali_main + ans_main +
            cond_ali + cong_ali + cong_ans + # ali_ans + cond_ans + cond_cong + 
            cong_ali_ans # + cond_ali_ans + cond_cong_ans + cond_cong_ali + 
            | Participant)); # cond_cong_ali_ans
	
glmm_resp_etd = fit(
    MixedModel, 
    f_resp_etd, 
    df_contra,
    Bernoulli(),
    ProbitLink(),
    contrasts = contr)

issingular(glmm_resp_etd)

# Optimal model
MixedModels.likelihoodratiotest(glmm_resp_etd, glmm_resp_rdc)

mods = [glmm_resp_etd, glmm_resp_rdc];
DataFrame(
    model = [:glmm_resp_etd, :glmm_resp_rdc],
    npar = dof.(mods),
    deviance = deviance.(mods),
    AIC = aic.(mods),
    BIC = bic.(mods),
    AICc = aicc.(mods),
)

# glmm_resp_rdc is used as the optimal model
glmm_resp_opt = glmm_resp_rdc 

# bootstrap (to be finished)
bs_resp1 = parametricbootstrap(MersenneTwister(42), 2, glmm_resp_opt);
# bs_resp1.allpars

#################################
## Response times
df_rt = @subset(df_contra, :isCorrect == 1)

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
		zerocorr(cond_main + cong_main + ali_main + 
			cond_cong + cond_ali + cong_ali + 
			cond_cong_ali | Participant));
	
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
		zerocorr(cond_main + # cong_main + ali_main + 
			cond_ali + cong_ali # + cond_cong + 
			| Participant)); # cond_cong_ali
	
lmm_rt_rdc = fit(
    MixedModel, 
    f_rt_rdc, 
    df_rt,
    contrasts = contr)

issingular(lmm_rt_rdc)

# Extended model
f_rt_etd = @formula(logRT ~ Condition * Congruency * Alignment +
        (cond_main + # cong_main + ali_main + 
            cond_ali + cong_ali # + cond_cong + 
            | Participant)); # cond_cong_ali

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
        (cond_main + # cong_main + ali_main + 
            cong_ali # + cond_cong + cond_ali + 
            | Participant)); # cond_cong_ali

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
bs_rt = parametricbootstrap(MersenneTwister(42), 2000, lmm_rt_opt);

bs_rt.allpars