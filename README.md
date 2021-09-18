Four approaches are used to analyze the same data set (although not all of them work). 

<br>

# Brief introduction
The data are a subset of <a href="https://psyarxiv.com/j8g6z/" target="_blank">one previous study</a> where data were analyzed with repeated-measures ANOVA. The independent variables are:

Independent variables:
+ Condition: `monocular` (**O**) vs. `CFS` (**F**)  
+ Congruency: `congruent` (**C**) vs. `incongruent` (**I**)  
+ Alignment: `aligned` (**A**) vs. `misaligned` (**M**)  
+ answer: `same` (**S**) vs. `different` (**D**) (used as `signal` and `noise` in Signal Detection Theory; combined when analyzing response times)

Dependent variables:  
+ Sensitivity d' (logistic regression with `probit` link)
+ Correct response times (`lognormal` transformation)

As discussed <a href="https://psyarxiv.com/yhmzg/" target="_blank">earlier</a>, we may claim observing the composite face effect in a particular condition (e.g., `CFS` or `monocular`) only when (1) the performance for aligned faces is better in the congruent relative to the incongruent condition (i.e., `congruent_aligned` > `incongruent_aligned`) and (2) the increased performance in congruent relative to incongruent condition is larger for aligned compared to misaligned faces (i.e., (`congruent_aligned` - `incongruent_aligned`) > (`congruent_misaligned` - `incongruent_misaligned`)). 

Moreover, for the third question, we also need to examine whether the composite effect in the `monocular` condition is larger than that in the `CFS` condition (i.e., [(`monocular_congruent_aligned` - `monocular_incongruent_aligned`) - (`monocular_congruent_misaligned` - `monocular_incongruent_misaligned`)] > [(`CFS_congruent_aligned` - `CFS_incongruent_aligned`) - (`CFS_congruent_misaligned` - `CFS_incongruent_misaligned`)]). 

<br> 
In summary, the effects of interest in this study are:

1. `congruent` vs. `incongruent` for aligned faces in `monocular` and `CFS` condition, respectively (E1). 
2. The two-way interaction between `Congruency` and `Alignment` for the monocular and CFS condition separately (E2). 
3. The three-way interaction between `Condition`, `Congruency` and `Alignment` (E3).

<br>

# Analyses
There are three RMarkdown files that analyzed the composite face task with (generalized) linear mixed-effects models (using `library(lme4)`):

- [lmm_CF_emmeans.Rmd](./lmm_CF_emmeans.Rmd) ([output](https://haiyangjin.github.io/Mixed-Model-CF/lmm_CF_emmeans.html))
- [lmm_CF_nested.Rmd](./lmm_CF_nested.Rmd) ([output](https://haiyangjin.github.io/Mixed-Model-CF/lmm_CF_nested.html))
- [lmm_CF_priori.Rmd](./lmm_CF_priori.Rmd)

The same data set was also analyzed with [MixedModels.jl](https://github.com/JuliaStats/MixedModels.jl) in Julia:

- [lmm_CF_Pluto.jl](./lmm_CF_Pluto.jl)  
- [lmm_CF.jl](./lmm_CF.jl)
    - The bootstrapping samples were further anallyzed in [lmm_CF_jl.Rmd](./lmm_CF_jl.Rmd) ([output](https://haiyangjin.github.io/Mixed-Model-CF/lmm_CF_jl.html))

## emmeans
`lmm_CF_emmeans.rmd` creates one lmm model and then uses `library(emmeans)` to examine the effects of interest.

## nested
`lmm_CF_nested.rmd` creates three models with different nested effects (e.g., `DV ~ A / B * C`) to examine the above three kinds of effects of interest, respectively. 

## priori
I intended to use priori contrasts to examine all the three kinds of effects in the same model but realize that it is impossible to do so due to "singular contrast matrix" issue. Therefore, this approach is not feasible for this study. 

## jl
Codes in both .jl files are similar and the main difference is that `lmm_CF_Pluto.jl` runs in [Pluto notebooks](https://github.com/fonsp/Pluto.jl) whereas `lmm_CF.jl` runs in Julia REPL (better in VS code). Critically, the uncertainty of the estimated parameters were calculated with the parametric bootstrap approach (`parametricbootstrap()`).
