This repository stores different approaches to analyze data for the composite face task. 

<br>

# Brief introduction
The data was from a subset of [one previous study](https://psyarxiv.com/j8g6z/). 

Independent variabls:
+ Condition: `monocular` (**O**) vs. `CFS` (**F**)  
+ Congruency: `congruent` (**C**) vs. `incongruent` (**I**)  
+ Alignment: `aligned` (**A**) vs. `misaligned` (**M**)  
+ answer: `same` (**S**) vs. `different` (**D**) (used as `signal` and `noise` in Signal Detection Theory; combined when analyzing response times)

Dependent variables:  
+ Sensitivity d' (logistic regression with `probit` link)
+ Correct response times (`lognormal` transformation)

As discussed [earlier](https://psyarxiv.com/yhmzg/), we may claim observing the composite face effect in a particular condition (e.g., `CFS` or `monocular`) only when (1) the performance for aligned faces is better in the congruent relative to the incongruent condition (i.e., `congruent_aligned` > `incongruent_aligned`) and (2) the increased performance in congruent relative to incongruent condition is larger for aligned compared to misaligned faces (i.e., (`congruent_aligned` - `incongruent_aligned`) > (`congruent_misaligned` - `incongruent_misaligned`)). 

Moreover, for the third question, we also need to examine whether the composite effect in the `monocular` condition is larger than that in the `CFS` condition (i.e., [(`monocular_congruent_aligned` - `monocular_incongruent_aligned`) - (`monocular_congruent_misaligned` - `monocular_incongruent_misaligned`)] > [(`CFS_congruent_aligned` - `CFS_incongruent_aligned`) - (`CFS_congruent_misaligned` - `CFS_incongruent_misaligned`)]). 

<br> 
In summary, the effects of interst in this study are:

1. `congruent` vs. `incongruent` for aligned faces in `monocular` and `CFS` condition, respectively (E1). 
2. The two-way interaction between `Congruency` and `Alignment` for the monocular and CFS condition separately (E2). 
3. The three-way interaction between `Condition`, `Congruency` and `Alignment` (E3).

<br>

# Analyses
There are three RMarkdown files that analyzed the composite face task with (generalized) linear mixed-effects models (using `library(lme4)`):

- lmm_CF_emmeans.rmd
- lmm_CF_nested.rmd
- lmm_CF_priori.rmd

## emmeans
`lmm_CF_emmeans.rmd` creates one lmm model and then uses `library(emmeans)` to examine the effects of interests.

## nested
`lmm_CF_nested.rmd` creates three models with different nested effects (e.g., `DV ~ A / B * C`) to examine the above three kinds of effects of interest, respectively. 

## priori
I inteded to use priori contrasts to examine all the three kinds of effects in the same model but realized that it is impossible due to "singular contrast matrix" issue. Also, I'm not sure how I should set the other contrasts to explain as much variances as possible. Still trying to figure out what I should do...

<br> 

# My questions
1. What are the potentional issues (disadvantages) with the `emmeans` approach?
2. The `nested` approach seems to be a bit cumbersome. Is there any way to imporve this? Except for this, is there any other  disadvantages for this approach?
3. Is it possible to use priori contrasts to test all the three effects of interests? If not, what is the (more) optimal approach (potentially this is relating to the second question)?
4. Results of `emmeans` vs. `nested`: I notice that the results for response times (lognormal) obtained from `nested` match those from `emmeans` (at least for the effects of interest). Is this a coincidence? (But the results for behavioral responses do not match; logistic regression with `probit` link) 