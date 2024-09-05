This R script performs QTL mapping of variation in *Daphnia magna*
resistance to hindgut attachment by a *Pasteuria ramosa* isolate (P21).
The and analysis are based on Broman and Sen 2009 [“A Guide to QTL
Mapping with
R/qtl”](https://link.springer.com/book/10.1007/978-0-387-92125-9)

R/qtl version 1.42-8; 19 February 2018

### PREPARE WORK ENVIRONMENT

``` r
library(qtl)    #loads R/qtl package version 1.42-8, the package must be installed before

library(help=qtl)
help(read.cross) #Use this to read the documentation about how the qtl data are read by R 
```

### IMPORT AND LOOK AT DATA

Import QTL mapping data into R by reading the csv file. The file
contains clone names, observed phenotype (binary attachment test data;
0=resistant, 1=susceptible) and all marker regions with the genotypes.

``` r
susc <- read.cross('csv', file='QTL_SNP_map.csv', genotypes=c('AA', 'AB', 'BB'))
```

    ##  --Read the following data:
    ##   340  individuals
    ##   1324  markers
    ##   2  phenotypes
    ##  --Cross type: f2

``` r
summary(susc)
```

    ##     F2 intercross
    ## 
    ##     No. individuals:    340 
    ## 
    ##     No. phenotypes:     2 
    ##     Percent phenotyped: 100 100 
    ## 
    ##     No. chromosomes:    10 
    ##         Autosomes:      1 2 3 4 5 6 7 8 9 10 
    ## 
    ##     Total markers:      1324 
    ##     No. markers:        165 161 134 138 113 124 122 111 99 157 
    ##     Percent genotyped:  87.5 
    ##     Genotypes (%):      AA:25.9  AB:47.3  BB:26.9  not BB:0.0  not AA:0.0

``` r
plot(susc)#summary plots of data
```

![](Figs/import-1.png) - Genetic map shows where all the markers are on
the 10 chromosomes - Histogram of susceptible (1) and resistant (0)
shows slightly more susceptible

\#\#\#¶ INTERVAL MAPPING

Calculate conditional genotype probabilities (predicts genotypes at
specified locations between markers. This is required before interval
mapping). The step argument indicates the distance (in cM) between
positions at which the genotype probabilities are calculated (step=1
default), error.prob allows the probabilities to be calculated assuming
a given rate of genotyping errors.

``` r
dat <- calc.genoprob(susc, step=1, error.prob=0.01)
summary(dat)
```

    ##     F2 intercross
    ## 
    ##     No. individuals:    340 
    ## 
    ##     No. phenotypes:     2 
    ##     Percent phenotyped: 100 100 
    ## 
    ##     No. chromosomes:    10 
    ##         Autosomes:      1 2 3 4 5 6 7 8 9 10 
    ## 
    ##     Total markers:      1324 
    ##     No. markers:        165 161 134 138 113 124 122 111 99 157 
    ##     Percent genotyped:  87.5 
    ##     Genotypes (%):      AA:25.9  AB:47.3  BB:26.9  not BB:0.0  not AA:0.0

``` r
#Perform a genome scan with a single-QTL model. Because my phenotype data is binary, 
#I specify model='binary'
out.bin <- scanone(dat, pheno.col=2, model='binary')
summary(out.bin, threshold=3) #Find the QTL above LOD=3
```

    ##           chr pos  lod
    ## c3.loc157   3 157 63.5

``` r
plot(out.bin, ylab="LOD score", xlab="linkage group") #plots LOD scores across genome
```

![](Figs/intervalmapping-1.png)

``` r
plot(out.bin, chr=3) #plots LOD scores only on chromosome 3
```

![](Figs/intervalmapping-2.png)

### ESTABLISH STATISTICAL SIGNIFICANCE OF QTL

Calculate the significance thresholds. Again I use the binary model.
n.perm defines the number of permutation replicates. I reduced n.perm to
1000 because binary method is more computationally demanding than
Haley-Knott method.

``` r
operm.bin <- scanone(dat, pheno.col=2, model='binary', n.perm=1000, verbose=FALSE) 
plot(operm.bin) #plots a histogram of the permutation results
```

![](Figs/significance%20thresh-1.png)

``` r
summary(operm.bin, alpha=c(0.20,0.05)) #estimates genome-wide LOD thresholds for significance levels
```

    ## LOD thresholds (1000 permutations)
    ##      lod
    ## 20% 3.01
    ## 5%  3.76

``` r
#20% and 5%. Note that the output will be slightly different each time the permutations are performed.

plot(out.bin, ylab="LOD score", xlab="Linkage group", lwd=4, cex.lab=1.5) #plots LOD scores across genome
add.threshold(out.bin, perms=operm.bin, alpha=0.05, lty=2, lwd=2)
```

![](Figs/significance%20thresh-2.png)

``` r
#Pick out LOD peaks (max 1 per chromosome) that meet 5% significance level
#and obtain genome-scan-adjusted p-value for each LOD peak
summary(out.bin, perms=operm.bin, alpha=0.05, pvalues=TRUE) 
```

    ##           chr pos  lod pval
    ## c3.loc157   3 157 63.5    0

``` r
#Get upper confidence limit on true p-value
binom.test(0, 1000)$conf.int 
```

    ## [1] 0.000000000 0.003682084
    ## attr(,"conf.level")
    ## [1] 0.95

### ESTIMATE INTERVAL OF QTL LOCATION

``` r
lodint(out.bin, 3, 1.8)  #1.8 refers to the LOD support interval (1.8 recommended for intercrosses)
```

    ##           chr pos      lod
    ## c3.loc155   3 155 61.01275
    ## c3.loc157   3 157 63.47529
    ## c3.loc160   3 160 59.55919

``` r
bayesint(out.bin, 3, 0.95) #0.95 refers to 95% Bayes credible interval
```

    ##           chr pos      lod
    ## c3.loc156   3 156 62.83319
    ## c3.loc157   3 157 63.47529
    ## c3.loc158   3 158 63.20316

### ESTIMATE QTL EFFECTS

Create an effect plot, which plots phenotype averages for genotype
groups at an inferred QTL. This function uses the multiple imputation
method to obtain estimates of genotpye-specific phenotype averages, so
sim.geno must be performed before.

To perform multiple imputations we use sim.geno. It is similar to
calc.genoprob, but has an additional argument (n.draws which is the
number of imputations).

``` r
dat2 <- sim.geno(dat, step=1, n.draws=16, error.prob=0.001)

#find the marker closest to QTL
find.marker(dat, 3, 157)
```

    ## [1] "scaffold00288_965"

``` r
marker_eff <- effectplot(dat2, pheno.col=2, mname1="scaffold00288_965")
```

![](Figs/QTL%20effects-1.png)

``` r
marker_eff #inspect effectplot object
```

    ## $Means
    ## scaffold00288_965.AA scaffold00288_965.AB scaffold00288_965.BB 
    ##          0.000024918          0.163826936          0.909517857 
    ## 
    ## $SEs
    ## scaffold00288_965.AA scaffold00288_965.AB scaffold00288_965.BB 
    ##           0.03819401           0.02881949           0.02283731

``` r
#Now using a pseudomarker (the actual QTL position) instead of an actual marker
pseudo_eff <- effectplot(dat2, pheno.col=2, mname1="c3.loc157", ylab="Proportion", xlab="", main="Susceptibility")
```

![](Figs/QTL%20effects-2.png)

``` r
pseudo_eff #inspect effectplot object
```

    ## $Means
    ##   3@157.0.AA   3@157.0.AB   3@157.0.BB 
    ## 0.0004014068 0.1361939119 0.8972022136 
    ## 
    ## $SEs
    ## 3@157.0.AA 3@157.0.AB 3@157.0.BB 
    ## 0.03770233 0.03664313 0.02204613

From this plot and output we can tell that resistance is dominant
(because heterozygotes are (mostly) resistant. 13.6193912% +/- 3.6643128
% of AB are susceptible)

### FIT A DEFINED QTL MODEL

``` r
#create qtl object with makeqtl
qtl <- makeqtl(dat, what='prob', chr=3, pos=157)
qtl #shows summary of qtl object
```

    ##   QTL object containing genotype probabilities. 
    ## 
    ##       name chr pos n.gen
    ## Q1 3@157.0   3 157     3

``` r
plot(qtl) #plots location of QTL on genetic map
```

![](Figs/fitqtl-1.png)

``` r
#refine QTL location
rqtl <- refineqtl(dat, pheno.col=2, qtl=qtl, formula=y~Q1, method="hk", model="binary")
```

    ## pos: 157 
    ## Iteration 1 
    ##  Q1 pos: 157 -> 157
    ##     LOD increase:  0 
    ## all pos: 157 -> 157 
    ## LOD increase at this iteration:  0 
    ## overall pos: 157 -> 157 
    ## LOD increase overall:  0

``` r
#no change in position 

#use fitqtl to fit the model 
out.fq <- fitqtl(dat, pheno.col=2, model='binary', qtl=qtl)
summary(out.fq) 
```

    ## 
    ##      fitqtl summary
    ## 
    ## Method: Haley-Knott regression 
    ## Model:  binary phenotype
    ## Number of observations : 340 
    ## 
    ## Full model result
    ## ----------------------------------  
    ## Model formula: y ~ Q1 
    ## 
    ##       df     LOD     %var Pvalue(Chi2)
    ## Model  2 62.8827 57.33205            0

Results show that 0 % of the variability can be explained by the found
QTL.

### MARKERS AS COVARIATES

Because the LOD for the QTL on chr3 is so high, I will set the marker
near <a href="mailto:chr3@157" class="email">chr3@157</a> as a covariate
to help find other QTL. I will also set it as an interactive covariate
to help find QTL that may interact with it.

``` r
#find marker closest to QTL, pull out genotype data, assure none is missing
mar <- find.marker(dat, 3, 157)
g <- pull.geno(dat)[,mar]
sum(is.na(g))
```

    ## [1] 24

``` r
#[1] 26

#impute missing data at the marker and use imputed genotypes as if they were observed
g.imp <- pull.geno(fill.geno(dat))[,mar]
sum(is.na(g.imp))
```

    ## [1] 0

``` r
#because this is an intercross, create a 2-column numeric matrix encoding the genotype data:
#one column indicating one of the homozygotes and the other indicating the heterozygotes
g.matrix <- cbind(as.numeric(g.imp==1), as.numeric(g.imp==2))
```

##### 1. Perform a genome scan with the marker as an **additive** covariate:

``` r
out.covar <- scanone(dat, pheno.col=2, model='binary', addcovar=g.matrix)
#15 warnings: Jacobian matrix is singular (this usually happens in regions with 
#little evidence for QTL and can generally be ignored)

#run permutation test with QTL as covariate (omit chr 3 from analysis)
operm.covar <- scanone(dat, pheno.col=2, model='binary', addcovar=g.matrix, chr=-3, n.perm=1000)
```

    ## Permutation 20 
    ## Permutation 40 
    ## Permutation 60 
    ## Permutation 80 
    ## Permutation 100 
    ## Permutation 120 
    ## Permutation 140 
    ## Permutation 160 
    ## Permutation 180 
    ## Permutation 200 
    ## Permutation 220 
    ## Permutation 240 
    ## Permutation 260 
    ## Permutation 280 
    ## Permutation 300 
    ## Permutation 320 
    ## Permutation 340 
    ## Permutation 360 
    ## Permutation 380 
    ## Permutation 400 
    ## Permutation 420 
    ## Permutation 440 
    ## Permutation 460 
    ## Permutation 480 
    ## Permutation 500 
    ## Permutation 520 
    ## Permutation 540 
    ## Permutation 560 
    ## Permutation 580 
    ## Permutation 600 
    ## Permutation 620 
    ## Permutation 640 
    ## Permutation 660 
    ## Permutation 680 
    ## Permutation 700 
    ## Permutation 720 
    ## Permutation 740 
    ## Permutation 760 
    ## Permutation 780 
    ## Permutation 800 
    ## Permutation 820 
    ## Permutation 840 
    ## Permutation 860 
    ## Permutation 880 
    ## Permutation 900 
    ## Permutation 920 
    ## Permutation 940 
    ## Permutation 960 
    ## Permutation 980 
    ## Permutation 1000

``` r
summary(operm.covar, alpha=c(0.20, 0.05))
```

    ## LOD thresholds (1000 permutations)
    ##      lod
    ## 20% 3.00
    ## 5%  3.76

``` r
summary(out.covar, perms=operm.covar, alpha=0.05, pvalues=TRUE)
```

    ##           chr pos  lod pval
    ## c3.loc162   3 162 7.29    0

``` r
#           chr pos  lod pval
#c3.loc162   3 162 6.73    0
```

The QTL found in this genome scan is the next closest marker to the
original QTL and will not be considered further.

###### 2. Perform a genome scan with the marker as an **interactive** covariate

this allows us to detect loci that show an interaction with the chr3
locus

``` r
out.inter <- scanone(dat, pheno.col=2, model='binary', addcovar=g.matrix, intcovar=g.matrix)

#run another permutation test, this time with both addcovar and intercovar
operm.inter <- scanone(dat, pheno.col=2, model='binary', addcovar=g.matrix, intcovar=g.matrix, chr=-3, n.perm=1000)
```

    ## Permutation 20 
    ## Permutation 40 
    ## Permutation 60 
    ## Permutation 80 
    ## Permutation 100 
    ## Permutation 120 
    ## Permutation 140 
    ## Permutation 160 
    ## Permutation 180 
    ## Permutation 200 
    ## Permutation 220 
    ## Permutation 240 
    ## Permutation 260 
    ## Permutation 280 
    ## Permutation 300 
    ## Permutation 320 
    ## Permutation 340 
    ## Permutation 360 
    ## Permutation 380 
    ## Permutation 400 
    ## Permutation 420 
    ## Permutation 440 
    ## Permutation 460 
    ## Permutation 480 
    ## Permutation 500 
    ## Permutation 520 
    ## Permutation 540 
    ## Permutation 560 
    ## Permutation 580 
    ## Permutation 600 
    ## Permutation 620 
    ## Permutation 640 
    ## Permutation 660 
    ## Permutation 680 
    ## Permutation 700 
    ## Permutation 720 
    ## Permutation 740 
    ## Permutation 760 
    ## Permutation 780 
    ## Permutation 800 
    ## Permutation 820 
    ## Permutation 840 
    ## Permutation 860 
    ## Permutation 880 
    ## Permutation 900 
    ## Permutation 920 
    ## Permutation 940 
    ## Permutation 960 
    ## Permutation 980 
    ## Permutation 1000

``` r
summary(operm.inter, alpha=c(0.20, 0.05))
```

    ## LOD thresholds (1000 permutations)
    ##      lod
    ## 20% 3.69
    ## 5%  4.32

``` r
summary(out.inter, perms=operm.inter, alpha=0.05, pvalues=TRUE)
```

    ##           chr pos  lod  pval
    ## c9.loc104   9 104 5.17 0.005

``` r
#plot differences in LOD scores from analysis with marker as interactive
#covariate from analysis with marker as just additive covariate 
plot(out.inter - out.covar, ylab="interaction LOD score")
```

![](Figs/interactive-1.png)

``` r
plot(out.inter - out.covar, chr=c(3,4), ylab="interaction LOD score") 
```

![](Figs/interactive-2.png)

Conclusion: there is possibly a QTL on chr4 that interacts with the QTL
on chr3

### FIT A NEW QTL MODEL

Based on the evidence from the previous genome scan, create a QTL model
that includes QTL on chr3 and chr4, with an interaction between the two
QTL

``` r
qtl.c3c4 <- makeqtl(dat, what='prob', chr=c(3,4), pos=c(157, 20.4))

#fit new model with interaction
out.fqc3intc4 <- fitqtl(dat, pheno.col=2, qtl=qtl.c3c4, formula=y~Q1*Q2, method = "hk", model="binary")
summary(out.fqc3intc4)
```

    ## 
    ##      fitqtl summary
    ## 
    ## Method: Haley-Knott regression 
    ## Model:  binary phenotype
    ## Number of observations : 340 
    ## 
    ## Full model result
    ## ----------------------------------  
    ## Model formula: y ~ Q1 + Q2 + Q1:Q2 
    ## 
    ##       df      LOD     %var Pvalue(Chi2)
    ## Model  8 67.07792 59.68896            0
    ## 
    ## 
    ## Drop one QTL at a time ANOVA table: 
    ## ----------------------------------  
    ##                df    LOD   %var Pvalue(Chi2)    
    ## 3@157.0         6 66.601 59.045      < 2e-16 ***
    ## 4@20.4          6  4.195  2.357      0.00366 ** 
    ## 3@157.0:4@20.4  4  3.346  1.869      0.00393 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

The LOD scores for the other terms (QTL on chr4 and interaction term)
are not huge (will be compared to penalized LOD scores later in
multi-qtl analysis), and the p-values are not to be trusted because they
don’t account for the scan across the genome.

### TWO-DIMENSIONAL GENOME SCAN USING TWO-QTL MODEL

checks for linked QTL and QTL with limited marginal effects

2D scan using binary model (this takes a long time since it looks at
each pair of chromosomes separately. May consider running on a server.)

``` r
out2.bin <- scantwo(dat, pheno.col=2, model='binary')
save(out2.bin, file="out2bin.RData")
```

load the data generated from the 2D scan

``` r
load("out2bin.RData")
```

Perform permutation test so we can make sense of results. Because this
is very computationally challenging and extremely time consuming, we do
the significance tests clustered and on a server.  
To avoid using the same seed in each set of permutations, we call
set.seed separately for each set.

The script QTL\_permtests.R includes code for loading the R/qtl package,
reading in the csv file, and calculating conditional genotype
probablities, since this info is needed to run the permutation tests.

Run the script using the following bash code

``` bash
Rscript QTL_permtests.R
```

The QTL\_permtests.R script is included here:

``` r
library(qtl)
susc <- read.cross('csv', file='QTL_SNP_map.csv', genotypes=c('AA', 'AB', 'BB'))
summary(susc)

dat <- calc.genoprob(susc, step=1, error.prob=0.01)

set.seed(85842518)
operm2A <- scantwo(dat, model='binary', pheno.col=2, n.perm=100)
save(operm2A, file="perm2A.RData")

set.seed(85842519)
operm2B <- scantwo(dat, model='binary', pheno.col=2, n.perm=100)
save(operm2B, file="perm2B.RData")

set.seed(85842520)
operm2C <- scantwo(dat, model='binary', pheno.col=2, n.perm=100)
save(operm2C, file="perm2C.RData")

set.seed(85842521)
operm2D <- scantwo(dat, model='binary', pheno.col=2, n.perm=100)
save(operm2D, file="perm2D.RData")

set.seed(85842522)
operm2E <- scantwo(dat, model='binary', pheno.col=2, n.perm=100)
save(operm2E, file="perm2E.RData")

set.seed(85842523)
operm2F <- scantwo(dat, model='binary', pheno.col=2, n.perm=100)
save(operm2F, file="perm2F.RData")

set.seed(85842524)
operm2G <- scantwo(dat, model='binary', pheno.col=2, n.perm=100)
save(operm2G, file="perm2G.RData")

set.seed(85842525)
operm2H <- scantwo(dat, model='binary', pheno.col=2, n.perm=100)
save(operm2H, file="perm2H.RData")

set.seed(85842526)
operm2I <- scantwo(dat, model='binary', pheno.col=2, n.perm=100)
save(operm2I, file="perm2I.RData")

set.seed(85842527)
operm2J <- scantwo(dat, model='binary', pheno.col=2, n.perm=100)
save(operm2J, file="perm2J.RData")
```

After running the previous script, load the files produced from the
permutation tets.

``` r
load("perm2A.RData")
load("perm2B.RData")
load("perm2C.RData")
load("perm2D.RData")
load("perm2E.RData")
load("perm2F.RData")
load("perm2G.RData")
load("perm2H.RData")
load("perm2I.RData")
load("perm2J.RData")


#Combine the 10 batches of permutation tests into one object
operm2 <- c(operm2A, operm2B, operm2C, operm2D, operm2E, operm2F, operm2G, operm2H, operm2I, operm2J)
```

Use the permutation tests to obtain thresholds (expected LOD scores
under null model), and then compare with observed LOD scores

``` r
summary(operm2)
```

    ## susc (1000 permutations)
    ##     full  fv1  int  add  av1  one
    ## 5%  9.25 7.16 5.98 6.40 3.61 3.68
    ## 10% 8.76 6.73 5.64 5.93 3.38 3.32

``` r
#Look at pairs of QTL from 2D scan that exceed 5% thresholds
summary(out2.bin, perms=operm2, alpha=0.05) 
```

    ##     There were no pairs of loci meeting the criteria.

``` r
#Look at pairs of QTL from 2D scan that exceed 10% thresholds
summary(out2.bin, perms=operm2, alpha=0.1)
```

    ##     There were no pairs of loci meeting the criteria.

### STUDY EFFECTS OF PUTATIVE LINKED LOCI

Study effect of putative linked loci on chr3 and chr4 using plotPXG and
effectplot. effectplot requires imputed genotype data. Will rerun
sim.geno, this time with more n.draws (settings taken from Rqtl guide)
and only on chr3 and chr4.

``` r
datc3c4 <- sim.geno(subset(dat, chr=c(3,4)), step=2.5, error.prob=0.001, n.draws=256)

#find markers
marc3 <- find.marker(datc3c4, "3", 157)
marc4 <- find.marker(datc3c4, "4", 20.4)
par(mfrow=c(1,2)) #this line makes it so both plots appear together

#plot pxg
plotPXG(datc3c4, pheno.col=2, c(marc3, marc4))

#make effect plot using closest markers
effectplot(datc3c4, pheno.col=2, mname1=marc3, mname2=marc4)
```

![](Figs/linked-1.png) Plots show that these putative QTL have effects
of opposite sign (a given genotype (BB in this case) has opposite
effects at different loci)–this means they are linked in repulsion…loci
in repulsion will exhibit little marginal effect, but if considered
together they may stand out.

### MULTIPLE QTL MODEL

brings together all putative QTL and interactions

##### STEP 1: MODEL SELECTION

based on scanone and scantwo results

**Model \#1:** single qtl on chr3 y\~Q1 (based on initial scanone
results) **Model \#2:** interactive with qtl on chr 3,4 y\~Q1:Q2 (based
on scanone with <a href="mailto:chr3@157" class="email">chr3@157</a> as
covariate)

##### STEP 2: MODEL FIT

already performed previously

**Model \#1:** “qtl” **Model \#2:** “qtl.c3c4”

##### STEP 3: MODEL SEARCH

explore and fit additional models using several tools–manual approach

**Model 1**

------------------------------------------------------------------------

scan for additional (additive) QTL to add to model

``` r
out.aq1 <- addqtl(dat, pheno.col=2, qtl=qtl, formula=y~Q1, method="hk", model="binary")
plot(out.aq1, ylab="LOD score")
```

![](Figs/model%201%20addqtl-1.png)

``` r
max(out.aq1)
```

    ##          chr pos  lod
    ## c9.loc97   9  97 2.32

Use permutation tests to find LOD thresholds. While our permutation
results from scanone do not formally apply in the present case,in which
we are controlling for the locus on chr3, they nevertheless provide a
reasonable guide:

``` r
summary(operm.bin, alpha=c(0.20, 0.10, 0.05))
```

    ## LOD thresholds (1000 permutations)
    ##      lod
    ## 20% 3.01
    ## 10% 3.35
    ## 5%  3.76

The LOD score of 2.3196999 is not even above the 20% significance
threshold, so I will \#not consider this locus further.

Scan for additional loci that interact with chr3 locus

``` r
out.aqi <- addqtl(dat, pheno.col=2, qtl=qtl, formula=y~Q1*Q2, method="hk", model="binary")

plot(out.aqi, ylab="LOD score")
```

![](Figs/model1%20interact-1.png)

``` r
max(out.aqi)
```

    ##                    chr pos  lod
    ## scaffold02116_1534   9 113 5.64

``` r
summary(operm2)
```

    ## susc (1000 permutations)
    ##     full  fv1  int  add  av1  one
    ## 5%  9.25 7.16 5.98 6.40 3.61 3.68
    ## 10% 8.76 6.73 5.64 5.93 3.38 3.32

For thresholds, I use numbers in column fv1 (compares interactive model
to single QTL model) to analyze out.aqi. None of the peaks are above
even the 10% thresholds, so I do not consider them further.

**Model 2**

------------------------------------------------------------------------

scan for additional (additive) QTL to add to model

``` r
out.aq2 <- addqtl(dat, pheno.col=2, qtl=qtl.c3c4, formula=y~Q1*Q2, method="hk", model="binary")
plot(out.aq2, chr=4, ylab="LOD score")
```

![](Figs/model2%20addqtl-1.png)

``` r
max(out.aq2)
```

    ##          chr pos  lod
    ## c9.loc98   9  98 2.23

same interpretation as with model \#1. The LOD score of 2.2300964 is not
even above the 20% significance threshold from scanone.

Scan for additional loci that interact with chr3 locus

``` r
out.aqi2.1 <- addqtl(dat, pheno.col=2, qtl=qtl.c3c4, formula=y~Q1*Q2+Q1*Q3, method="hk", model="binary")

plot(out.aqi2.1, chr=4, ylab="LOD score")
```

![](Figs/model2%20interact-1.png)

``` r
max(out.aqi2.1)
```

    ##                    chr pos  lod
    ## scaffold02116_1534   9 113 5.63

No peaks above fv1 threshold

Scan for additional loci that interact with chr4 locus

``` r
out.aqi2.2 <- addqtl(dat, pheno.col=2, qtl=qtl.c3c4, formula=y~Q1*Q2+Q2*Q3, method="hk", model="binary")
plot(out.aqi2.2, ylab="LOD score")
```

![](Figs/model2%20intereact%20chr4-1.png)

``` r
max(out.aqi2.2)
```

    ##                    chr  pos  lod
    ## scaffold01702_4063   3 43.9 5.43

No peaks above fv1 threshold

##### STEP 4: MODEL COMPARISON

based on automated search algorithm stepwiseqtl

\#STEP 4: MODEL COMPARISON \#(based on automated search algorithm
stepwiseqtl)

\#first calculate LOD penalties using results from first permutation
test print(pen.hk \<- calc.penalties(operm2)) \# main heavy light
\#3.678864 5.984202 3.479083 calc.penalties(operm2, alpha=c(0.05, 0.10,
0.20)) \# main heavy light \#5% 3.678864 5.984202 3.479083 \#10%
3.324984 5.641420 3.407786 \#20% 3.011405 5.272679 3.314256

first calculate LOD penalties using results from first permutation test

``` r
print(pen.hk <- calc.penalties(operm2)) #default threshold is alpha = 0.05
```

    ##     main    heavy    light 
    ## 3.678864 5.984202 3.479083

``` r
#view penalties at multiple thresholds
calc.penalties(operm2, alpha=c(0.05, 0.10, 0.20))
```

    ##         main    heavy    light
    ## 5%  3.678864 5.984202 3.479083
    ## 10% 3.324984 5.641420 3.407786
    ## 20% 3.011405 5.272679 3.314256

Perform a stepwise qtl analysis…

…starting from **null model**

``` r
outsw.null.hk <- stepwiseqtl(dat, pheno.col=2, method="hk", model="binary", max.qtl=8, penalties=pen.hk, verbose=FALSE, keeplodprofile=TRUE, keeptrace=TRUE)
outsw.null.hk
```

    ##   QTL object containing genotype probabilities. 
    ## 
    ##       name chr pos n.gen
    ## Q1 3@157.0   3 157     3
    ## 
    ##   Formula: y ~ Q1 
    ## 
    ##   pLOD:  59.204

``` r
plotLodProfile(outsw.null.hk, ylab="Profile LOD score")
```

![](Figs/stepwise%20null-1.png)

``` r
#compares, at each position for a given QTL, the model with the QTL of interest
#at that particular position (and with the positions of all other QTL fixed at
#their maximum likelihood estimates) to the model with the QTL of interest omitted

#plot graphical representation of models visited
thetrace.nhk <- attr(outsw.null.hk, "trace")
par(mfrow=c(3,5))
for(i in seq(along=thetrace.nhk))
  plotModel(thetrace.nhk[[i]], chronly=TRUE,main=paste(i, ": pLOD =",round(attr(thetrace.nhk[[i]], "pLOD"), 2)))
```

![](Figs/stepwise%20null-2.png)![](Figs/stepwise%20null-3.png)

…starting from **Model 1**

``` r
outsw.sw1.hk <- stepwiseqtl(dat, pheno.col=2, method="hk", model="binary", qtl= qtl, formula=y~Q1, max.qtl=8, penalties=pen.hk, verbose=FALSE, keeplodprofile=TRUE, keeptrace=TRUE)
outsw.sw1.hk
```

    ##   QTL object containing genotype probabilities. 
    ## 
    ##       name chr pos n.gen
    ## Q1 3@157.0   3 157     3
    ## 
    ##   Formula: y ~ Q1 
    ## 
    ##   pLOD:  59.204

``` r
#same result as before

plotLodProfile(outsw.sw1.hk, ylab="Profile LOD score")
```

![](Figs/stepwise%20model%201-1.png)

``` r
thetrace.1hk <- attr(outsw.sw1.hk, "trace")
par(mfrow=c(3,5))
for(i in seq(along=thetrace.1hk))
  plotModel(thetrace.1hk[[i]], chronly=TRUE,main=paste(i, ": pLOD =",round(attr(thetrace.1hk[[i]], "pLOD"), 2)))
```

![](Figs/stepwise%20model%201-2.png)![](Figs/stepwise%20model%201-3.png)

…starting from **Model 2**

``` r
outsw.sw2.hk <- stepwiseqtl(dat, pheno.col=2, method="hk", model="binary", qtl= qtl.c3c4, formula=y~Q1*Q2, max.qtl=8, penalties=pen.hk, verbose=FALSE, keeplodprofile=TRUE, keeptrace=TRUE)
outsw.sw2.hk
```

    ##   QTL object containing genotype probabilities. 
    ## 
    ##       name chr pos n.gen
    ## Q1 3@157.0   3 157     3
    ## 
    ##   Formula: y ~ Q1 
    ## 
    ##   pLOD:  59.204

``` r
plotLodProfile(outsw.sw2.hk, ylab="Profile LOD score")
```

![](Figs/stepwise%20model%202-1.png)

``` r
thetrace.2hk <- attr(outsw.sw2.hk, "trace")
par(mfrow=c(3,5))
for(i in seq(along=thetrace.2hk))
  plotModel(thetrace.2hk[[i]], chronly=TRUE,main=paste(i, ": pLOD =",round(attr(thetrace.2hk[[i]], "pLOD"), 2)))
```

![](Figs/stepwise%20model%202-2.png)![](Figs/stepwise%20model%202-3.png)
So the final model is the single QTL model–with one large-effect QTL on
chr3

Get estimated effects of final model

``` r
out.fq.ests <- fitqtl(dat, pheno.col=2, model='binary', qtl=qtl, get.ests=TRUE, dropone=FALSE)
summary(out.fq.ests)
```

    ## 
    ##      fitqtl summary
    ## 
    ## Method: Haley-Knott regression 
    ## Model:  binary phenotype
    ## Number of observations : 340 
    ## 
    ## Full model result
    ## ----------------------------------  
    ## Model formula: y ~ Q1 
    ## 
    ##       df     LOD     %var Pvalue(Chi2)
    ## Model  2 62.8827 57.33205            0
    ## 
    ## 
    ## Estimated effects:
    ## -----------------
    ##               est      SE      t
    ## Intercept -1.7907  0.4558 -3.928
    ## 3@157.0a   4.0148  0.8851  4.536
    ## 3@157.0d  -0.8597  0.9854 -0.873

**Conclusion: The best model includes 1 QTL with LOD score 63, which
explains 57% of the phenotypic variance**

------------------------------------------------------------------------

Test for epistasis between loci D and F
=======================================

### FIT A NEW QTL MODEL

Create a QTL model that includes F locus and D locus with an interaction
between the two QTL. D locus is between markers scaffold01547\_75 (101.6
cM) and scaffold02269\_730 (105.3 cM). Bento et al 2020 uses the latter
marker to test for epistasis, so I will use that position

``` r
qtl.DandF <- makeqtl(dat, what='prob', chr=c(3,9), pos=c(157, 105.3))

#examine QTL object
qtl.DandF
```

    ##   QTL object containing genotype probabilities. 
    ## 
    ##       name chr    pos n.gen
    ## Q1 3@157.0   3 157.00     3
    ## Q2 9@105.3   9 105.31     3

``` r
plot(qtl.DandF)
```

![](Figs/DandF-1.png)

Fit model with interaction

``` r
out.DandF<- fitqtl(dat, pheno.col=2, qtl=qtl.DandF, formula=y~Q1*Q2, method = "hk", model="binary")
summary(out.DandF)
```

    ## 
    ##      fitqtl summary
    ## 
    ## Method: Haley-Knott regression 
    ## Model:  binary phenotype
    ## Number of observations : 340 
    ## 
    ## Full model result
    ## ----------------------------------  
    ## Model formula: y ~ Q1 + Q2 + Q1:Q2 
    ## 
    ##       df      LOD     %var Pvalue(Chi2)
    ## Model  8 67.36413 59.84493            0
    ## 
    ## 
    ## Drop one QTL at a time ANOVA table: 
    ## ----------------------------------  
    ##                 df    LOD   %var Pvalue(Chi2)    
    ## 3@157.0          6 65.115 56.844      < 2e-16 ***
    ## 9@105.3          6  4.481  2.513      0.00213 ** 
    ## 3@157.0:9@105.3  4  2.478  1.371      0.02231 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Make effect plot to visualize putative interaction

``` r
library(devEMF)#version 3.6
emf(file = "DandF.emf")
DandF_eff <- effectplot(dat2, pheno.col=2,  mname1 = "c9.loc105", mname2="c3.loc157", ylab="Proportion susccetible", xlab="", legend.lab = "", main="")
DandF_eff
```

    ## $Means
    ##             3@157.0.AA 3@157.0.AB 3@157.0.BB
    ## 9@105.0.AA 4.56941e-05  0.1181821  0.7362739
    ## 9@105.0.AB 0.00000e+00  0.1217341  0.9791759
    ## 9@105.0.BB 0.00000e+00  0.1700226  0.8982706
    ## 
    ## $SEs
    ##            3@157.0.AA 3@157.0.AB 3@157.0.BB
    ## 9@105.0.AA 0.06708215 0.07630757 0.04416324
    ## 9@105.0.AB 0.05785191 0.04443653 0.03111537
    ## 9@105.0.BB 0.06862411 0.05605997 0.04937382

``` r
dev.off()
```

    ## png 
    ##   2
