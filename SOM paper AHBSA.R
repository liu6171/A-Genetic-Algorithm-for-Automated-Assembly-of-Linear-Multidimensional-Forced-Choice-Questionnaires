# Supplemental Online Material: R Codes for the Real Form Assemblies in "A Genetic Algorithm for 
# Automated Assembly of General Multidimensional Linear Tests"

# set the working directory where the 'SOM paper AHBSA.Rdata' file is
# saved, for example, 'c:\GA'
setwd("c:/GA")

# load functions and data (workspace)
load("SOM paper AHBSA.Rdata")

library(eatATA)
library(mvtnorm)

# In the code below, the functions I wrote are explicitly specified or noted;
# all the other functions are from the eatATA package or R.  You can always use
# the 'dump' function to see the source code of a function. For example,

dump("AHBSA.sampling", file = "", control = NULL)

#All the objective funcitons and constraints used in the real forms 
#are described in the paper.

###################
#Real data
####################

#Real statement pool: SSlope=discriminate, SLoc=difficulty (2PLM)

head(statement.pool.real)

#  ItemID Direction Domain.ID Skill.ID Domain.Skill.ID SSlope   SLoc
#1  ASS01        NA         3        2               8  1.786 -0.500
#2  ASS02        NA         3        2               8  1.608 -0.841
#3  ASS03        NA         3        2               8  0.833 -1.157
#4  ASS04        NA         3        2               8  1.765 -0.538
#5  ASS05        -1         3        2               8  0.900  0.283
#6  ASS06        NA         3        2               8  2.307 -0.517

#Assemble real pair forms 

##################################################################
#Function to generate eligible pairs from a statement pool as input to 
#AHBSA function--AHBSA.sampling.  
#
#Description:
#Generate eligible pairs from a statement pool to be used as input to the AHBSA 
#function--AHBSA.sampling. 
#
#Usage:
#eligible.item.selection(statement.feature, quad=c(-2,0,2), neg.reverse=T) 
#
#Arguments:
#statement.feature -- a data frame including a statement pool
#quad -- a vector including the quadratures per dimension
#neg.reverse -- logic, T = negative statements are reversed scored in the statement pool. 
#
#Return:
#a 3-element list: (1) a Nstatement by Nstatement matrix with 1 indicating an eligible item 
#and 0 otherwise; (2) a Nstatement by Nstatement by 17 array including item features for each 
#eligible item; (3)  a Nstatement by Nstatement by Nquad by Nquad array including item information 
#at each quadrature point across two dimensions for each eligible  item.   
#
#Note: The funciton assumes the variable names in statement.pool.real. The criterion for
#an eligible item is that the two statements measure different domains.  
#
##########################################################################
eligible.pair.AHBSA <- eligible.item.selection(statement.pool.real, 
            quad = c(-1,0,1), neg.reverse = T)


########################################################################
#Function to assemble test forms using AHBSA. 
#
#Description:
#Assemble test forms with any block sizes using AHBSA. It can also assemble 
#random forms that only meet constraits. 
#
#Usage:
#AHBSA.sampling(eligible.item, constraint.fun, population.size, 
#bias.ratio=2^-4,theta.density.array, theta.type="ML" ,
#theta.cor.inv.array=NULL, converge.crit="CP", 
#random.draw=F, random.draw.time=NULL, ntrait=5, 
#item.direction.pattern.limit=NULL, skill.limit=rep(4,15), 
#skill.dir.limit=NULL, random.N.gen=NULL,domain.comb.limit=array(2, dim=rep(5,3)))
#
#Arguments:
#eligible.item -- an output from eligible.item.selection()
#constraint.fun -- function implementing constraits: constraint.fun.Likert for Likerts,
#			constraint.fun for pairs, constraint.fun.triplet for triplets 
#population.size -- size of populaiton 
#bias.ratio -- bias ratio in AHBSA sampling
#theta.density.array -- an array with Ntrait+1 dimensions (e.g., [3,3,3,3,3,5] for 5 traits each with 3 quadratures) storing normalized density of
#	each quadrature point across traits. The last dimension repeats the denisty Ntrait times.  
#theta.type -- type of latent score estimates: "MAP"=Maximum A Posteriori, "ML"= Maximum Likelihood
#theta.cor.inv.array -- an array with Ntrait+2 dimensions (e.g., [3,3,3,3,3,5,5] for 5 traits each with 3 quadratures) storing the inverse of intertrait
#	correlation matrix. The matrix is repeated for each quadrature point across traits. 			
#converge.crit -- convergence criterion: "CP" = all forms in a population are the same; numeric value = maximum numbers of generations
#random.draw -- logic, T = assemble random forms
#random.draw.time -- seconds before stoping random assembly 
#ntrait -- number of traits measured (dimensions)
#item.direction.pattern.limit -- define item direction constrait passed to constraint.fun() and constraint.fun.triplet(): 1 = 
#		each combination of Domain.ID.1, Domain.ID.2, and Dir.item must appear once (only used in constraint.fun()); 
#		a numeric vector defining the max number of items in each category of Dir.item. 
#skill.limit -- a numeric vector defining the max number of items in each skill level as used in the paper, passed to constraint.fun.triplet() and constraint.fun.Likert()
#skill.dir.limit -- a Nskill by 2 matrix defining the max number of items in each combination of skill and Dir.item, passed to constraint.fun.Likert()
#random.N.gen -- Number of random forms to be assembled
#domain.comb.limit -- a Ndomain by 3 matrix defining the max number of items in each combination of three domains, passed to constraint.fun.triplet()
# 
#Return:
#a list including all generations, and each generation is a 4-element list  
#1. a population.size-element list: each element is an decision array with block.size dimension, and each dimension includes Nstatement elements with 
#	1 = a statement or combination is selected, and 0 otherwise.
#2. a population.size-element vector including mean reliability across domains for each assembled form.
#3. a population.size by Ndomain matrix including test reliability for each domain in each assembled form. 
#4. cummulative running time (seconds) 
##########################################################################

#Assemble real pair forms by AHBSA
AHBSA.out.pair <- AHBSA.sampling(eligible.item=eligible.pair.AHBSA, 
            constraint.fun=constraint.fun, population.size = 50, test.length = 30, 
             bias.ratio = 2^-4, theta.density.array=theta.density.array,
             theta.type="MAP", theta.cor.inv.array=theta.cor.inv.array, 
            item.direction.pattern.limit = NULL, converge.crit = "CP")

#test reliability of each domain in the final real pair form 
AHBSA.out.rel.pair <- AHBSA.out.pair[[length(AHBSA.out.pair)]][[3]][1, ]

#Assemble 1000 random real pair forms
random.out.pair <- AHBSA.sampling(eligible.item=eligible.pair.AHBSA, 
            constraint.fun=constraint.fun, population.size = 1, test.length = 30, 
            bias.ratio = 2^-4,theta.density.array=theta.density.array,
             theta.type="MAP", theta.cor.inv.array=theta.cor.inv.array,
            random.draw = T, random.N.gen = 1000)
   
#mean test reliability of each domain across 1000 random real pair forms 
random.out.rel.pair<- apply(sapply(random.out.pair, function(x) x[[3]]), 1, mean)

#Assemble real pair forms by MIP

##################################################################
#Function to generate eligible pairs from the output of eligible.item.selection() 
#for MIP.  
#
#Description:
#Generate eligible pairs from the output of eligible.item.selection() to use as   
#an item pool for MIP. 
#
#Usage:
#MIP.item.pool.gen(eligible.pair.sim)
#
#Arguments:
#eligible.pair.sim -- a output from eligible.item.selection()
#
#Return:
#a data frame with the following item feature and dummy variables (see two example cases below in head(eligible.pair.MIP, 2)): 
#Dir.1, Dir.2 = Statements 1 and 2's directions (1=positive, -1=negative)
#Dir.item = item’s direction (1 = two positive statements, 2 = two negative statements, 3 = one positive and one negative statements)
#ItemID.1, ItemID.2 = Statements 1 and 2's labels
#StatementID.1, StatementID.2 = Statements 1 and 2’s IDs 
#Domain.ID.1, Domain.ID.2 = Statements 1 and 2’s domain IDs (1–5) 
#Domain.Skill.ID.1, Domain.Skill.ID.2 = Statements 1 and 2’s skill IDs (1–15, with Skills 1–3 belonging to Domain 1, Skills 4–6 belonging to Domain 2, and so on)
#SSlope.1, SSlope.2 = Statements 1 and 2’s discrimination parameters 
#SLoc.1, SLoc.2 = Statements 1 and 2’s difficulty parameters
#Info.Sum = sum of the five domains’ expected item information
#pair.ID = item’s ID 
#Statex = dummy variable for Statement ID x (x = 1–91)
#pair.skill.domain1 = combination of Domain.Skill.ID.1 and Domain.ID.2
#pair.skill.domain2 = combination of Domain.Skill.ID.2 and Domain.ID.1
#pair.domain.dir = combination of Domain.ID.1, Domain.ID.2, and Dir.item
#DSx = dummy variable for a domain-skill combination in a pair (x = 1–60)
#
##########################################################################

eligible.pair.MIP <- MIP.item.pool.gen(eligible.pair.AHBSA)

head(eligible.pair.MIP, 2)

#  Dir.1 Dir.2 Dir.item ItemID.1 Domain.ID.1 Skill.ID.1 Domain.Skill.ID.1
#1     1    -1        3    EMP01           1          1                 1
#2     1    -1        3    EMP03           1          1                 1
#  SSlope.1 SLoc.1 StatementID.1 ItemID.2 Domain.ID.2 Skill.ID.2
#1     1.27 -2.742             1    RES01           2          1
#2    1.045 -2.622             2    RES01           2          1
#  Domain.Skill.ID.2 SSlope.2 SLoc.2 StatementID.2   info.sum pair.ID State1
#1                 4   -0.842 -0.738            19 0.06698692       1      1
#2                 4   -0.842 -0.738            19 0.08727068       2      0
#  State2 State3 State4 State5 State6 State7 State8 State9 State10 State11
#1      0      0      0      0      0      0      0      0       0       0
#2      1      0      0      0      0      0      0      0       0       0
#  State12 State13 State14 State15 State16 State17 State18 State19 State20
#1       0       0       0       0       0       0       0       1       0
#2       0       0       0       0       0       0       0       1       0
#  State21 State22 State23 State24 State25 State26 State27 State28 State29
#1       0       0       0       0       0       0       0       0       0
#2       0       0       0       0       0       0       0       0       0
#  State30 State31 State32 State33 State34 State35 State36 State37 State38
#1       0       0       0       0       0       0       0       0       0
#2       0       0       0       0       0       0       0       0       0
#  State39 State40 State41 State42 State43 State44 State45 State46 State47
#1       0       0       0       0       0       0       0       0       0
#2       0       0       0       0       0       0       0       0       0
#  State48 State49 State50 State51 State52 State53 State54 State55 State56
#1       0       0       0       0       0       0       0       0       0
#2       0       0       0       0       0       0       0       0       0
#  State57 State58 State59 State60 State61 State62 State63 State64 State65
#1       0       0       0       0       0       0       0       0       0
#2       0       0       0       0       0       0       0       0       0
#  State66 State67 State68 State69 State70 State71 State72 State73 State74
#1       0       0       0       0       0       0       0       0       0
#2       0       0       0       0       0       0       0       0       0
#  State75 State76 State77 State78 State79 State80 State81 State82 State83
#1       0       0       0       0       0       0       0       0       0
#2       0       0       0       0       0       0       0       0       0
#  State84 State85 State86 State87 State88 State89 State90 State91
#1       0       0       0       0       0       0       0       0
#2       0       0       0       0       0       0       0       0
#  pair.skill.domain1 pair.skill.domain2 pair.domain.dir1 DS1 DS2 DS3 DS4 DS5
#1                1_2                4_1            1_2_3   1   0   0   1   0
#2                1_2                4_1            1_2_3   1   0   0   1   0
#  DS6 DS7 DS8 DS9 DS10 DS11 DS12 DS13 DS14 DS15 DS16 DS17 DS18 DS19 DS20 DS21
#1   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0
#2   0   0   0   0    0    0    0    0    0    0    0    0    0    0    0    0
#  DS22 DS23 DS24 DS25 DS26 DS27 DS28 DS29 DS30 DS31 DS32 DS33 DS34 DS35 DS36
#1    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
#2    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
#  DS37 DS38 DS39 DS40 DS41 DS42 DS43 DS44 DS45 DS46 DS47 DS48 DS49 DS50 DS51
#1    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
#2    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
#  DS52 DS53 DS54 DS55 DS56 DS57 DS58 DS59 DS60 pair.domain.dir
#1    0    0    0    0    0    0    0    0    0           1_2_3
#2    0    0    0    0    0    0    0    0    0           1_2_3

#MIP objective function: maximizing Info.Sum
 test.info.obj.pair <- maxObjective(
		nForms = 1, 
		itemValues = eligible.pair.MIP$info.sum, 
            weight = 1, 
		itemIDs = eligible.pair.MIP$pair.ID)

#Constraint 1: a statement appears at most once
statement.con.pair <- 
combineConstraints(lapply(1:91, function(s){
itemValuesRangeConstraint(
  nForms=1,
  itemValues=eligible.pair.MIP[, paste0("State",s)],
  range=c(0,1),
  itemIDs = eligible.pair.MIP$pair.ID
)}))

#Constraint 2: a skill pairs with a domain at most once
skill.domain.con.pair  <- combineConstraints(lapply(1:60, function(s){
itemValuesRangeConstraint(
  nForms=1,
  itemValues=eligible.pair.MIP[, paste0("DS",s)],
  range=c(0,1),
  itemIDs = eligible.pair.MIP$pair.ID
)}))

#MIP solver in eatATA
MIP.out.time.pair <- system.time(MIP.out.pair <- useSolver(
	list(test.info.obj.pair, statement.con.pair, skill.domain.con.pair), 
      solver = "GLPK"))

#inspect the final form
MIP.out.pair.data <- inspectSolution ( MIP.out.pair , items = eligible.pair.MIP , idCol = "pair.ID")
MIP.out.pair.data[[1]]

########################################################################
#Function to calculate reliabilities of test forms output by useSolver(). 
#
#Description:
#calculate test reliabilities of ML or MAP score estimates for forms output 
#by useSolver().
#
#Usage:
#MIP.reliability( MIP.out, item.pool, Nitem.test, Nstate, eligible.item, 
#theta.density.array, theta.type="MAP", theta.cor.inv.array, block.size=2, 
#idCol="pair.ID")
#
#Arguments:
#MIP.out -- an output from useSolver()
#item.pool -- item pool used in useSolver()
#Nitem.test -- number of items in a form
#Nstate -- number of statement in the statement pool
#eligible.item -- an output from eligible.item.selection()
#theta.density.array -- an array with Ntrait+1 dimensions (e.g., [3,3,3,3,3,5] for 5 traits each with 3 quadratures) storing normalized density of
#	each quadrature point across traits. The last dimension repeats the denisty Ntrait times.  
#theta.type -- type of latent score estimates: "MAP"=Maximum A Posteriori, "ML"= Maximum Likelihood
#theta.cor.inv.array -- an array with Ntrait+2 dimensions (e.g., [3,3,3,3,3,5,5] for 5 traits each with 3 quadratures) storing the inverse of intertrait
#	correlation matrix. The matrix is repeated for each quadrature point across traits. 
#block.size -- block size of test form: 1, 2, or 3
#idCol -- character, ID column name passed to inspectSolution()
#Return:
#a three-element list: 
#1. a vector including reliability estimate for each trait;
#2. an array with Ntrait+1 dimensions (e.g., [3,3,3,3,3,5] for 5 traits each with 3 quadratures) storing the variance of each trait at
#	each quadrature point across traits.
#3. an array with Ntrait+2 dimensions (e.g., [3,3,3,3,3,5,5] for 5 traits each with 3 quadratures) storing the test information matrix at
#	each quadrature point across traits. 
##########################################################################

MIP.rel.pair <- MIP.reliability(MIP.out = MIP.out.pair, 
            item.pool = eligible.pair.MIP, Nitem.test = 30, 
            Nstate=91, eligible.item = eligible.pair.AHBSA, theta.density.array=theta.density.array, 
            theta.type="MAP", theta.cor.inv.array=theta.cor.inv.array)

#Assemble real triplet forms

#Assemble real triplet forms by AHBSA

##################################################################
#Function to generate eligible triplets from a statement pool as input to 
#AHBSA function--AHBSA.sampling.  
#
#Description:
#Generate eligible triplets from a statement pool to be used as input to the AHBSA 
#function--AHBSA.sampling. 
#
#Usage:
#eligible.triplet.selection(statement.feature, quad=c(-2,0,2), neg.reverse=T) 
#
#Arguments:
#statement.feature -- a data frame including a statement pool
#quad -- a vector including the quadratures per dimension
#neg.reverse -- logic, T = negative statements are reversed scored in the statement pool. 
#
#Return:
#a 3-element list: (1) a Nstatement by Nstatement by Nstatement array with 1 indicating an eligible item 
#and 0 otherwise; (2) a Nstatement by Nstatement by Nstatement by 25 array including item features for each 
#eligible item; (3)  a Nstatement by Nstatement by Nstatement by Nquad by Nquad by Nquad array including item infomration 
#at each quadrature point across three dimensions for each eligible item.   
#
#Note: The funciton assumes the variable names in statement.pool.real. The criterion for
#an eligible item is that the three statements measure three different domains.  
#
##########################################################################


eligible.triplet.AHBSA <- eligible.triplet.selection(statement.pool.real, 
            quad = c(-1,0,1), neg.reverse = T)

#Assemble real triplet forms by AHBSA
AHBSA.out.triplet <- AHBSA.sampling(eligible.item=eligible.triplet.AHBSA, 
            constraint.fun=constraint.fun.triplet, population.size = 50, 
            test.length = 20, bias.ratio = 2^-4, theta.density.array=theta.density.array,
             theta.type="MAP", theta.cor.inv.array=theta.cor.inv.array, item.direction.pattern.limit = NULL, 
            skill.limit = rep(4, 15), domain.comb.limit = array(2, 
                dim = rep(5, 3)), converge.crit = "CP")

#test reliability of each domain in the final real triplet form 
AHBSA.out.rel.triplet <- AHBSA.out.triplet[[length(AHBSA.out.triplet)]][[3]][1, ]

#Assemble 1000 random real triplet forms      
      
random.out.triplet <- AHBSA.sampling(eligible.item=eligible.triplet.AHBSA, 
            constraint.fun=constraint.fun.triplet, population.size = 1, 
            test.length = 20, bias.ratio = 2^-4, theta.density.array=theta.density.array,
             theta.type="MAP", theta.cor.inv.array=theta.cor.inv.array,  
            random.draw = T, random.N.gen = 1000, item.direction.pattern.limit = NULL, 
            skill.limit = rep(4, 15), domain.comb.limit = array(2, dim = rep(5, 3)))

#mean test reliability of each domain across 1000 random real triplet forms    
random.out.rel.triplet <- apply(sapply(random.out.triplet, 
            function(x) x[[3]]), 1, mean)

#Assemble real triplet forms by MIP

##################################################################
#Function to generate eligible triplets from the output of eligible.triplet.selection() 
#for MIP.  
#
#Description:
#Generate eligible triplets from the output of eligible.triplet.selection() to use as   
#an item pool for MIP. 
#
#Usage:
#MIP.item.pool.gen.triplet(eligible.triplet)
#
#Arguments:
#eligible.triplet -- a output from eligible.triplet.selection()
#
#Return:
#a data frame with the following item feature and dummy variables (see two example cases below in head(eligible.triplet.MIP, 2)): 
#Dir.1, Dir.2, Dir.3 = Statements 1, 2, and 3's directions (1=positive, -1=negative)
#Dir.item = item’s direction (1 = three positive statements, 2 = three negative statements, 3 = two positive and one negative statements, 4 = one positive and two negative statements)
#ItemID.1, ItemID.2, ItemID.3 = Statements 1, 2, and 3's labels
#StatementID.1, StatementID.2, StatementID.3 = Statements 1, 2, and 3’s IDs 
#Domain.ID.1, Domain.ID.2, Domain.ID.3 = Statements 1, 2, and 3’s domain IDs (1–5) 
#Domain.Skill.ID.1, Domain.Skill.ID.2, Domain.Skill.ID.3 = Statements 1, 2, and 3’s skill IDs (1–15, with Skills 1–3 belonging to Domain 1, Skills 4–6 belonging to Domain 2, and so on)
#SSlope.1, SSlope.2, SSlope.3 = Statements 1, 2, and 3’s discrimination parameters 
#SLoc.1, SLoc.2,  SLoc.3 = Statements 1, 2, and 3’s difficulty parameters
#Info.Sum = sum of the five domains’ expected item information
#triplet.ID = item’s ID 
#Statex = dummy variable for Statement ID x (x = 1–91)
#DSx = dummy variable for a combination of Domain.ID.1, Domain.ID.2, Domain.ID.3, 
#	and one of the 9 skills in the previous three domains (x = 1–90)
#triplet.domain = combination of Domain.ID.1, Domain.ID.2, and Domain.ID.3
#
##########################################################################

eligible.triplet.MIP <- MIP.item.pool.gen.triplet(eligible.triplet.AHBSA)

head(eligible.triplet.MIP, 2)
#  Dir.1 Dir.2 Dir.3 Dir.item ItemID.1 Domain.ID.1 Skill.ID.1 Domain.Skill.ID.1
#1     1    -1     1        3    EMP01           1          1                 1
#2     1    -1     1        3    EMP03           1          1                 1
#  SSlope.1 SLoc.1 StatementID.1 ItemID.2 Domain.ID.2 Skill.ID.2
#1     1.27 -2.742             1    RES01           2          1
#2    1.045 -2.622             2    RES01           2          1
#  Domain.Skill.ID.2 SSlope.2 SLoc.2 StatementID.2 ItemID.3 Domain.ID.3
#1                 4   -0.842 -0.738            19    ENE01           3
#2                 4   -0.842 -0.738            19    ENE01           3
#  Skill.ID.3 Domain.Skill.ID.3 SSlope.3 SLoc.3 StatementID.3  info.sum
#1          1                 7    1.007 -2.591            38 0.5892980
#2          1                 7    1.007 -2.591            38 0.5629433
#  triplet.ID State1 State2 State3 State4 State5 State6 State7 State8 State9
#1          1      1      0      0      0      0      0      0      0      0
#2          2      0      1      0      0      0      0      0      0      0
#  State10 State11 State12 State13 State14 State15 State16 State17 State18
#1       0       0       0       0       0       0       0       0       0
#2       0       0       0       0       0       0       0       0       0
#  State19 State20 State21 State22 State23 State24 State25 State26 State27
#1       1       0       0       0       0       0       0       0       0
#2       1       0       0       0       0       0       0       0       0
#  State28 State29 State30 State31 State32 State33 State34 State35 State36
#1       0       0       0       0       0       0       0       0       0
#2       0       0       0       0       0       0       0       0       0
#  State37 State38 State39 State40 State41 State42 State43 State44 State45
#1       0       1       0       0       0       0       0       0       0
#2       0       1       0       0       0       0       0       0       0
#  State46 State47 State48 State49 State50 State51 State52 State53 State54
#1       0       0       0       0       0       0       0       0       0
#2       0       0       0       0       0       0       0       0       0
#  State55 State56 State57 State58 State59 State60 State61 State62 State63
#1       0       0       0       0       0       0       0       0       0
#2       0       0       0       0       0       0       0       0       0
#  State64 State65 State66 State67 State68 State69 State70 State71 State72
#1       0       0       0       0       0       0       0       0       0
#2       0       0       0       0       0       0       0       0       0
#  State73 State74 State75 State76 State77 State78 State79 State80 State81
#1       0       0       0       0       0       0       0       0       0
#2       0       0       0       0       0       0       0       0       0
#  State82 State83 State84 State85 State86 State87 State88 State89 State90
#1       0       0       0       0       0       0       0       0       0
#2       0       0       0       0       0       0       0       0       0
#  State91 DS1 DS2 DS3 DS4 DS5 DS6 DS7 DS8 DS9 DS10 DS11 DS12 DS13 DS14 DS15
#1       0   1   0   0   1   0   0   1   0   0    0    0    0    0    0    0
#2       0   1   0   0   1   0   0   1   0   0    0    0    0    0    0    0
#  DS16 DS17 DS18 DS19 DS20 DS21 DS22 DS23 DS24 DS25 DS26 DS27 DS28 DS29 DS30
#1    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
#2    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
#  DS31 DS32 DS33 DS34 DS35 DS36 DS37 DS38 DS39 DS40 DS41 DS42 DS43 DS44 DS45
#1    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
#2    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
#  DS46 DS47 DS48 DS49 DS50 DS51 DS52 DS53 DS54 DS55 DS56 DS57 DS58 DS59 DS60
#1    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
#2    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
#  DS61 DS62 DS63 DS64 DS65 DS66 DS67 DS68 DS69 DS70 DS71 DS72 DS73 DS74 DS75
#1    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
#2    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
#  DS76 DS77 DS78 DS79 DS80 DS81 DS82 DS83 DS84 DS85 DS86 DS87 DS88 DS89 DS90
#1    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
#2    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
#  triplet.domain
#1          1_2_3
#2          1_2_3

#MIP objective function: maximizing Info.Sum
test.info.obj.triplet <- maxObjective(nForms = 1, itemValues = eligible.triplet.MIP$info.sum, 
            weight = 1, itemIDs = eligible.triplet.MIP$triplet.ID)

#Constraint 1: a statement appears at most once
statement.con.triplet <- 
combineConstraints(lapply(1:91, function(s){
itemValuesRangeConstraint(
  nForms=1,
  itemValues=eligible.triplet.MIP[, paste0("State",s)],
  range=c(0,1),
  itemIDs = eligible.triplet.MIP$triplet.ID
)}))

#Constraint 2: a skill appears in a combination of three domains at most once
skill.domain.con.triplet <- 
combineConstraints(lapply(1:90, function(s){
itemValuesRangeConstraint(
  nForms=1,
  itemValues=eligible.triplet.MIP[, paste0("DS",s)],
  range=c(0,1),
  itemIDs = eligible.triplet.MIP$triplet.ID
)}))

#Constraint 3: a combination of three domains appears twice in a form
domain.nlevles.triplet <- length(levels(eligible.triplet.MIP$triplet.domain))
domain.con.triplet <- 
itemCategoryRangeConstraint(
  nForms=1,
  itemCategories=eligible.triplet.MIP$triplet.domain,
  range=matrix(2, nc=2, nrow=domain.nlevles.triplet), 
  itemIDs = eligible.triplet.MIP$triplet.ID
)

#MIP solver in eatATA
MIP.out.time.triplet <- system.time(MIP.out.triplet <- useSolver(list(test.info.obj.triplet, 
            statement.con.triplet, skill.domain.con.triplet, domain.con.triplet), 
            solver = "GLPK"))

#inspect the final triplet form
MIP.out.triplet.data <- inspectSolution ( MIP.out.triplet , items = eligible.triplet.MIP , idCol = "triplet.ID")
MIP.out.triplet.data[[1]]

#calculate test reliabilities of the final triplet form 
MIP.rel.triplet <- MIP.reliability(MIP.out = MIP.out.triplet, 
            item.pool = eligible.triplet.MIP, Nitem.test = 20, 
            Nstate=91, eligible.item = eligible.triplet.AHBSA,  theta.density.array=theta.density.array, 
            theta.type="MAP", theta.cor.inv.array=theta.cor.inv.array, block.size=3, idCol="triplet.ID")


#Assemble real Likert forms 

##################################################################
#Function to generate statement feature data from a statement pool to be used in
#AHBSA function--AHBSA.sampling.  
#
#Description:
#Generate statement feature data from a statement pool to be used in
#AHBSA function--AHBSA.sampling.
#
#Usage:
#Likert.item.pool.gen(statement.pool.feature)
#
#Arguments:
#statement.pool.feature -- a data frame including a statement pool
#
#Return:
#a data frame containing the following variables: ItemID, Direction, Domain.ID, 
Skill.ID, Domain.Skill.ID, SSlope, SLoc, StatementID, Dir.item       
#
#Note: The funciton assumes the variable names in statement.pool.real.  
#
##########################################################################
Likert.feature <- Likert.item.pool.gen(statement.pool.real)

head(Likert.feature, 2)
#     ItemID  Direction Domain.ID Skill.ID Domain.Skill.ID SSlope  SLoc    
#[1,] "EMP01" "1"       "1"       "1"      "1"             "1.27"  "-2.742"
#[2,] "EMP03" "1"       "1"       "1"      "1"             "1.045" "-2.622"
#     StatementID Dir.item
#[1,] "1"         "1"     
#[2,] "2"         "1"     

#calculate item information for each statement  
Likert.info <- t(apply(matrix(as.numeric(Likert.feature[, 
            c("SSlope", "SLoc")]), nc = 2), 1, info.2PL.fun, 
            theta = c(-1,0,1), simplify = T))
 
#create a list as the input of the eligible.item argument in AHBSA.sampling()
eligible.Likert.AHBSA<- list(array(rep(1,91), dim=91),Likert.feature, Likert.info)

#Assemble real Likert form by AHBSA
AHBSA.out.Likert <-  AHBSA.sampling(eligible.item=eligible.Likert.AHBSA, constraint.fun=constraint.fun.Likert, population.size=50, 
	test.length=60,bias.ratio=2^-4,theta.density.array=theta.density.array,
      theta.type="MAP", theta.cor.inv.array=theta.cor.inv.array, skill.limit=rep(4,15), converge.crit = "CP")

#test reliability of each domain in the final real Likert form  
AHBSA.out.rel.Likert <- AHBSA.out.Likert[[length(AHBSA.out.Likert)]][[3]][1,]

#Assemble 1000 random real Likert forms
random.out.Likert <-  AHBSA.sampling(eligible.item=eligible.Likert.AHBSA, constraint.fun=constraint.fun.Likert, population.size=1, 
	test.length=60,bias.ratio=2^-4,theta.density.array=theta.density.array,
      theta.type="MAP", theta.cor.inv.array=theta.cor.inv.array,skill.limit=rep(4,15), random.draw=T, random.N.gen=1000)

#mean test reliability of each domain across 1000 random real Likert forms 
random.out.Likert <- apply(sapply(random.out.Likert, function(x) x[[3]]), 1, mean)

#Assemble real Likert forms by MIP

##################################################################
#Function to generate a Likert statement pool from the input of the eligible.item argument  
#in AHBSA.sampling() for MIP.  
#
#Description:
#Generate a Likert statement pool from the input of the eligible.item argument  
#in AHBSA.sampling() for MIP.  
#
#Usage:
#MIP.Likert.item.pool.gen(eligible.Likert, quad)
#
#Arguments:
#eligible.Likert -- an input of the eligible.item argument in AHBSA.sampling()
#quad -- a vector including the quadratures per dimension
#
#Return:
#a data frame containing the following variables (see two example cases below in head(eligible.Likert.MIP, 2)): 
#ItemID, Direction, Domain.ID, Skill.ID, Domain.Skill.ID, SSlope, SLoc, StatementID, Dir.item, Skill.dir, info.sum     
#
##########################################################################

eligible.Likert.MIP <- MIP.Likert.item.pool.gen(eligible.Likert.AHBSA, quad=c(-1,0,1))

head(eligible.Likert.MIP,2)
#  ItemID Direction Domain.ID Skill.ID Domain.Skill.ID SSlope   SLoc StatementID
#1  EMP01         1         1        1               1   1.27 -2.742           1
#2  EMP03         1         1        1               1  1.045 -2.622           2
#  Dir.item Skill.dir   info.sum
#1        1       1_1 0.11627348
#2        1       1_1 0.08090184

#MIP objective function: maximizing Info.Sum
test.info.obj.Likert <- maxObjective(
  	nForms=1,
  	itemValues=eligible.Likert.MIP$info.sum,
  	weight = 1,
  	itemIDs = eligible.Likert.MIP$StatementID
)

#Constraint 1: a statement appears at most once
Likert.StatementID.factor <- as.factor(eligible.Likert.MIP$StatementID)
Likert.StatementID.factor.nlevels <- length(levels(Likert.StatementID.factor))

statement.con.Likert <- 
itemCategoryMaxConstraint(
  nForms=1,
  itemCategories=Likert.StatementID.factor,
  max=rep(1,Likert.StatementID.factor.nlevels),
  itemIDs = eligible.Likert.MIP$StatementID
)

#Constraint 2: a skill appears 4 times in a form
Likert.skill.domain.factor <- as.factor(eligible.Likert.MIP$Domain.Skill.ID)
Likert.skill.domain.nlevles <- length(levels(Likert.skill.domain.factor))

skill.domain.con.Likert <- 
itemCategoryRangeConstraint(
  nForms=1,
  itemCategories=Likert.skill.domain.factor,
  range=matrix(4, nc=2, nrow=Likert.skill.domain.nlevles), 
  itemIDs = eligible.Likert.MIP$StatementID
)

#MIP solver in eatATA
MIP.out.time.Likert <- system.time(MIP.out.Likert <- useSolver(list(test.info.obj.Likert, statement.con.Likert, skill.domain.con.Likert), 
	solver="GLPK"))

#inspect the final Likert form
MIP.out.Likert.data <- inspectSolution (MIP.out.Likert, items = eligible.Likert.MIP, idCol = "StatementID")
MIP.out.Likert.data[[1]]

#calculate test reliabilities of the final Likert form 
MIP.rel.Likert <- MIP.reliability(MIP.out = MIP.out.Likert, 
            item.pool = eligible.Likert.MIP, Nitem.test = 60, 
            Nstate=91, eligible.item = eligible.Likert.AHBSA, theta.density.array=theta.density.array, 
            theta.type="MAP", theta.cor.inv.array=theta.cor.inv.array, block.size=1, idCol="StatementID")