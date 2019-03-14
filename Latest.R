###### DATA SETUP #######

# Get the Delija Wahba sample
nsw_dw <- read_dta("nsw_dw.dta")
#View(nsw_dw)

# Get the CPS-1 Data
cps_controls <- read_dta("cps_controls.dta")
#View(cps_controls)

# Get a name of the dimensions for creating the new_df
Full_dims <-c("treat","re78","age","education","black","hispanic",
              "married","nodegree","re74","re75")

# create the new_df
new_df_1 <- nsw_dw[,Full_dims]
new_df_2 <- cps_controls[,Full_dims]
new_df <- rbind(new_df_1, new_df_2)


# Attach the data
attach(new_df)

# Get treatment and outcome
Tr <- treat
Y <- re78/1000 #from the thing
pop = 100 # not really used yet

# X is a matrix of all dimensions except re78 and treat
X <- cbind(age, education, black, hispanic, married, nodegree, re75, re74)

# BalanceMat is matrix of all dimensions except re78 and treat, along with 
# a new dim, called v11, which is re74 * re75
BalanceMatrix <- cbind(age, education, black, hispanic, married, nodegree, re74, re75,
                       I(age*education), I(age*black), I(age*hispanic), I(age*married), I(age*nodegree),
                       I(age*re74), I(age*re75), I(education*black), I(education*hispanic),
                       I(education*married), I(education*nodegree), I(education*re74), I(education*re75),
                       I(black*hispanic), I(black*married), I(black*nodegree), I(black*re74),
                       I(black*re75), I(hispanic*married), I(hispanic*nodegree), I(hispanic*re74),
                       I(hispanic*re75), I(married*nodegree), I(married*re74), I(married*re75),
                       I(nodegree*re74), I(nodegree*re75), I(re74*re75))

#Balance formula for use in the Match Balance function
balance_formul <- treat~ age + education + black + hispanic + married + nodegree + re74 + re75 +
                      I(age*education) + I(age*black) + I(age*hispanic) + I(age*married) +
                      I(age*nodegree) + I(age*re74) + I(age*re75) + I(education*black) +
                      I(education*hispanic) + I(education*married) + I(education*nodegree) +
                      I(education*re74) + I(education*re75) + I(black*hispanic) + I(black*married) +
                      I(black*nodegree) + I(black*re74) + I(black*re75) + I(hispanic*married) +
                      I(hispanic*nodegree) + I(hispanic*re74) + I(hispanic*re75) + I(married*nodegree) +
                      I(married*re74) + I(married*re75) + I(nodegree*re74) + I(nodegree*re75) + I(re74*re75)


# Execute in Parrallel

# Use 1 less core than computer has, and make a cluster
detectCores()
no_cores <- detectCores()-1
cl <- makeCluster(no_cores)


##### FUNCTION RUNNING

output <- extractDataFromGenoud(X, BalanceMatrix, Y, Tr, new_df, cl, balance_formul, pop)
graphGold(output, known_treat_effect = 1.794)



