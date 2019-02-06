#extractDataFromGenoud <- function(x, balancematrix, pop, y, treatment, df){
extractDataFromGenoud <- function(X, BalanceMatrix, Y, Tr, df, cl, balance_formul, pop){
  
  # Create temp directory for the genoud.pro
  temp <- toString(tempdir())
  
  # Do genout, mout and mb for one set
  genout <- GenMatch(Tr=Tr, X=X, BalanceMatrix=BalanceMatrix, estimand="ATT", M=1,
                     pop.size=50, max.generations=10, wait.generations=5, cluster = cl, print.level = 3,
                     project.path = cat(temp,"genoud.pro"), nboots = 500)
  # print level >3 means it prints all gens, and file is made
  # path is the temp dir where genoud.pro goes
  
  # genout <- GenMatch(Tr=Tr, X=X, BalanceMatrix = BalanceMatrix,
  #                    estimand="ATT",
  #                    pop.size=10, max.generations=10, wait.generations=5, cluster = cl)
  
  mout <- Match(Y=Y, Tr=Tr, X=X, estimand="ATT", Weight.matrix=genout)
  summary(mout)
  
  mb <- MatchBalance(balance_formul,
                     match.out=mout, nboots=500)
  
  # Overwrite the previous genoud, write the full output of genout
  # also, is taking genoud.pro and turning to genoud.txt
  file.copy(paste(gsub("//", "/", temp),"/genoud.pro", sep = ""),"./genoud.txt",overwrite = TRUE)
  
  # load the file, named as con
  con <- file("genoud.txt", "r")
  
  # output reads the lines of the file
  output <- readLines(con)
  
  #### Output looks like this
  
  # first line is "Generation: 0 \t Population Size: 20 \t Fit Values: 22 \t Variables: 10"
  # then ""
  # then 1-20 (b/c one for each member of pop)
  # then 32 numbers -
  # first 22 are the $value of genout, then the last 10 are the weights
  # this repeats for each child
  # then repeat for each generation
  
  ####
  
  ## drops all of the rows that are spaces or the generation
  for(row in output){
    if(substr(row,1,1) == "G" || substr(row,1,1) == ""){
      index <- match(row, output)
      output <- output[-index]
    }
  }
  
  # delets all the tabs
  output <- gsub(" \t ", ",",output)
  # then gets a list for of each line
  outputlist <- strsplit(output,"[,]")
  
  # create a new dataframe
  # The first col will be the lowest p val, the last col will be the treatment effect
  # and the middle cols will be the weights of the weight matrix
  cols = ncol(BalanceMatrix)
  df <- matrix(data = NA, ncol = cols+1, nrow = length(outputlist))
  #This was updated (ncol should be 1 more than cols)
  
  for(candidate in c(1:length(outputlist))){
    
    # 1st # is the member of the pop
    # 2nd number is the lowest p val ( the next 21 are the other p vals (2 for each) )
    # 23rd is the 1st of weights, till 32
    
    # Get the lowest p val and the weight matrix
    selected <- outputlist[[candidate]][c(2,(cols*2+2):(cols*2+ncol(X)+1))]
    
    # Convert them to numeric
    selection <- as.numeric(selected)
    
    # set the first 11 to what we have selected
    df[candidate,c(1:cols)] <- selection
  }
  
  for(i in c(1:nrow(df))){
    # matching using the weight matrix
    # diag makes the list turn into diag of matrix
    mout <- Match(Y=Y, Tr=Tr, X=X, estimand="ATT", Weight.matrix=diag(df[i,c(2:(ncol(df)-1))]))
    # changed the y and x and treat
    ### CHANGE TO have the correct inputs
    
    # put it in the last col
    df[i,ncol(df)] <- mout$est
    
  }
  return (df)
}

###### DATA SETUP

# Get the Delija Wahba sample
nsw_dw <- read_dta("nsw_dw.dta")
View(nsw_dw)

# Get the CPS-1 Data
cps_controls <- read_dta("cps_controls.dta")
View(cps_controls)

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


output1 <- extractDataFromGenoud(X, BalanceMatrix, Y, Tr, new_df, cl, balance_formul, pop)
graphGold(output1, known_treat_effect = 1.794)

output2 <- extractDataFromGenoud(X, BalanceMatrix, Y, Tr, new_df, cl, balance_formul, pop)
graphGold(output2)

output3 <- extractDataFromGenoud(X, BalanceMatrix, Y, Tr, new_df, cl, balance_formul, pop)
graphGold(output3)

##### TESTING

