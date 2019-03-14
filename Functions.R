# Installing Packages and Imports

install.packages('Matching')
install.packages('foreign')
install.packages('rgenoud')

library(Matching)
library(foreign)
library(rbounds)
library(dplyr)
library(haven)
library(parallel)
library(rgenoud)

#####

# Extract Data From Genoud
# This function runs Genmatch, (with the specified genmatch params)
# and with the specified data

# It uses the info from Genoud to get the treatment effect and min p value
# for each child in each generation that were run. Thus, the final number of
# data points is popsize * generations run (only find this out after is run)

#####

extractDataFromGenoud <- function(X, BalanceMatrix, Y, Tr, df, cl, balance_formul, pop){
  
  # Create temp directory for the genoud.pro
  temp <- toString(tempdir())
  
  # Do genout, mout and mb for one set
  genout <- GenMatch(Tr=Tr, X=X, BalanceMatrix=BalanceMatrix, estimand="ATT", M=1,
                     pop.size=25, max.generations=10, wait.generations=5, cluster = cl, print.level = 3,
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
  cols = ncol(X)
  df <- matrix(data = NA, ncol = cols+2, nrow = length(outputlist)) # changed
  #This was updated (ncol should be 1 more than cols)
  end <- length(outputlist[[1]])
  start <- end - cols+1
  
  for(candidate in c(1:length(outputlist))){
    
    # 1st # is the member of the pop
    # 2nd number is the lowest p val ( the next 21 are the other p vals (2 for each) )
    # 23rd is the 1st of weights, till 32
    
    # Get the lowest p val and the weight matrix
    
    selected <- outputlist[[candidate]][c(2,start:end)]
    
    # Convert them to numeric
    selection <- as.numeric(selected)
    
    # set the first 11 to what we have selected
    for_sel = cols +1
    df[candidate,c(1:for_sel)] <- selection
  }
  
  for(i in c(1:nrow(df))){
    # matching using the weight matrix
    # diag makes the list turn into diag of matrix
    mout <- Match(Y=Y, Tr=Tr, X=X, estimand="ATT", Weight.matrix=diag(df[i,c(2:(ncol(df)-1))])) # its here i think
    # changed the y and x and treat
    ### CHANGE TO have the correct inputs
    
    # put it in the last col
    df[i,ncol(df)] <- mout$est
    
  }
  return (df)
}

#####

# Given as input a dataframe of points with P-value and Treatment effects
# (Note: Works with both Storage.df type data frames, as
# well as 2 col data frames with pval then TE)
# Prints a scatter plot of the points
# Can add the true treatment effect line
# Normally doesn't show y-axis (conceals the treatment effect), but can be shown if wanted

#####

graphGold <- function(df, known_treat_effect = NA, show_y_axis = FALSE, title = NA){
  
  # The last col is where the treatment effect is stored
  # and the first col is the lowest p val
  
  if (show_y_axis){
    plot(x=df[c(1:nrow(df)),1], y=df[c(1:nrow(df)),ncol(df)],
         ylab = "Treatment Effect", xlab = "Lowest P-value of Balance Test", main = title)
    
  # Don't show the y axis
  } else if(!show_y_axis) {
    plot(x=df[c(1:nrow(df)),1], y=df[c(1:nrow(df)),ncol(df)],
         ylab = "Treatment Effect", xlab = "Lowest P-value of Balance Test", yaxt = 'n', main = title)
  }
  
  # adds the treatment effect line
  if (!(is.na(known_treat_effect))){
    abline(h=known_treat_effect) 
  }
  
}

#####

# Executes one GenMatch function, returns the lowest p-value and treatment effect
# can be run in parrallel

#####


getATET_min_pval<-function(X, BalanceMatrix, Y, Tr, df, cl, balance_formul){
  
  # Input the X, Balance Matrix, Y Tr and df of the data
  # get the cluster that you have made
  # make the balance frmula match the balance matrix ## could do this automatically
  
  ## Should include changing the other hyper params
  
  # Do genmatch, mout, and mb
  genout <- GenMatch(Tr=Tr, X=X, BalanceMatrix = BalanceMatrix,
                     estimand="ATT",
                     pop.size=50, max.generations=20, wait.generations=5, cluster = cl)
  
  mout <- Match(Y=Y, Tr=Tr, X=X, estimand="ATT", Weight.matrix=genout)
  summary(mout)
  
  
  mb  <- MatchBalance(balance_formul,
                      data=new_df, match.out=mout, nboots=500)
  
  # Get the smallest p_value and the ATET from the matching
  smallest_p_value <- mb$AMsmallest.p.value
  ATET <- mout$est
  
  #return(c(smallest_p_value, ATET, mb, mout)) ## Use this to return the full results
  return(c(smallest_p_value, ATET))
  
}

