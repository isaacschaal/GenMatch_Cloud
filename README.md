# GenMatch_Cloud
Use for testing and tutorials for Genetic Matching executed in AWS EC2.

Clone this repo to your EC2 server

## Instructions

First install the packages you dont have installed, and then load the packages. Run all the code in the functions folder. To create the figure, we are only using the ``extractDataFromGenoud`` function and the ``graphGold`` functions. 

The first extracts all of the children form genmatch and finds the lowest p-value and treatment effect at each child. The population size and generation setting can be changed. The second graphs the output, with optional y axis and optional true treatment effect lines.

To get the figure, run the Latest.r file, which shoud have the correct settings of ``extractDataFromGenoud`` and ``graphGold`` to reproduce the figure, though it may take a few tries. 

## TODO

Add the function that randomally populates bad matches. It currently has hardcoded balance matrix etc, and needs to be made universalisable. 
