library(ExcessMass)

#very simple example to test the work of localmax & excessm. 
#Excess mass should be 0.4 (taking end point correction into account)
x=c(0,3,4,5,8)
excessm(x,0.1)

#Another easy to calculate example
#excess mass increases from 0.45 under allowing for one mode to 0.75 under two modes
x=c(1,2,2,3,3,9,9,10,11,11)
excessm(x,0.05, M=2)

#Testing exmplot using the same data as above.
#Getting the corresponding values by hand or by using the excessm function

Lambda=c(0.025,0.05,0.075,0.01)
exmplot(x, M=2, Lambda=Lambda)

##Note that there is no need in testing exmsilhouette & mexmsilhouette,
#as they use the data of excessm and visualize it
#The validity of search max Lambda can be seen in the plot representation