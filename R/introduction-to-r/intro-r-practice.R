#Dive right into the R
#######################################################



# Using R as Calculator
2 + 8
sum(2,8)

2*4+6-5/2.5
(2*4)+6-5/2.5

sqrt(25)
5^2
5**2


#Logical operations
4==2
4==4
4>=3
4> 5

a = 5
b = sqrt(a*runif(1)*10)

if (a > b) {print("a is greater than b")} else {print("a is less than b")}

#Objects and Variables

a = 5
a = 67
a + b
a = a + b 


name <- c("Patrick", "Gregory", "Bernard", "Lesla", "Bridget", "Rico", "Temwa", "Andrew", "Cecil","Martha","Merriam", "William", "Martha","Mada","Sara",NA)
age  <- c(23,29,31,21,34,38,28,33,25,NA,30,35,33,30,29,NA)
sex  <- c(1,1,1,2,2,1,2,1,2,2,2,1,2,2,2,1)

#Note that you can just type the object without the print function
print(name) 
age
sex

#We can get number of elements in each vector/object/variable
length(name)
length(age)
length(sex)

#Data frames
#Possibly, interest would be to see corresponding values of each individual side by side. These may not be the best, but they at least save some notable purpose.


cbind(name, age, sex) #Combines the variables column-wise
rbind(name, age, sex) #Combines the variables low-wise. 
data.frame(name, age, sex) #Combines the variables column-wise. This is the best.

#Saving the data frame to an object called mydata

mydata <- data.frame(name, age, sex)
mydata

# In R, "NA" implies missing value (Not applicable). Currently, we would say that we don't yet have Martha's age. And there is a certain male with no name and age.
# During analysis of data with missing values, a decision has to be made either to remove or impute them in the analysis (not necessarily the dataset) 

#creating and inserting new variables in a data frame

#Note now that we have two Marthas, and we may want to distinguish between the two. We may want to give them sequential id numbers. Before we do that, let's look at this

#1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16
c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
1:16
1:length(name)
seq(1,16)
seq(1,16,2)


# Now, we know how we can create the ids, and we need to attach them to the corresponding names. We can first create id <- 1: length(name) and then attach it to the data frame mydata, ie as mydata <- data.frame(id,mydata). Alternatively, this process can be done directly, i.e
mydata <- data.frame(id = 1:length(name), mydata)
mydata


# Getting a feel of the data
head(mydata)      #print the first default 6 rows
head(mydata,5)    #print the first n=5 rows
tail(mydata)
tail(mydata,3)

#Indexing elements
#There are times we want information for specific individuals

mydata[4,1]
mydata[4,2]
mydata[c(4),]

#What if you wanted to infor for person on row number 1 and 2. What's your guess?

mydata[c(5,4),c(2,4)]
mydata[c(5,4,5),c(2,4)]   #Everything is possible ;)
mydata[c(4),2] 
mydata[c(4,1),2]
mydata[4,]
mydata[c(4),]       #If you don't specify any column after the comma, it selects all available columns
mydata[,c(2,4)]     #Similarly, if you don't specify any row, it selects all available rows


# And we can specify name of the column, we don't always have to remember the column number of our variable of interest

mydata[1:3,"name"]
mydata[1:3,c("name", "sex")]
mydata[,"name"]


mydata[mydata$sex==1,]
mydata[mydata$sex==2,]
mydata[mydata$sex==2,c("name","sex","age")]
mydata$age[mydata$sex==1]   #Print ages for males

#Some Summaries of the data

mean(mydata$age[mydata$sex==1], na.rm = T)  #Avearge age for males
mean(mydata$age[mydata$sex==2], na.rm = T)  #Avearge age for females
aggregate(age~sex, data = mydata, mean)     #Print average for sex categories

#Sorting/Ordering data.

a <- c(7,4,9,1)
sort(a) #This sorts individual elements in ascending order
order(a)


mydata[c(2,6,3),]
mydata[order(mydata$sex),] #Sort data by sex
mydata[order(mydata$sex, mydata$age),] #Sort data by sex and age

#Working inside the data frame
#Some data summaries

rm(age, name, sex, a, b)
mydata$age
mean(mydata[,"age"])
mean(mydata[,"age"], na.rm = T) 
mean(mydata["age"], na.rm = T)  #Why won't this work?


head(mydata$age,7)
mydata$age[6]
mydata$age[c(6,3)]
mean(mydata$age, na.rm = T)
table(mydata$sex)

#Frequencies and proportions


1/sum(c(1,1,2)) #First element
1/sum(c(1,1,2)) #Second elemenet
2/sum(c(1,1,2)) #Third element

c(1,1,2)/4      #How cool!!
c(1,1,2)/sum(c(1,1,2))
a <- c(1,1,2)
a/sum(a)

prop.table(1)
prop.table(c(1,1))
prop.table(c(1,1,2))  #I hope now you get the idea. Check this

a <- c(5,8,6)
prop.table(a) 

#What do you make of these below?
a <- c(1,1,2)
prop.table(a) 
prop.table(table(a))
#Hopefully, we get the idea! 
#Now, we get back to our data,

prop.table(table(mydata$sex))

#Basic graphs

hist(mydata$age)
barplot(table(mydata$sex))  #Barplot of frequencies of sex categories

barplot(prop.table(table(mydata$sex)))  #Barplot of proportions of sex categories

#Merging data
#Lets create some data to illustrate the concept

id <- sample(1:20, 7, replace = F )
ed <- factor(sample(1:3,7,replace = T), levels = c(1,2,3), labels = c("none","formal","informal"))

edu.data <- data.frame(id, ed)

#Merge using common variable "id"

merge(mydata,edu.data, by = "id")
merge(mydata,edu.data, by = "id", all.x = T)
merge(mydata,edu.data, by = "id", all.y = T)
merge(mydata,edu.data, by = "id", all = T)



#Working with Characteristics

#grepl("characters","string")

grepl("test","testing")
grepl("test","est")
grepl("test","Test")


#Practice 1

candidates <- read.csv("candidates.csv")     #Candidates data
names(candidates) <- c("time","sex", "nationality","status","function_use","tidyverse","need_analyse","like_learn","learning_mode")

#Q01: What is the distribution of nationality?
#Q02: What is the proportion of males opting for online mode of learning?
#Q03: What is the proportion of candidates who uses R function by sex?

cbind(table(candidates$origin_country))




table(cbind(table(candidates$nationality)))

table(candidates$nationality)

candidates$origin_country <- candidates$nationality
candidates$origin_country[grepl("ambia",candidates$nationality)] <- "Zambian"

cbind(table(candidates$origin_country))

candidates$origin_country[grepl("alawi",candidates$nationality)] <- "Malawian"
candidates$origin_country[grepl("ALAWI",candidates$nationality)] <- "Malawian"


#Candidates, Diamonds Dataset

#Model building: Regression: Linear Regression, ANOVA, 

diamonds <- read.csv("diamond.csv") 

#Section: More useful functions

#Section: More about graphs (ggplot)

#Section: User-defined functions


#Section: Tidyverse

#Merging data
#Lets create some data to illustrate the concept
#########################################################################


#Section: More useful functions
#Section: More about graphs
#Section: User-defined functions

#Modelling

diamond <- read.csv("diamond.csv")
View(diamond)
summary(diamond)
table(diamond$cut)
table(diamond$color)
table(diamond$clarity)
table(diamond$clarity, diamond$cut)
rowSums(table(diamond$clarity, diamond$cut))
row_tot = rowSums(table(diamond$clarity, diamond$cut))
col_tot = colSums(table(diamond$clarity, diamond$cut))
colSums(table(diamond$clarity, diamond$cut))
cbind(table(diamond$clarity, diamond$cut))
cbind(table(diamond$clarity, diamond$cut),col_tot)
cbind(table(diamond$clarity, diamond$cut),row_tot)
rbind(cbind(table(diamond$clarity, diamond$cut),row_tot),row_tot)
colSums(cbind(table(diamond$clarity, diamond$cut),row_tot),row_tot)
aggregate(price~color, data = diamond, mean)
aggregate(price~color+cut, data = diamond, mean)
aggregate(price~color+cut, data = diamond, length)
aggregate(price~color, data = diamond, length)
table(diamond$color)

#Correlation matrix
cor(diamond$price,diamond$depth)
cor(data.frame(diamond$x,diamond$y, diamond$z))
cor(data.frame(x = diamond$x,y = diamond$y, z = diamond$z))
cor(data.frame(carat = diamond$CARAT, x = diamond$x,y = diamond$y, z = diamond$z))


#Covariance matrix
cov(data.frame(diamond$x,diamond$y, diamond$z))
cov(data.frame(x = diamond$x,y = diamond$y, z = diamond$z))
cov(data.frame(carat = diamond$CARAT, x = diamond$x,y = diamond$y, z = diamond$z))

#Scatter Plots
plot(diamond$x,diamond$y)
plot(diamond$x,diamond$y, main = "Plot of Y against X")
plot(diamond$x,diamond$y, main = "Plot of Y against X", ylab = "y", xlab = "x")
plot(diamond$x,diamond$y, main = "Plot of Y against X", ylab = "y", xlab = "x", ylim = c(0,10))
plot(diamond$x,diamond$y, main = "Plot of Y against X", ylab = "y", xlab = "x", ylim = c(0,10), xlim = c(0,10))

plot(diamond$x,diamond$y, main = "Plot of Y against X", pch = 1)


#Anova

aov(price ~ cut, data = diamond)
summary(aov(price~cut, data = diamond))

my_anova2 <- aov(price~cut + color, data = diamond)
summary(my_anova2)

lm(price~cut, data = diamond)
summary(lm(price~cut, data = diamond))
table(diamond$cut)
summary(lm(price~cut+color, data = diamond))
summary(lm(price~cut+color+cut*color, data = diamond))


#Linear Regression using least squares


my_linear1 <- lm(price~cut, data = diamond)
summary(my_linear1)

my_linear2 <- lm(price~x+y+z + cut, data = diamond)
summary(my_linear2)

#Linear Regression with interaction
my_linear3 <- lm(price~x+y+z + cut + x*y, data = diamond) 
summary(my_linear3)



#Linear Regression using generalized linear models


my_linear1 <- glm(price~cut, data = diamond)
summary(my_linear1)

my_linear2 <- glm(price~x+y+z + cut, data = diamond)
summary(my_linear2)

#Linear Regression with interaction
my_linear3 <- glm(price~x+y+z + cut + x*y, data = diamond) 
summary(my_linear3)


#Logistic regresion

#Is job category associated with high/low salary? 

employee <- read.csv("employee.csv")
View(employee)
str(employee)
summary(employee)
employee$high_sal = 1
employee$high_sal[employee$salary < 28875] = 0



table(employee$jobcat,employee$high_sal)
prop.table(table(employee$jobcat,employee$high_sal))
prop.table(table(employee$jobcat,employee$high_sal),1)

colSums(prop.table(table(employee$jobcat,employee$high_sal),1))

chisq.test(table(employee$jobcat,employee$high_sal))


#By default, glm outputs linear regression
glm(high_sal~jobcat, data = employee)

#Binary logistic regression specified family
glm(high_sal~jobcat, data = employee, family = binomial)




# Tidyverse

#The tidyverse package actually contains other packages (dplyr, ggplot2, etc.) and you’ll see that when you load the tidyverse package using library(). Remember the package must be installed to your device before it can be loaded into your libraries! For help on installing packages, refer to Section 


library(tidyverse)


# pipes

#The pipe operator, (%>%), feeds the results of one operation into the next operation. It is more handy when there is a sequence of operations on a data frame. The advantage of using the pipe operator is that it makes code extremely easy to read.

# mutate, group_by, summarize, filter, select, arrange



set.seed(356863)
region <- sample(1:3,1000, prob = c(0.2,0.5,0.3), replace = T)
distr <- c(sample(1:7,200,replace = T),sample(1:10,504,replace = T),sample(1:15,296,replace = T))
resid <- c(rep(c(1,2),c(39,161)),rep(c(1,2),c(203,301)),rep(c(1,2),c(75,221)))
bdywgt <- rnorm(1000,33,9)

#Tibble (a data frame version)
sample.data <- tibble(region,distr, resid, bdywgt)
rm(region,distr,resid,bdywgt)

sample.data <- sample.data %>% 
  mutate(dist_code = region*100+distr)

#Summarize
sample.data %>% 
  group_by(region) %>% 
  summarise(meanwgt = mean(bdywgt), 
            freqn   = NROW(bdywgt))


seed <- read.csv("seed.csv")
str(seed)
seed[sample(1:nrow(seed),10),]

table(seed$Blocks)
table(seed$cultivar)
table(seed$seedchem)

mean(seed$response)
aggregate(response~Blocks,data=seed, mean)


my.aov <- aov(response~Blocks,data=seed)
summary(my.aov)

#Boxplot
boxplot(response~Blocks, data=seed)

#ggplot version
my.blocks.resp <-  ggplot(seed, aes(x = factor(Blocks), y = response)) 
my.blocks.resp + geom_boxplot()
my.blocks.resp + geom_boxplot() + coord_flip() + labs(x = "plot blocks", y = "yield", title = "Maize yield over plot blocks")

# ggplot2


set.seed(356863)
region <- sample(1:3,1000, prob = c(0.2,0.5,0.3), replace = T)
distr <- c(sample(1:7,200,replace = T),sample(1:10,504,replace = T),sample(1:15,296,replace = T))
resid <- c(rep(c(1,2),c(39,161)),rep(c(1,2),c(203,301)),rep(c(1,2),c(75,221)))
bdywgt <- rnorm(1000,33,9)

#Tibble (a data frame version)
sample.data <- tibble(region,distr, resid, bdywgt)
rm(region,distr,resid,bdywgt)

sample.data <- sample.data %>% 
  mutate(dist_code = region*100+distr)

#Summarize
sample.data %>% 
  group_by(region) %>% 
  summarise(meanwgt = mean(bdywgt), 
            freqn   = NROW(bdywgt))



seed <- read.csv("seed.csv")
str(seed)
seed[sample(1:nrow(seed),10),]

table(seed$Blocks)
table(seed$cultivar)
table(seed$seedchem)

mean(seed$response)
aggregate(response~Blocks,data=seed, mean)

my.aov <- aov(response~Blocks,data=seed)
summary(my.aov)

#Graphs
boxplot(response~Blocks, data=seed)


my.blocks.resp <-  ggplot(seed, aes(x = factor(Blocks), y = response)) 
my.blocks.resp + geom_boxplot()
my.blocks.resp + geom_boxplot() + coord_flip() + labs(x = "plot blocks", y = "yield", title = "Maize yield over plot blocks")


#With ggplot2, data and aesthetic mappings are supplied in ggplot(), then layers are added on with +. This is an important pattern, and as you learn more about ggplot2 you’ll construct increasingly sophisticated plots by adding on more types of components.

#Almost every plot maps a variable to x and y, so naming these aesthetics is tedious, so the first two unnamed arguments to aes() will be mapped to x and y. This means that the following code is identical to the example above:

#Linear regression

salaries <- read.csv("WorkSalaries.csv")
View(salaries)



mylm.plot <- ggplot(salaries, aes(yrs.service, salary))
mylm.plot + geom_point()



my.lm <- lm(salary~yrs.service, data = salaries)
summary(my.lm)

ls(my.lm)



mylm.plot +
  geom_smooth(method = lm, se = F, fullrange = T)











#Tidyverse applications

#Renaming varables
candidates <- read.csv("candidates.csv")
View(candidates)

candidates <- candidates %>% 
  rename( time = Timestamp,
          sex  = Gender,
          nationality = Nationality.Country.of.Origin,
          status = Which.of.the.following.qualifications.best.describes.you.,
          function_use = Do.you.know.how.to.use.functions.in.R.,
          tidyverse = Do.you.already.use.the.tidyverse.packages.in.R.such.as.dplyr..tidyr..tibble.and.ggplot2..,
          need_analyse = Do.you.need.to.analyse.large.collections.of.related.time.series.,
          like_learn = Would.you.like.to.learn.how.to.use.some.tidy.tools.for.time.series.analysis.including.visualization..decomposition.and.forecasting..,
          learning_mode = Which.of.the.two.modes.of.learning.will.be.suitable.for.you.)


#Select only variables of interest
candidates <- candidates %>% 
  select(-time) %>% 
  mutate(ind = 1,
         pers_id = 1:nrow(candidates))


#Calculate number of observations (hhs) by team
candidates <- candidates %>% 
  group_by(nationality) %>% 
  mutate(total_nat = sum(ind))



#Plotting with ggplot2
###########################################################



#Plotting sf objects with ggplot2 
############################################################





#Manipulate data

age_int <- c(0,1,seq(5,95,5))
nqx <- c(0.02592,0.0042,0.00232,0.00201,0.00443,0.00611,0.00632,0.00654,0.01098,0.01765,0.02765,0.04387,0.05987,0.09654,0.13654,0.18765,0.25439,0.37887,0.47898,0.57908,1)
lx <- c(100000)

for (i in 2:length(nqx))
{
  lx[i] <- round(lx[i-1] - lx[i-1]*nqx[i-1])  
}

ndx <- round(nqx * lx)


#Another example

x   <- c(0,1,seq(5,75,5))
n   <- c(1,4,rep(5,(length(x)-2)))
nMx <- c(0.1072,0.0034,0.0010,0.0007,0.0017,0.0030,0.0036,0.0054,0.0054,0.0146,0.0128,0.0269,0.0170,0.0433,0.0371,0.0785,0.0931)
nkx <- c(0.33,1.56,rep(2.5,length(nMx)-2))

nqx <- round((n*nMx)/(1 + (n - nkx)*nMx),4)

lx <- c(100000)
for (i in 2:length(nqx))
{
  lx[i] <- round(lx[i-1] - lx[i-1]*nqx[i-1]) 
  #Lx[i] <- (lx[i-1] + lx[i])*2.5
}

ndx <- round(nqx*lx)
nLx <- n * lx - ndx*(n - nkx)

Tx <- NA

for (i in 1:length(nqx))
{
  Tx[i] <- sum(nLx[i:length(nqx)])
  #Lx[i] <- (lx[i-1] + lx[i])*2.5
}

ex <- Tx/lx
as.data.frame(cbind(x,n,nMx,nkx, nqx,lx,ndx,nLx,Tx,ex))


