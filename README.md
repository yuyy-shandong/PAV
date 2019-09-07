# PAV

PAV(Post Auxiliary Variables regression),is an R package for efficient causal effect estimation of environmental exposure on health outcome. Just taking advantage of a post-outcome exposure as an auxiliary variable, the PAV regression can obtain an unbiased and robust causal effect estimation of exposure on outcome. 


# Installation
It is easy to install the development version of PAV package using the 'devtools' package. The typical install time on a "normal" desktop computer is less than one minute.

```
# install.packages("devtools")
library(devtools)
install_github("yuyy-shandong/PAV")
```


# Usage
There is only one functions in PAV package.
You can find the instructions by '?PAV'. 

library(PAV)

?PAV



# Example


```
u <- rnorm(1000,0,1)
x1 <- 0.5*u +rnorm(1000,0,1)
x3 <- 0.5*u +rnorm(1000,0,1)
y <- 2*x1 + 1*u +rnorm(1000,0,1)
data <- data.frame (u,x1,x3,y)
model <- PAV (y~x1+x3, data = data, sdmethod ="normal",x1.name = "x1",x2.name = "x3",boots.no = NULL)
model
```


# Development
This R package is developed by Yuanyuan Yu and Fuzhong Xue.


