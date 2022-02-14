# X-LDRShiny
The core algorithm of X-LD is implemented in C and C++ and combined with Rshiny to build interactive web apps straight from R.
You can easily deploy EigenGWAS directly at you PC/sever by running the following commands in RStudio/R.
~~~
# Make sure the following R packages are installed before running
install.packages(c("shiny","bsplus","RColorBrewer","corrplot"))
# Then, you can run it directly at you RStudio, if you have shiny package installed.
library(shiny)
runGitHub("X-LDRshiny", "huangxin0221")
~~~
All platforms support.

