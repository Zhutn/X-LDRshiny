
# X-LD Shiny
The core algorithm of X-LD is implemented in C and C++ and combined with Rshiny to build interactive web apps straight from R.
You can easily deploy X-LD directly at you PC/sever by running the following commands in RStudio/R.
## Quick Start Guide
### Run the commands below to initialize X-LD Shiny.
~~~
# This is an R console
# 1. First, make sure the following R packages are installed before running
install.packages(c("shiny","bsplus","ggplot2","reshape2"))
# 2. Then, you can run it directly at you RStudio, if you have shiny package installed.
library(shiny)
runGitHub("X-LDRshiny", "huangxin0221")
~~~
### Other ways to initialize X-LD Shiny.
Normally for most of the users working with **Windows/MacOS/Unbuntu etc.**, there should be a window or browser tab pops up, X-LD is then ready for the analysis.
~~~
# 1. Download the source codes of X-LD Shiny
Find the download link at the homepage Code button, download the zipped source codes.
# 2. Unzip the soucrce codes locally
# 3. Open R/Rstudio and run the commands below to initialize X-LD Shiny.
setwd("/home/your_name/X-LDRshiny-main/")
library(shiny)
runApp()
~~~
For some of the **linux OS without displays**, LAN remote access can be used for the X-LD Shiny.
~~~
# This is an R console

# Assume your linux IP is 100.100.100.1, and port 1234 is accessible
library(shiny)
runApp(host='100.100.100.1',port=1234)
# X-LD is then availble at computer in LAN, by visit 100.100.100.1:1234 in the browser.


# The way to open the specific port (for example 1234) in linux:
# This is a linux bash/terminal
# Root privileges are required
# make sure the firewall is active
systemctl status firewalld
systemctl start firewalld

# open port 1234
firewall-cmd --zone=public --add-port=1234/tcp --permanent

# note
# --zone: scope zone
# --permanent: set options permanently, a change will only be part of the runtime configuration without this option

# reload firewall configuration
firewall-cmd --reload

~~~
### Some possible errors and their solutions 
~~~
# For Windows users
Error in utils::download.file(url, method = method, ...) : 
  cannot open URL 'https://github.com/huangxin0221/X-LDRshiny/archive/master.tar.gz' # Do not use runGitHub function to initialize X-LD Shiny 

# For Mac users

# For linux users
Error: GDK_BACKEND does not match available displays  # Normal error reporting without any adjustment
~~~
All platforms support.

