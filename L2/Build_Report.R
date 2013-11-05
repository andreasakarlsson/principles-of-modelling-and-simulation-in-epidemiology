# # File name
# fname <- 'L2'
# 
# # Set working directory
# setwd('/home/andkar/workspace/MEB/Modeling_in_Epidemilogy/L2')
# 
# # Load packages
# require(knitr)
# require(markdown)
# 
# # Create .md, .html, and .pdf files
# knit(paste(fname,'.Rmd',sep=""))
# markdownToHTML(paste(fname,'.md',sep=""), paste(fname,'.html',sep=""), options=c("use_xhml"))
# system(paste("pandoc -s ",fname,".html -o ",fname,".pdf",sep=""))

# Set working directory
setwd('/home/andkar/workspace/MEB/Modeling_in_Epidemilogy/L2')

# Load packages
require(knitr)
require(markdown)

# Create .md, .html, and .pdf files
knit("L2.Rmd")
markdownToHTML('L2.md', 'L2.html', options=c("use_xhml"))
system("pandoc -s L2.html -o L2.pdf")
system("pandoc -s --self L2.html -o L2_final.html")