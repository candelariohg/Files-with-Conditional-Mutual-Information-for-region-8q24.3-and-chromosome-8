################################################################################
#                                                                              #
# This program calculates the conditional mutual information of the entire     # 
# chromosome 8 with the entire chromosome 8 for expression values of the       # 
# luminalB cancer variant conditioned with the copy number variant of the q24-3#
# region belonging to the same chromosome.                                     #
# The result is saved in an external array compatible with R.                  #
#                                                                              #
# Made by: Candelario Hernandez Gomez                                          #
# Last modified: February 21, 2022                                             #
#                                                                              #
################################################################################


# Arrays are read from an external file.

x <- read.table(file = "Control.tsv", header = FALSE, sep = "\t" )
                             # Control.tsv is the expression array. 

Rc8 = nrow(x)                       # Number of rows in expression file. 
Cc8 = ncol(x)                       # Number of columns in expression file.

rownames(x) <- x[,1]                # Name of rows in expression file.

print("Rows: ")
print(Rc8)
print("Columns: ")
print(Cc8)

################################################################################
# The CNVs are read.                                                           # 
################################################################################

y <- read.table(file = "CNVS.tsv", header = FALSE, sep = " " )
Rrq = nrow(y)                       # Number of rows in region q24.3
Crq = ncol(y)                       # Numbero columns in region q24.3

################################################################################
# The three-dimensional array that will store the Mutual Conditional           # 
# Information is created.                                                      #
################################################################################

v1 <- replicate(Rc8,0)
v2 <- replicate(Rc8,0)
v3 <- replicate(Rrq,0)

IMC <- array(c(v1,v2,v3), dim = c(Rc8,Rc8,Rrq))

################################################################################
# The matrix is filled.                                                        # 
################################################################################

print("Filling the matrix.")

library(infotheo)
for (i in 1:Rc8)                      # R counts from 1
	
{                                     # Open the firs for

	fila1 = unname(unlist(x[i,]))
	fila1d <- discretize(fila1)   # Infotheo calculates the MIC 
	j <- i+1

        while( j <= Rc8)
	{ # Open while
                fila2 = unname(unlist(x[j,]))
		fila2d <- discretize(fila2)

		for(k in 1:1)
		{   # Open third cicle
			fila3 = unname(unlist(y[k,]))
			fila3d <- discretize(fila3)
			imc  = condinformation(fila1d, fila2d, fila3d, method="emp")
		        IMC[i, j, k] = imc	
			IMC[j, i, k] = imc # IMC(i,j|k) == IMC(j,k|i)
		}   # Close third cicle
	j <- j+1
       } #Close while
}                                   # Close for.




rownames(IMC) = c(rownames(x))
colnames(IMC) = c(rownames(x))


################################################################################
# File is saved                                                                #
################################################################################

saveRDS(IMC, file = "IMC-Control.Rcompatible", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)


