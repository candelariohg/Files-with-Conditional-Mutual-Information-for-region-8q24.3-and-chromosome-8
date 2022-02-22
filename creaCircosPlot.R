################################################################################
#                                                                              #
# This program opens a file with N layers of MxM filled with CMI values for    #
# different genes, receives a number m from the user to create an adjacency    #
# matrix and a circos plot from it.                                            #
#                                                                              #
# You should run it (in bash) in this way:                                     #
#                                                                              #
# $ Rscript thisFile.R 300                                                     #
#                                                                              #
# where 300 is the number of links that you want in the network thah will be   #
# created.                                                                     #
#                                                                              #
# Made by: Candelario Hernandez Gomez                                          #
# Last modified: February 22, 2022                                             #
################################################################################


library(circlize)            # necessary to plot the circos.

################################################################################
# 1) R compatible tsv file is received.                                        #
################################################################################

valoresIMC <- readRDS("IMC-LumA.Rcompatible", refhook = NULL)

Tamano = dim(valoresIMC)[1]  # Dimension of layer.

print("The first dimension:")
print(Tamano)

dim2 = dim(valoresIMC)[2]
print("The second dimension:")
print(dim2)

dim3 = dim(valoresIMC)[3]
print("The third dimension: ")
print(dim3)

################################################################################
# 2) The file with Ensembl identifiers and regions is received.                #
################################################################################

nombres <- read.table(file = '447-IE-Region.txt', sep = ' ', header = FALSE)

                             # 447-IE-Region.txt is a file with two columns
                             # containing the gene IDs and the region to which
                             # they belong. 


################################################################################
# 3) A loop is created that runs through all the values of the conditionals,   #
# the third dimension of the received array.                                   #
################################################################################

i = 1                        # Initial conditional layer 
f = 1                        # Final conditional layer
lapply(
i:f,

function(k)
{                            # Opens function.

s <- valoresIMC[,,k]                

################################################################################
# 4) In this for the names of the rows are changed from numbers to ENSEMBL     #
# identifiers.                                                                 #
################################################################################

for (i in 1:Tamano)
       row.names(s)[i] <- as.character(nombres[i,2])

colnames(s) = c(rownames(s))        



################################################################################
# 5) The adjacency matrix is created.                                          #
################################################################################

v1 <- replicate(Tamano,0)
v2 <- replicate(Tamano,0)

mAd <- array(c(v1,v2), dim = c(Tamano,Tamano))


################################################################################
#6) The layer s is converted into a vector, in order to sort it with the sort  #
# command.                                                                     #
################################################################################

arreglo <- as.vector(s)              
                                    
col = dim(s)[1]              # Number of columns.

arreglo_ordenado = sort(arreglo, decreasing = TRUE) # The array is ordered.

################################################################################
# 7) m is received from the console. The value m determines the number of pairs#
# that will be taken. The values are equal in pairs, because the matrix is     #
# symmetric.                                                                   #
#                                                                              #
# Execution should be like this:                                               #
# Rscript creaCircosPlot.R 300                                                 #
# where m = 300. A network is created with the first m links.                  # 
#                                                                              #
################################################################################

args <- commandArgs(TRUE)
m <- as.integer(args[1])

n = 2*m

limite = arreglo_ordenado[n]

print("The lower limit value of the CMI:")
print(limite)


################################################################################
# 8) Adjacency matrix is filled.                                               #
################################################################################
                   
                             # The matrix has zeros in all the places, they are
                             # turned into 1 if the interaction is higher than 
                             # limite.

mAd <- array(c(v1,v2), dim = c(Tamano,Tamano))
contador = 0
for (i in 1:Tamano)
{
        for (j in 1:Tamano)
        {
        if (s[i,j] >= limite)
        {
        mAd[i,j] = 1
        }

        }
}


################################################################################
# 9)  circos is drawn.                                                         #
################################################################################

################################################################################
# IMPORTANT: Change the first part of the name depending on the type of        #
# file: LumB or Control.                                                       # 
################################################################################

pdf(file= paste("circosLumB-", m,"capa-", k, ".pdf", sep=""))



################################################################################
# Bed is read.                                                                 #
################################################################################
                             # The bed ins needed to create the circos plot.

cama <- read.csv(file = '/home/cando/programasR/segundoArticulo/datitos.txt', sep =' ', header =FALSE)

print("Read the bed")




################################################################################
# tracks are created.                                                          #
################################################################################


circos.initializeWithIdeogram(chromosome.index = paste0("chr", c(8)),
                              plotType = c("ideogram", "axis", "labels"),
                              track.height = 0.2,
                              ideogram.height = convert_height(3, "mm"),
)

################################################################################
# labels are added.                                                            #
################################################################################

circos.genomicLabels(cama, labels.column = 4, side = "outside")

print("genomicLabels")

################################################################################
# links are added.                                                             #
################################################################################
contador = 0
for (i in 1:Tamano)
{ # Abre for

     j <- i + 1
while(j <= Tamano)
     { # opens while 
                   if (1 == mAd[i,j])
                  {                        # opnes if belonging to the link. 

	          vertice1 = cama[i,2]
                  vertice2 = cama[j,2]	     

	          if (vertice1 <= 45200000 && vertice2 <= 45200000)  # arm p
	          {
			  
	          if ((vertice1 >= 1) && (vertice2 <= 2300000) || (vertice1 >= 2300001 && vertice2 <= 6300000) ||
		      (vertice1 >= 6300001) && (vertice2 <= 12800000) || (vertice1 >= 12800001) && (vertice2 <= 19200000)  || 
		      (vertice1 >= 19200001) && (vertice2 <= 23500000) || (vertice1 >= 23500001) && (vertice2 <= 27500000) || 
		      (vertice1 >= 27500001) && (vertice2 <= 29000000) || (vertice1 >= 29000001) && (vertice2 <= 36700000) ||
		      (vertice1 >= 36700001) && (vertice2 <= 38500000) || (vertice1 >= 38500001) && (vertice2 <= 39900000) ||
		      (vertice1 >= 39900001) && (vertice2 <= 43200000) || (vertice1 >= 43200001) && (vertice2 <= 45200000) )
		  {   # color if start and end belong to the same citoband
                  circos.link("chr8", cama[i,2], "chr8", cama[j,2], col = "#13346E", lwd = 0.7 )        #arm P #
                  contador = contador + 1
		  }
             
	          else                                                   #intercitoband conections en arm p
	          {    # Not considered cases (in arm p but intercitoband).
                  circos.link("chr8", cama[i,2], "chr8", cama[j,2], col = "#70552B", lwd = 0.7 )        # ARM P #
                  contador = contador + 1
	          }    # Closes not considered cases (in arm p but intercitoband).

	          }                                                                                     # Brazo p


		  #########################################################################
		  # Now we consider that the points are not in arm p, but are both in q   #
		  # or 1 in p and the other in q.                                         #  
		  ######################################################################### 

	          else 
                  {    # Opens 'if not completely in arm p' 

		  if(vertice1 > 45200000 & vertice2 > 45200000)   
		  {                                                                 #  If both genes are in q. 

	          if ((vertice1 >= 45200001 && vertice2 <= 47200000) || (vertice1 >= 47200001 && vertice2 <= 51300000) ||
                      (vertice1 >= 51300001 && vertice2 <= 51700000) || (vertice1 >= 51700001 && vertice2 <= 54600000) ||
                      (vertice1 >= 54600001 && vertice2 <= 60600000) || (vertice1 >= 60600001 && vertice2 <= 61300000) ||
                      (vertice1 >= 61300001 && vertice2 <= 65100000) || (vertice1 >= 65100001 && vertice2 <= 67100000) ||
                      (vertice1 >= 67100001 && vertice2 <= 69600000) || (vertice1 >= 69600001 && vertice2 <= 72000000) ||
                      (vertice1 >= 72000001 && vertice2 < 74600000) || (vertice1 >= 74600001 && vertice2 <= 74700000) ||
                      (vertice1 >= 74700001 && vertice2 < 83500000) || (vertice1 >= 83500001 && vertice2 <= 85900000) ||
                      (vertice1 >= 85900001 && vertice2 < 92300000) || (vertice1 >= 92300001 && vertice2 <= 97900000) ||
                      (vertice1 >= 97900001 && vertice2 < 100500000) || (vertice1 >= 100500001 && vertice2 <= 105100000) ||
                      (vertice1 >= 105100001 && vertice2 < 109500000) || (vertice1 >= 109500001 && vertice2 <= 111100000) ||
                      (vertice1 >= 111100001 && vertice2 < 116700000) || (vertice1 >= 116700001 && vertice2 <= 118300000) ||
                      (vertice1 >= 118300001 && vertice2 < 121500000) || (vertice1 >= 121500001 && vertice2 <= 126300000) ||
                      (vertice1 >= 126300001 && vertice2 < 130400000) || (vertice1 >= 130400001 && vertice2 <= 138900000) || 
                      (vertice1 >= 138900001) && (vertice2 <= 145138636 ))
                      

		  {        # Opens 'If both genes are in one of the bands of arm q' 
                  
		  if((vertice1 >= 138900001) && (vertice2 <= 145138636))  
	          { # If both genes are in the region q24.3
                  circos.link("chr8", cama[i,2], "chr8", cama[j,2], col = "#CA8E17", lwd = 0.7 )        
                  contador = contador + 1
		  }
		  
		  else
		  {   # If both genes are in a q band but not in q24.3
			  
                  circos.link("chr8", cama[i,2], "chr8", cama[j,2], col = "#13346E", lwd = 0.7 )        
                  contador = contador + 1
		  }  # Closes If both genes are in q band but not in q24.3 

		  }        # Closes 'If both genes are in one of the bands of arm q' 


	          else 
		  {   # Opens 'If both genes are in q arm but in different bands' 
	          
	          if (vertice1 >= 138900001 || vertice2 >= 138900001)
                  {                                 # Opens 'If one gene is in region q24.3'
                  circos.link("chr8", cama[i,2], "chr8", cama[j,2], col = "#CA8E17", lwd = 0.7)
		  contador = contador + 1
		  }                                 # Closes 'If one gene is in region q24.3'

		  else
		  { # Opens 'If genes are in different q bands but no one is in region q24.3'
                  circos.link("chr8", cama[i,2], "chr8", cama[j,2], col = "#70552B", lwd = 0.7)
		  contador = contador + 1
		  } # Closes 'If genes are in different q bands but no one is in region q24.3'
		  }  # Closes 'If both genes are in q arm but in different bands' 


		  }  # Closes 'If completely in arm q ' 

		  

                  
		   
		  

		  ############################################################################
		  # If one gen is in one arm and the other in q                              #
		  ############################################################################

		  else                                    # Inter arm links 
	          {

	          if((vertice1 >= 138900001) || (vertice2 >= 138900001))
	          {    # Opens if one gene is in region q24.3  #
		  circos.link("chr8", cama[i,2], "chr8", cama[j,2], col="#CA8E17", lwd = 0.7)
                  contador = contador + 1
		  }    # Closes if one gene is in region 24.3
		  
		  else
		  {   #If one gene is in p, other in q but no one in region q24.3
		  circos.link("chr8", cama[i,2], "chr8", cama[j,2], col="#BF1A51", lwd = 0.7)
                  contador = contador + 1
		  }   #Closes If one gene is in p, other in q but no one in region q24.3 

		  }                                      # Inter arm links 


                  } # Closes 'If not completely in arm p' 
		  

                  }	               # Closes if corresponding to the link 
     j <- j + 1

     } # Cierra while 

} # Cierra for

print("Number of links: ")
print(contador)
circos.clear()

dev.off()                                        # Figure is finished
} # Closes function
)   # Closes lapply



