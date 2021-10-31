################################################################################
#                                                                              #
# Este programa calcula la informacion mutua condicional de todo el cromosoma  #
# 8 con todo el cromosoma 8 para valores de expresion de la                    #
# la variante de cancer luminalB condicionada con la variante del numero de    #
# copias de la region q24-3 perteneciente al mismo cromosoma.                  #
# Se guarda el resultado en una matriz externa compatible con R.               #
#                                                                              #
# Hecho por: Candelario Hernandez Gomez                                        #
# Ultima modificacion: 31 de octubre de 2021                                   # 
################################################################################

# Se leen las matrices de un archivo externo.

x <- read.table(file = "archivoOrdenadoControl.tsv", header = FALSE, sep = "\t" )
     # El archivo recibido debe estar ordenado de acuerdo al inicio de cada gen.
Rc8 = nrow(x)                       # Numero de renglones.
Cc8 = ncol(x)

rownames(x) <- x[,1]                # Nombre de renglones archivo de cromosoma 8.

print("Renglones cromosoma 8:")
print(Rc8)
print("Columnas c8")
print(Cc8)

# Se leen las variantes de numero de copias
#cnvs-q24-3-OrdenadosPorCercania.tsv

y <- read.table(file = "pseudoCNVS.tsv", header = FALSE, sep = " " )
Rrq = nrow(y)                       # Numero de renglones region q24_3.
Crq = ncol(y)

# Se crea el arreglo tridimensional que guardara los valores de IMC
v1 <- replicate(Rc8,0)
v2 <- replicate(Rc8,0)
v3 <- replicate(Rrq,0)

IMC <- array(c(v1,v2,v3), dim = c(Rc8,Rc8,Rrq))

###############################################################################
# Se llena la matriz.                                                         #
###############################################################################

print("Ya voy a llenar la matriz.")

library(infotheo)
for (i in 1:Rc8)                      # R cuenta desde 1
{                                     # Abre primer for

        fila1 = unname(unlist(x[i,]))
        fila1d <- discretize(fila1)   # Infotheo trabaja con numeros enteros.
        j <- i+1

        while( j <= Rc8)
        { # Abre while
                fila2 = unname(unlist(x[j,]))
                fila2d <- discretize(fila2)

                for(k in 1:1)
                {   # Abre tercer ciclo
                        fila3 = unname(unlist(y[k,]))
                        fila3d <- discretize(fila3)
                        imc  = condinformation(fila1d, fila2d, fila3d, method="emp")
                        IMC[i, j, k] = imc
                        IMC[j, i, k] = imc # IMC(i,j|k) == IMC(j,k|i)
                }   # Cierra tercer ciclo
        j <- j+1
       } #Cierra while
}                                   # Cierra for.




rownames(IMC) = c(rownames(x))
colnames(IMC) = c(rownames(x))

################################################################################
# Se guarda el archivo.                                                        #
################################################################################

saveRDS(IMC, file = "IMC-Control.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)






