# AUTHOR: Jesús Antonio Rocamontes Morales
# E-MAIL: jesus.rocamontes@uaz.edu.mx
# INSTITUTION: El Colegio de la Frontera Sur / Universidad Autonoma de Zacatecas
# SCRIPT DESCRIPTION: Un script en R para calcular proporción de nucleotidos en secuencias. 
# a simple R script for calculating nucleotide composition in sequences. 

# Instalando Biostrings, que es parte del paquete BiocManager, el
# cual no es parte del CRAN. Usar este script:
# Install the following packages for the script.

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


BiocManager::install("Biostrings")

library(Biostrings)
library(openxlsx)

# Cargamos las secuencias con la función readDNAStringSet.
# No es necesario que es estas secuencias esten alineadas unas con otras, pues
# estamos calculando los valores independientemente para cada una para colocarlos
# en una tabla. 

# Load the sequences and create an object for each with readDNAStringSet
# Sequences do not have to be aligned, values are calculated for each separately.

# Incluye el número de secuencias que quieras analizar.
# Include any number of sequences. 

Sequence1 <- readDNAStringSet("/Your_file.fas")
Sequence2 <- readDNAStringSet("/Your_file.fas")
Sequence3 <- readDNAStringSet("/Your_file.fas")

# Creamos un set con todas las secuencias.
# Create a larger set with all the sequences.

Sequences <- DNAStringSet(c(morenoi, mutica_c, mutica_I))

# Calculamos la composición de nucleotidos. Aqui buscamos la frecuencia de los
#nucleotidos que indiquemos.Se cuenta la frecuencia de cada nucleotido.
# Nucleotide composition is calculated here.We calculate frequency for reach
# letter.
base_composition <- letterFrequency(Sequences, letters = c("A", "T", "C", "G"))

# Calculamos el total de nucleotidos para cada secuencia.
# Here we calculate the total for each nucleotide.
total_nucleotides <- rowSums(base_composition)

# Porcentajes para cada nucleotido. Multiplicamos por 100 para obtener porcentajes.
# Percentage for reach nucleotide is calculated here. 

A_percentage <- base_composition[,"A"] / total_nucleotides * 100
T_percentage <- base_composition[,"T"] / total_nucleotides * 100
C_percentage <- base_composition[,"C"] / total_nucleotides * 100
G_percentage <- base_composition[,"G"] / total_nucleotides * 100
AT_percentage <- (base_composition[,"A"] + base_composition[,"T"]) / total_nucleotides * 100
GC_percentage <- (base_composition[,"G"] + base_composition[,"C"]) / total_nucleotides * 100

# Curvas para AT y GC. Esta formula es de Grigoriev, 1998 - Analyzing genomes with cumulative skew diagrams
# AT and GC skew are calculated here. The formulae is from Grigoriev 1998.

AT_skew <- (base_composition[,"A"] - base_composition[,"T"]) / (base_composition[,"A"] + base_composition[,"T"])
GC_skew <- (base_composition[,"G"] - base_composition[,"C"]) / (base_composition[,"G"] + base_composition[,"C"])

# Convertimos los resultados en un data frame.
# Convert your results into a data frame.

results <- data.frame(Sequences = names(Sequences), 
                      A = A_percentage,
                      T = T_percentage,
                      C = C_percentage,
                      G = G_percentage,
                      AT = AT_percentage,
                      GC = GC_percentage,
                      AT_Skew = AT_skew,
                      GC_Skew = GC_skew)

# Creamos la tabla, así bien perrona. No olvides cambiar el destino donde 
# se guardará tu tabla.La puedes editar en LibreCalc, Excel, etc.
# Create a table and save it in your prefered location. Edit with pertinent 
# software. 

print(results)

write.xlsx(results, "/home/documents/file.xlsx")
