if (!require(Rmisc)) { install.packages("Rmisc") }
if (!require(prospectr)) { install.packages('prospectr')}
if (!require(tidyr)) { install.packages('tidyr')}
if (!require(readxl)) { install.packages('readxl')}
if (!require(dplyr)) { install.packages('dplyr')}
if (!require(purrr)) { install.packages('purrr')}
if (!require(ggplot2)) { install.packages('ggplot2')}


library(tidyr)
library(readxl)
library(Rmisc) # pour la commande summarySE
# Attention Rmisc charge aussi plyr qui entre en conflit avec dplyr
#You have loaded plyr after dplyr - this is likely to cause problems.
#If you need functions from both plyr and dplyr, please load plyr first, then dplyr:library(plyr); library(dplyr) 
#dplyr:library(plyr) # ca marche pas
library(dplyr)
library(prospectr)
library(purrr)
library(ggplot2)
library(forcats)


# Définir des variables pour l'ensemble du script

## Le répertoire de travail
WORKDIR <- "Path_to_your_data"

## Nombre de point qui vont être retiré par la smoothing
# Cette valeur doit toujours être paire (elle est égale au parametre w de la fonction savitzkyGolay -1)
MovingPoint <- 10

## Lignée de référence
RefLine <- "wt" 

## nb de réplicats techniques
RepTechNb <- 6

#===================================================================================================================================
# Fonctions pour le script
#===================================================================================================================================
# Charge les données à partir du csv
#la fonction ne charge que les lignes qui commencent par le motif suivant :
#un ou plusieurs chiffres suivi d'un point puis de un ou plusieurs chiffres puis d'une virgule 
#ou d'un ou plusiuers chiffre suivi d'une virgule, 
#suivi par des blocs de un ou plusieurs chiffres suivi d'un point suivi de un ou plusieurs chiffres suivi d'une virgule
# Retire les colonnes qui portent le même nom
# Retire la colonne sur numéraire 8
# Ne conserve que 2 chiffre après la virgule pour les valeurs numérique
# Retourne un data frame contenant les données
#===================================================================================================================================
f_load_data_csv <- function(csvFile, colName) {
  lines <- readLines(csvFile)
  linesNbToImport <- grep("^\\d+(\\.\\d+)?,\\d+\\.\\d+,", lines, perl = TRUE)
  df <- read.table(text = lines[linesNbToImport], header = FALSE, sep = ",", dec = ".", check.names = FALSE)
  colnames(df) <- c(colName)
  duplicated_names <- duplicated(colnames(df))
  df <- df[!duplicated_names]
  df[,8] <- NULL
  df <- df %>% mutate(across(where(is.numeric), ~ round(., 2)))
  return(df)
}
#===================================================================================================================================

#===================================================================================================================================
# Prépare les données
# Smoothe les données selon la méthode de Savitzky Golay et tiens compte de la supression des points au début et à la fin
# Identifie le maximun pour chaque lignée dans la zone 650 - 700 nm et normalise les données par cette valeur
# Fait la moyenne en ligne des réplicats techniques (cad par groupe de n colonnes)
# Retourne une liste des data frame contenant les étapes de préparatoin des données : smoothées, normalisées et moyenne des réplicats
# techniques
#===================================================================================================================================
f_process_data <- function(DATA, LINE){
  # n point moving average filter supress the (n-1)/2 first points and the (n-1)/2 last point. It is necessary 
  # to remove those points in the vector of the waves length
  lambda <- DATA$nm[-(1:MovingPoint/2)]
  lambda <- head(lambda, - MovingPoint/2)
  
  df_smooth <- savitzkyGolay( t(as.matrix(DATA[,-1])), p = 3, w = MovingPoint + 1, m = 0)
  # Function savitzkyGolay() works on matrix with samples in row and waves length in column, so use t() to convert the data to proper format
  df_smooth <- merge(t(df_smooth), lambda, by.x='row.names',by.y='row.names', all.x=T)
  df_smooth <- df_smooth %>% select(y, everything())
  df_smooth$Row.names <- NULL
  df_smooth <- df_smooth[order(df_smooth$y),]
  
  ## check result
  plot(DATA$nm,DATA[,2],"l")
  plot(df_smooth$y,df_smooth[,2],"l")
  
  # Normalisation
  # Pour normaliser il faut diviser l'ensemble des données par la hauteur du premier pic
  # On va rechercher ce pic dans l'interval 650 nm - 700 nm
  
  max_vec <- apply(df_smooth[(df_smooth$y>=650) & (df_smooth$y <700),][,-1], 2, max)
  
  df_norm <- as.data.frame(mapply('/', df_smooth[, -1], max_vec))
  
  plot(df_smooth$y,df_norm[,1],"l")
  
  # Moyenne des réplicats techniques
  # Calcul du nombre de replicat biologique
  nb_rep_bio <- ncol(df_norm)/RepTechNb
  select_col <- c(1:RepTechNb)
  
  df_mean <- data.frame(lambda)
  
  f_compute_repTech_mean <- function (select_col) {
    result <- as.data.frame(rowMeans(df_norm[,select_col]))
    return(cbind(df_mean, result))
  }
  
  for (repBioCount in c(1:nb_rep_bio)) {
    df_mean <- f_compute_repTech_mean(select_col)
    colnames(df_mean)[repBioCount +1] <- sprintf("%s_%d ", LINE, repBioCount)
    select_col <- select_col+RepTechNb
  }
  return(list(SMOOTH=df_smooth, NORM=df_norm, MEAN=df_mean))
}
#===================================================================================================================================


# Coprs du script
## Charger les données

# Charger les données pour la première lignée étudiée
data1 <- f_load_data_csv("g8_1.csv", c("nm", "g8_1_1", "g8_1_2", "g8_1_3", "g8_1_4", "g8_1_5", "g8_1_6"))
data2 <- f_load_data_csv("g8_2.csv", c("nm", "g8_2_1", "g8_2_2", "g8_2_3", "g8_2_4", "g8_2_5", "g8_2_6"))
data3 <- f_load_data_csv("g8_3.csv", c("nm", "g8_3_1", "nm", "g8_3_2", "nm", "g8_3_3", "nm", "g8_3_4", "nm", "g8_3_5", "nm","g8_3_6"))
data4 <- f_load_data_csv("g8_4.csv", c("nm", "g8_4_1", "nm", "g8_4_2", "nm", "g8_4_3", "nm", "g8_4_4", "nm", "g8_4_5", "nm","g8_4_6"))

#Assembler les données en un seul data frame
# RAssemble les data frame dans une liste pour les assembler en une seule commande
liste <- list(data1, data2, data3, data4)
# Crée une nouvelle liste vide pour contenir tous les data frame de toutes les lignées utilisées
Liste_data <- list()
# Ajoute à la liste un élément nommé par le nom de lapremière lignée et contenant 
#toutes les données assemblées en un seul data frame 
Liste_data[['g8']] <-purrr::reduce(.x = liste, merge, by = c('nm'), all = T)
# Si necessaire rempli les cases contenant un NA par la case la plus proche vers le haut
Liste_data[['g8']] <- fill(Liste_data[['g8']], !nm)


# Charger les données pour la seconde lignée étudiée
data1 <- f_load_data_csv("wt1.csv", c("nm", "wt_1_1", "nm", "wt_1_2", "nm", "wt_1_3", "nm", "wt_1_4", "nm", "wt_1_5", "nm", "wt_1_6"))
data2 <- f_load_data_csv("wt2.csv", c("nm", "wt_2_1", "nm", "wt_2_2", "nm", "wt_2_3", "nm", "wt_2_4", "nm", "wt_2_5", "nm", "wt_2_6"))
data3 <- f_load_data_csv("wt3.csv", c("nm", "wt_3_1", "nm", "wt_3_2", "nm", "wt_3_3", "nm", "wt_3_4", "nm", "wt_3_5", "nm", "wt_3_6"))
data4 <- f_load_data_csv("wt4.csv", c("nm", "wt_4_1", "nm", "wt_4_2", "nm", "wt_4_3", "nm", "wt_4_4", "nm", "wt_4_5", "nm", "wt_4_6"))

#Assembler les données en un seul data frame
liste <- list(data1, data2, data3, data4)
Liste_data[['wt']] <-purrr::reduce(.x = liste, merge, by = c('nm'), all = T)
Liste_data[['wt']] <- fill(Liste_data[['wt']], !nm)

# Charger les données pour la troisième lignée étudiée
data1 <- f_load_data_csv("h1_1.csv", c("nm", "h_1_1", "nm", "h_1_2", "nm", "h_1_3", "nm", "h_1_4", "nm", "h_1_5", "nm", "h_1_6"))
data2 <- f_load_data_csv("h1_2.csv", c("nm", "h_2_1", "nm", "h_2_2", "nm", "h_2_3", "nm", "h_2_4", "nm", "h_2_5", "nm", "h_2_6"))
data3 <- f_load_data_csv("h1_3.csv", c("nm", "h_3_1", "nm", "h_3_2", "nm", "h_3_3", "nm", "h_3_4", "nm", "h_3_5", "nm", "h_3_6"))
data4 <- f_load_data_csv("h1_4.csv", c("nm", "h_4_1", "nm", "h_4_2", "nm", "h_4_3", "nm", "h_4_4", "nm", "h_4_5", "nm", "h_4_6"))

#Assembler les données en un seul data frame
liste <- list(data1, data2, data3, data4)
Liste_data[['h']] <-purrr::reduce(.x = liste, merge, by = c('nm'), all = T)
Liste_data[['h']] <- fill(Liste_data[['h']], !nm)

# Re faire le même processus pour chaque lignée à étudier
#===================================================================================================================================

# Crée une liste vide qui contiendra les résulats du processing des données
resultat_final <- list()

# Exécute la fonction pour chaque lignée et ajoute le résultat dans la liste resultat_final en le nommant
for ( i in c(1:length(Liste_data)) ) {
  current_line <- names(Liste_data)[i]
  data <- Liste_data[[i]]
  resultat_final[[current_line]] <- f_process_data(DATA= data, LINE= current_line)
}


## Assembler les df_mean de toutes les lignées
# Rassemble les df dans une liste
Lmean <- list()
for ( i in c(1:length(resultat_final)) ) {
  df <- resultat_final[[i]][['MEAN']]
  Lmean <- c(Lmean, list(df))
}

#Joindre les éléments de la liste
complete_mean <-purrr::reduce(.x = Lmean, merge, by = c('lambda'), all = T)

##Calculer l'intervalle de confiance
# Rendre les données tidy : version longue
complete_mean_long <- pivot_longer(complete_mean,!lambda , names_to = "Line", values_to = "value")

# Séparer la colonne Line en Line et rep
complete_mean_long <- tidyr::separate(data = complete_mean_long, col = Line, c("Line", "rep"), sep = "_", remove=TRUE)


# Calculer l'intervalle de confiance
my_summary <- summarySE(complete_mean_long, measurevar="value", groupvars=c("lambda", "Line"))

# Créer un objet target_order pour forcer l'ordre dans les graph 
my_summary$Line <- as.factor(my_summary$Line)
L <- levels(my_summary$Line)
L <- L[!(L %in% RefLine)]
target_order <- c(RefLine, L)


## Faire les graph

# With 95% ci, all singing all dancing  FULL 650-800
p<- my_summary %>% mutate(Line = fct_relevel(Line, target_order)) %>%  
  ggplot(aes(x=lambda, y=value ))+
  geom_ribbon(aes(ymin=value - ci, ymax=value + ci, fill=Line),alpha=0.3)+
  scale_fill_manual(values = c("lightgrey", "lightgreen", "lightblue"))+
  geom_line(aes(x=lambda, y=value, color=Line), size=0.3 )+
  scale_color_manual(values = c("black", "green", "blue"))+
  scale_x_continuous(breaks=seq(650,800,100), limits=c(650,800))+
  theme_classic()

print(p)

ggsave("full77K.pdf",width = 8, height = 3)
ggsave("full77K.svg", width=8, height=3)

# With 95% ci, all singing all dancing ZOOMED 660-710
p<-my_summary %>% mutate(Line = fct_relevel(Line, target_order)) %>%
  ggplot(aes(x=lambda, y=value))+
  geom_ribbon(aes(ymin=value - ci, ymax=value + ci, fill=Line),alpha=0.3)+
  scale_fill_manual(values = c("lightgrey", "lightgreen", "lightblue"))+
  geom_line(aes(x=lambda, y=value, color=Line), size=0.3 )+
  scale_color_manual(values = c("black", "green", "blue"))+
  #scale_color_viridis_d(begin = 0, end = 1, option = 'viridis')+
  scale_x_continuous(breaks=seq(660,710,20), limits=c(660,710))+
  scale_y_continuous(limits = c(0, max(my_summary[660:710, 4])))+
  theme_classic()

print(p)

ggsave("zoom77K.pdf",width = 8, height = 3)
ggsave("zoom77K.svg", width=8, height=3)



































