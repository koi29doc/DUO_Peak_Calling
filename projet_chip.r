#################
### Chip project
#################

setwd("C:/Users/omics/Desktop/Data/Session 8")

region_a <- read.table(file = "ChIPseq_results_regionA.txt", header = T)
region_b <- read.table(file = "ChIPseq_results_regionB.txt", header = T)
region_c <- read.table(file = "ChIPseq_results_regionC.txt", header = T)

###
# Region A
###

View(region_a)
summary(region_a$IP > region_a$Control) # il faut tendre vers 94


### Visuisalisation des données
plot(x = region_a$Pos, y = region_a$IP + region_a$Control)
lines(x = region_a$Pos,y = region_a$IP, col = "red" )
lines(x = region_a$Pos,y = region_a$Control, col = "blue" )

### Création de la fenètre glissante

## Paramètres
taille_fenetre <- 10 #taille_fenetre désignle la taille de la fene^tre
taille_pas <-  6 # taille_pas désigne la valeur du pas
moy_ratio <- 1 # moy_ratio correspond au rapport des moyennes IP / CTRL

## Séquence de répétition 
seq_pas <- seq(from = 1, to = nrow(region_a),by = taille_pas)

## Function

sub_region_a_1 <- c()
sub_region_a_2 <- c()
sub_region_a_3 <- c()
sub_region_a_4 <- c()
sub_region_a_5 <- c()



## Version stringente
for (i in seq_pas){
  sub_region_a_1 <- region_a[i:(i+taille_fenetre), ]
  sub_region_a_2 <- region_a[i+taille_fenetre:(i+2*taille_fenetre), ]
  if ((mean(sub_region_a_1$IP, na.rm = T) > moy_ratio*mean(sub_region_a_1$Control, na.rm = T)) & (mean(sub_region_a_2$IP, na.rm = T) > moy_ratio*mean(sub_region_a_2$Control, na.rm = T))){
    sub_region_a_3 <- rbind(sub_region_a_3, sub_region_a_1)
    sub_region_a_4 <- sub_region_a_3[sub_region_a_3$IP>sub_region_a_3$Control,]
    rownames(sub_region_a_4) <- c()
    sub_region_a_5 <- sub_region_a_4[-which(duplicated(sub_region_a_4)), ]
  }
}

dim(sub_region_a_5)
sub_region_a_5

# Commentaires : la détectiono du pic et de sa taille dépens principalement du mean raio et varie peu selon le pas et les fenetres


## Version cool
for (i in seq_pas){
  sub_region_a_1 <- region_a[i:(i+taille_fenetre), ]
  if (mean(sub_region_a_1$IP, na.rm = T) > moy_ratio*mean(sub_region_a_1$Control, na.rm = T)){
    sub_region_a_2 <- rbind(sub_region_a_2, sub_region_a_1)
    sub_region_a_3 <- sub_region_a_2[sub_region_a_2$IP>sub_region_a_2$Control,]
    rownames(sub_region_a_3) <- c()
  }
}

head(sub_region_a_3)
dim(sub_region_a_3)

View(sub_region_a_3)

###
# Region B
###

View(region_b)
summary(region_b$IP > region_b$Control)


### Visuisalisation des données
plot(x = region_b$Pos, y = region_b$IP + region_b$Control)
lines(x = region_b$Pos,y = region_b$IP, col = "red" )
lines(x = region_b$Pos,y = region_b$Control, col = "blue" )

## Paramètres
taille_fenetre <- 4 #taille_fenetre désignle la taille de la fene^tre
taille_pas <-  3 # taille_pas désigne la valeur du pas
moy_ratio <- 1 # moy_ratio correspond au rapport des moyennes IP / CTRL

## Séquence de répétition 
seq_pas <- seq(from = 1, to = nrow(region_b),by = taille_pas)

## Function

sub_region_b_1 <- c()
sub_region_b_2 <- c()
sub_region_b_3 <- c()
sub_region_b_4 <- c()
sub_region_b_5 <- c()


## Version stringente
for (i in seq_pas){
  sub_region_b_1 <- region_b[i:(i+taille_fenetre), ]
  sub_region_b_2 <- region_b[i+taille_fenetre:(i+2*taille_fenetre), ]
  if ((mean(sub_region_b_1$IP, na.rm = T) > moy_ratio*mean(sub_region_b_1$Control, na.rm = T)) 
      & (mean(sub_region_b_2$IP, na.rm = T) > moy_ratio*mean(sub_region_b_2$Control, na.rm = T))){
    sub_region_b_3 <- rbind(sub_region_b_3, sub_region_b_1)
    sub_region_b_4 <- sub_region_b_3[sub_region_b_3$IP>sub_region_b_3$Control,]
    rownames(sub_region_b_4) <- c()
    sub_region_b_5 <- sub_region_b_4[-which(duplicated(sub_region_b_4)), ]
  }
}

dim(sub_region_b_5)
# View(sub_region_b_5)
barplot(region_b$Pos %in% sub_region_b_5$Pos)


