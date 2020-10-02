# displacement vectors for elements of interest

jaw.wire <- read.csv(file.path(fullpath, "DfA/jawlinks.csv"))


# SETUP 

# some packages
library(geomorph)
library(rstudioapi)

# set up the directory path for the session using ## rstudioapi::getSourceEditorContext()$path) ##
fullpath <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
datapath <- file.path(fullpath, "data")

# download data
skull.df <- read.csv(file = file.path(datapath, "skull_df.csv"), 
               row.names = 1, 
               na.strings = c("NA", ""))
vert.df <- read.csv(file = file.path(datapath, "vert_df_w_impute.csv"), 
                    na.strings = c("NA", ""))

##########################################################################################
##########################################################################################
##########################################################################################

# data subsets (elements)

# JAW
# columns
jaw.Rcols <- which(colnames(skull.df)=="J_R1_X"):which(colnames(skull.df)=="J_R20_Z")
jaw.Lcols <- which(colnames(skull.df)=="J_L1_X"):which(colnames(skull.df)=="J_L20_Z")

# rows (only for individuals with all 5 stages)
# using left side because it includes one more individual
jaw.pa <- rowSums(matrix(ifelse(!is.na(skull.df[,jaw.Lcols[1]]), 1, 0),
                           ncol=5, byrow=T))
names(jaw.pa) <- unique(skull.df$ID)
jaw.rows <- as.character(skull.df$ID) %in% names(jaw.pa)[jaw.pa==5]

# landmarks
jaw.lms <- skull.df[jaw.rows, jaw.Lcols]
# covariates
jaw.ages <- skull.df$Age[jaw.rows]
jaw.ID  <- droplevels(skull.df$ID[jaw.rows])
jaw.sex  <- skull.df$Sex[jaw.rows]
jaw.gh <-skull.df$GH[jaw.rows]

##########################################################################################

# FACE
# modules
module.df <- read.csv(file.path(datapath, "CranialLMGuide.csv"))
face.mod <- module.df[module.df$Module=="Face",]
base.mod <- module.df[module.df$Module=="Base",]
vault.mod <- module.df[module.df$Module=="Vault",]

cran.cols <- which(colnames(skull.df)=="MSI_L_X"):which(colnames(skull.df)=="PSA_R_Z")
cran.noms <- colnames(skull.df)[cran.cols]
cran.noms <- substr(cran.noms, 1, nchar(cran.noms)-2)

# vector logical for which columns apply to which modules
face.cols <- cran.noms %in% face.mod$LM_name
base.cols <- cran.noms %in% base.mod$LM_name
vault.cols <- cran.noms %in% vault.mod$LM_name

# lms by module
cran.pa <- rowSums(matrix(ifelse(!is.na(skull.df[,cran.cols[1]]), 1, 0),
                         ncol=5, byrow=T))
names(cran.pa) <- unique(skull.df$ID)
cran.rows <- as.character(skull.df$ID) %in% names(cran.pa)[cran.pa==5]

cran.lms <- skull.df[cran.rows, cran.cols]
face.lms <- cran.lms[,face.cols]
base.lms <- cran.lms[,base.cols]
vault.lms <- cran.lms[,vault.cols]

# covariates
cran.ages <- skull.df$Age[cran.rows]
cran.ID  <- droplevels(skull.df$ID[cran.rows])
cran.sex  <- skull.df$Sex[cran.rows]
cran.gh <-skull.df$GH[cran.rows]

##########################################################################################

# VERTEBRAE

# # columns
vert.cols <- which(colnames(vert.df)=="M1X"):which(colnames(vert.df)=="L14Z")

# rows (only for individuals with all 5 stages)
vert.pa <- names(which(table(vert.df$ID)==5))
vert.rows <- as.character(vert.df$ID) %in% vert.pa

# landmarks
vert.lms <- vert.df[vert.rows, vert.cols]
# covariates
vert.ages <- vert.df$Age[vert.rows]
vert.ID  <- droplevels(vert.df$ID[vert.rows])
vert.sex  <- vert.df$Sex[vert.rows]
vert.gh <-vert.df$GH[vert.rows]

##########################################################################################
##########################################################################################
##########################################################################################

# DATA ARRAYS

jaw.arr <- arrayspecs(A = jaw.lms, p = ncol(jaw.lms)/3, k = 3)
jaw.x <- colnames(jaw.lms)[seq(1,ncol(jaw.lms), 3)]
rownames(jaw.arr) <- substr(jaw.x, 1, nchar(jaw.x)-2)

face.arr <- arrayspecs(A = face.lms, p = ncol(face.lms)/3, k = 3)
face.x <- colnames(face.lms)[seq(1,ncol(face.lms), 3)]
rownames(face.arr) <- substr(face.x, 1, nchar(face.x)-2)

vault.arr <- arrayspecs(A = vault.lms, p = ncol(vault.lms)/3, k = 3)
vault.x <- colnames(vault.lms)[seq(1,ncol(vault.lms), 3)]
rownames(vault.arr) <- substr(vault.x, 1, nchar(vault.x)-2)

base.arr <- arrayspecs(A = base.lms, p = ncol(base.lms)/3, k = 3)
base.x <- colnames(base.lms)[seq(1,ncol(base.lms), 3)]
rownames(base.arr) <- substr(base.x, 1, nchar(base.x)-2)

vert.arr <- arrayspecs(A = vert.lms, p = ncol(vert.lms)/3, k = 3)
vert.x <- colnames(vert.lms)[seq(1,ncol(vert.lms), 3)]
rownames(vert.arr) <- substr(vert.x, 1, nchar(vert.x)-1)

##########################################################################################
##########################################################################################
##########################################################################################

# IDENTIFY MIDLINES
face.pairs <- cbind(grep("_L", rownames(face.arr)), grep("_R", rownames(face.arr)))
vault.pairs <- cbind(grep("_L", rownames(vault.arr)), grep("_R", rownames(vault.arr)))
base.pairs <- cbind(grep("_L", rownames(base.arr)), grep("_R", rownames(base.arr)))
vert.pairs <- cbind(grep("L", rownames(vert.arr)), grep("R", rownames(vert.arr)))

##########################################################################################
##########################################################################################
##########################################################################################

# SUPERIMPOSITIONS

# Jaw
jaw.gpa <- gpagen(jaw.arr)
jaw3D <- jaw.gpa$coords #Position of specimen after GPA
jaw2D <- two.d.array(jaw3D)
jaw.CS <- jaw.gpa$Csize

# Face
face.gpa <- gpagen(face.arr)
face3D <- face.gpa$coords #Position of specimen after GPA
face2D <- two.d.array(face3D)
face.CS <- face.gpa$Csize

# Vault
vault.gpa <- gpagen(vault.arr)
vault3D <- vault.gpa$coords #Position of specimen after GPA
vault2D <- two.d.array(vault3D)
vault.CS <- vault.gpa$Csize

# Base
base.gpa <- gpagen(base.arr)
base3D <- base.gpa$coords #Position of specimen after GPA
base2D <- two.d.array(base3D)
base.CS <- base.gpa$Csize

# Vertebrae
vert.gpa <- gpagen(vert.arr)
vert3D <- vert.gpa$coords #Position of specimen after GPA
vert2D <- two.d.array(vert3D)
vert.CS <- vert.gpa$Csize

#Calculate Group Means

JL.agg <- aggregate(two.d.array(coords),
                    by=list(jaw.gh, jaw.ages), FUN=mean)
JL.mus <- arrayspecs(JL.agg[,-(1:2)],
                     p=dim(coords)[[1]],k=dim(coords)[[2]])

dimnames(JL.mus)[[3]] <- c(paste(levels(jaw.gh), rep("21", times=4)), 
                           paste(levels(jaw.gh), rep("28", times=4)),
                           paste(levels(jaw.gh), rep("35", times=4)),
                           paste(levels(jaw.gh), rep("45", times=4)),
                           paste(levels(jaw.gh), rep("60", times=4)))

#Deficients
defmu21<- JL.mus[,,1]
defmu28 <- JL.mus[,,5]
defmu35 <- JL.mus[,,9]
defmu45 <- JL.mus[,,13]
defmu60 <- JL.mus[,,17]

#Sufficients
suffmu21<- JL.mus[,,4]
suffmu28 <- JL.mus[,,8]
suffmu35 <- JL.mus[,,12]
suffmu45 <- JL.mus[,,16]
suffmu60 <- JL.mus[,,20]

#Early
earlymu21<- JL.mus[,,2]
earlymu28<- JL.mus[,,6]
earlymu35<- JL.mus[,,10]
earlymu45<- JL.mus[,,14]
earlymu60<- JL.mus[,,18]

#Late
latemu21<- JL.mus[,,3]
latemu28<- JL.mus[,,7]
latemu35<- JL.mus[,,11]
latemu45<- JL.mus[,,15]
latemu60<- JL.mus[,,19]

# select colors for plotting
defmu.col <- "red"
suffmu.col <- "black"
earlymu.col <- "blue"
latemu.col <- "darkgreen"

#Deficient displacement vectors compared to sufficient (the magnifier is 4)
def21mag <- ((defmu21-suffmu21)*4)+suffmu21
def28mag <- ((defmu28-suffmu28)*4)+suffmu28
def35mag <- ((defmu35-suffmu35)*4)+suffmu35
def45mag <- ((defmu45-suffmu45)*4)+suffmu45
def60mag <- ((defmu60-suffmu60)*4)+suffmu60

#Early displacement vectors compared to sufficient (the magnifier is 4)
early21mag <- ((earlymu21-suffmu21)*4)+suffmu21
early28mag <- ((earlymu28-suffmu28)*4)+suffmu28
early35mag <- ((earlymu35-suffmu35)*4)+suffmu35
early45mag <- ((earlymu45-suffmu45)*4)+suffmu45
early60mag <- ((earlymu60-suffmu60)*4)+suffmu60

#Late displacement vectors compared to sufficient (the magnifier is 4)
late21mag <- ((latemu21-suffmu21)*4)+suffmu21
late28mag <- ((latemu28-suffmu28)*4)+suffmu28
late35mag <- ((latemu35-suffmu35)*4)+suffmu35
late45mag <- ((latemu45-suffmu45)*4)+suffmu45
late60mag <- ((latemu60-suffmu60)*4)+suffmu60

#######Defficient######

# plotting Deficient displacement over sufficient wireframe day 21 
open3d(windowRect=c(20,50,1000,700), box = F, asp=F) # L, T, R, B
spheres3d(suffmu21, col = "black", radius=0.003, alpha=1, main = "")
for(i in 1:nrow(wire)){
  suffmu21.2 <- rbind(suffmu21[wire[i,1],],
                      suffmu21[wire[i,2],])
  lines3d(suffmu21.2, col=suffmu.col, lwd=2)
}
for(i in 1:nrow(suffmu21)) arrow3d(suffmu21[i,], def21mag[i,], 
                                   type = "lines", col = "blue", lwd=3)

# plotting Deficient displacement over sufficient wireframe day 28 
open3d(windowRect=c(20,50,1000,700), box = F, asp=F) # L, T, R, B
spheres3d(suffmu28, col = "black", radius=0.003, alpha=1, main = "")
for(i in 1:nrow(wire)){
  suffmu28.2 <- rbind(suffmu28[wire[i,1],],
                      suffmu28[wire[i,2],])
  lines3d(suffmu28.2, col=suffmu.col, lwd=2)
}
for(i in 1:nrow(suffmu28)) arrow3d(suffmu28[i,], def28mag[i,], 
                                   type = "lines", col = "blue", lwd=3)

# plotting Deficient displacement over sufficient wireframe day 35 
open3d(windowRect=c(20,50,1000,700), box = F, asp=F) # L, T, R, B
spheres3d(suffmu35, col = "black", radius=0.003, alpha=1, main = "")
for(i in 1:nrow(wire)){
  suffmu35.2 <- rbind(suffmu35[wire[i,1],],
                      suffmu35[wire[i,2],])
  lines3d(suffmu35.2, col=suffmu.col, lwd=2)
}
for(i in 1:nrow(suffmu35)) arrow3d(suffmu35[i,], def35mag[i,], 
                                   type = "lines", col = "blue", lwd=3)

# plotting Deficient displacement over sufficient wireframe day 45 
open3d(windowRect=c(20,50,1000,700), box = F, asp=F) # L, T, R, B
spheres3d(suffmu45, col = "black", radius=0.003, alpha=1, main = "")
for(i in 1:nrow(wire)){
  suffmu45.2 <- rbind(suffmu45[wire[i,1],],
                      suffmu45[wire[i,2],])
  lines3d(suffmu45.2, col=suffmu.col, lwd=2)
}
for(i in 1:nrow(suffmu45)) arrow3d(suffmu45[i,], def45mag[i,], 
                                   type = "lines", col = "blue", lwd=3)

# plotting Deficient displacement over sufficient wireframe day 60 
  open3d(windowRect=c(20,50,1000,700), box = F, asp=F) # L, T, R, B
  spheres3d(suffmu60, col = "black", radius=0.003, alpha=1, main = "")
  for(i in 1:nrow(wire)){
    suffmu60.2 <- rbind(suffmu60[wire[i,1],],
                        suffmu60[wire[i,2],])
    lines3d(suffmu60.2, col=suffmu.col, lwd=2)
  }
  for(i in 1:nrow(suffmu60)) arrow3d(suffmu60[i,], def60mag[i,], 
                                     type = "lines", col = "blue", lwd=3)

######Early##########
# plotting early displacement over sufficient wireframe day 21
open3d(windowRect=c(20,50,1000,700), box = F, asp=F) # L, T, R, B
spheres3d(suffmu21, col = "black", radius=0.003, alpha=1, main = "")
for(i in 1:nrow(wire)){
  suffmu21.2 <- rbind(suffmu21[wire[i,1],],
                      suffmu21[wire[i,2],])
  lines3d(suffmu21.2, col=suffmu.col, lwd=2)
}
for(i in 1:nrow(suffmu21)) arrow3d(suffmu21[i,], early21mag[i,], 
                                   type = "lines", col = "red", lwd=3)

# plotting early displacement over sufficient wireframe day 28
open3d(windowRect=c(20,50,1000,700), box = F, asp=F) # L, T, R, B
spheres3d(suffmu28, col = "black", radius=0.003, alpha=1, main = "")
for(i in 1:nrow(wire)){
  suffmu28.2 <- rbind(suffmu28[wire[i,1],],
                      suffmu28[wire[i,2],])
  lines3d(suffmu28.2, col=suffmu.col, lwd=2)
}
for(i in 1:nrow(suffmu28)) arrow3d(suffmu28[i,], early28mag[i,], 
                                   type = "lines", col = "red", lwd=3)

# plotting early displacement over sufficient wireframe day 35
open3d(windowRect=c(20,50,1000,700), box = F, asp=F) # L, T, R, B
spheres3d(suffmu35, col = "black", radius=0.003, alpha=1, main = "")
for(i in 1:nrow(wire)){
  suffmu35.2 <- rbind(suffmu35[wire[i,1],],
                      suffmu35[wire[i,2],])
  lines3d(suffmu35.2, col=suffmu.col, lwd=2)
}
for(i in 1:nrow(suffmu35)) arrow3d(suffmu35[i,], early35mag[i,], 
                                   type = "lines", col = "red", lwd=3)

# plotting early displacement over sufficient wireframe day 45
open3d(windowRect=c(20,50,1000,700), box = F, asp=F) # L, T, R, B
spheres3d(suffmu45, col = "black", radius=0.003, alpha=1, main = "")
for(i in 1:nrow(wire)){
  suffmu45.2 <- rbind(suffmu45[wire[i,1],],
                      suffmu45[wire[i,2],])
  lines3d(suffmu45.2, col=suffmu.col, lwd=2)
}
for(i in 1:nrow(suffmu45)) arrow3d(suffmu45[i,], early45mag[i,], 
                                   type = "lines", col = "red", lwd=3)

# plotting early displacement over sufficient wireframe day 60
  open3d(windowRect=c(20,50,1000,700), box = F, asp=F) # L, T, R, B
  spheres3d(suffmu60, col = "black", radius=0.003, alpha=1, main = "")
  for(i in 1:nrow(wire)){
    suffmu60.2 <- rbind(suffmu60[wire[i,1],],
                        suffmu60[wire[i,2],])
    lines3d(suffmu60.2, col=suffmu.col, lwd=2)
  }
  for(i in 1:nrow(suffmu60)) arrow3d(suffmu60[i,], early60mag[i,], 
                                     type = "lines", col = "red", lwd=3)

#####Late#####
# plotting late displacement over sufficient wireframe day 21
open3d(windowRect=c(20,50,1000,700), box = F, asp=F) # L, T, R, B 
spheres3d(suffmu21, col = "black", radius=0.003, alpha=1, main = "")
for(i in 1:nrow(wire)){
  suffmu21.2 <- rbind(suffmu21[wire[i,1],],
                      suffmu21[wire[i,2],])
  lines3d(suffmu21.2, col=suffmu.col, lwd=2)
}
for(i in 1:nrow(suffmu21)) arrow3d(suffmu21[i,], late21mag[i,], 
                                   type = "lines", col = "chartreuse3", lwd=3)

# plotting late displacement over sufficient wireframe day 28
open3d(windowRect=c(20,50,1000,700), box = F, asp=F) # L, T, R, B
spheres3d(suffmu28, col = "black", radius=0.003, alpha=1, main = "")
for(i in 1:nrow(wire)){
  suffmu28.2 <- rbind(suffmu28[wire[i,1],],
                      suffmu28[wire[i,2],])
  lines3d(suffmu28.2, col=suffmu.col, lwd=2)
}
for(i in 1:nrow(suffmu28)) arrow3d(suffmu28[i,], late28mag[i,], 
                                   type = "lines", col = "chartreuse3", lwd=3)

# plotting late displacement over sufficient wireframe day 35
open3d(windowRect=c(20,50,1000,700), box = F, asp=F) # L, T, R, B
spheres3d(suffmu35, col = "black", radius=0.003, alpha=1, main = "")
for(i in 1:nrow(wire)){
  suffmu35.2 <- rbind(suffmu35[wire[i,1],],
                      suffmu35[wire[i,2],])
  lines3d(suffmu35.2, col=suffmu.col, lwd=2)
}
for(i in 1:nrow(suffmu35)) arrow3d(suffmu35[i,], late35mag[i,], 
                                   type = "lines", col = "chartreuse3", lwd=3)

# plotting late displacement over sufficient wireframe day 45
open3d(windowRect=c(20,50,1000,700), box = F, asp=F) # L, T, R, B
spheres3d(suffmu45, col = "black", radius=0.003, alpha=1, main = "")
for(i in 1:nrow(wire)){
  suffmu45.2 <- rbind(suffmu45[wire[i,1],],
                      suffmu45[wire[i,2],])
  lines3d(suffmu45.2, col=suffmu.col, lwd=2)
}
for(i in 1:nrow(suffmu45)) arrow3d(suffmu45[i,], late45mag[i,], 
                                   type = "lines", col = "chartreuse3", lwd=3)
  
# plotting late displacement over sufficient wireframe day 60
open3d(windowRect=c(20,50,1000,700), box = F, asp=F) # L, T, R, B
spheres3d(suffmu60, col = "black", radius=0.003, alpha=1, main = "")
for(i in 1:nrow(wire)){
  suffmu60.2 <- rbind(suffmu60[wire[i,1],],
                      suffmu60[wire[i,2],])
  lines3d(suffmu60.2, col=suffmu.col, lwd=2)
}
for(i in 1:nrow(suffmu60)) arrow3d(suffmu60[i,], late60mag[i,], 
                                   type = "lines", col = "chartreuse3", lwd=3)



########Day 21 composite########
# plotting Deficient displacement over sufficient wireframe day 21 
open3d(windowRect=c(20,50,1000,700), box = F, asp=F) # L, T, R, B
spheres3d(suffmu21, col = "black", radius=0.003, alpha=1, main = "")
for(i in 1:nrow(wire)){
  suffmu21.2 <- rbind(suffmu21[wire[i,1],],
                      suffmu21[wire[i,2],])
  lines3d(suffmu21.2, col=suffmu.col, lwd=2)
}
for(i in 1:nrow(suffmu21)) arrow3d(suffmu21[i,], def21mag[i,], 
                                   type = "lines", col = "blue", lwd=3)
for(i in 1:nrow(wire)){
  suffmu21.2 <- rbind(suffmu21[wire[i,1],],
                      suffmu21[wire[i,2],])
  lines3d(suffmu21.2, col=suffmu.col, lwd=2)
}
for(i in 1:nrow(suffmu21)) arrow3d(suffmu21[i,], early21mag[i,], 
                                   type = "lines", col = "red", lwd=3)
for(i in 1:nrow(wire)){
  suffmu21.2 <- rbind(suffmu21[wire[i,1],],
                      suffmu21[wire[i,2],])
  lines3d(suffmu21.2, col=suffmu.col, lwd=2)
}
for(i in 1:nrow(suffmu21)) arrow3d(suffmu21[i,], late21mag[i,], 
                                   type = "lines", col = "chartreuse3", lwd=3)


#######Day 28 composite########
# plotting Deficient displacement over sufficient wireframe day 28 
open3d(windowRect=c(20,50,1000,700), box = F, asp=F) # L, T, R, B
spheres3d(suffmu28, col = "black", radius=0.003, alpha=1, main = "")
for(i in 1:nrow(wire)){
  suffmu28.2 <- rbind(suffmu28[wire[i,1],],
                      suffmu28[wire[i,2],])
  lines3d(suffmu28.2, col=suffmu.col, lwd=2)
}
for(i in 1:nrow(suffmu28)) arrow3d(suffmu28[i,], def28mag[i,], 
                                   type = "lines", col = "blue", lwd=3)
for(i in 1:nrow(wire)){
  suffmu28.2 <- rbind(suffmu28[wire[i,1],],
                      suffmu28[wire[i,2],])
  lines3d(suffmu28.2, col=suffmu.col, lwd=2)
}
for(i in 1:nrow(suffmu28)) arrow3d(suffmu28[i,], early28mag[i,], 
                                   type = "lines", col = "red", lwd=3)
for(i in 1:nrow(wire)){
  suffmu28.2 <- rbind(suffmu28[wire[i,1],],
                      suffmu28[wire[i,2],])
  lines3d(suffmu28.2, col=suffmu.col, lwd=2)
}
for(i in 1:nrow(suffmu28)) arrow3d(suffmu28[i,], late28mag[i,], 
                                   type = "lines", col = "chartreuse3", lwd=3)

########Day 35 composite########
# plotting Deficient displacement over sufficient wireframe day 35 
open3d(windowRect=c(20,50,1000,700), box = F, asp=F) # L, T, R, B
spheres3d(suffmu35, col = "black", radius=0.003, alpha=1, main = "")
for(i in 1:nrow(wire)){
  suffmu35.2 <- rbind(suffmu35[wire[i,1],],
                      suffmu35[wire[i,2],])
  lines3d(suffmu35.2, col=suffmu.col, lwd=2)
}
for(i in 1:nrow(suffmu35)) arrow3d(suffmu35[i,], def35mag[i,], 
                                   type = "lines", col = "blue", lwd=3)
for(i in 1:nrow(wire)){
  suffmu35.2 <- rbind(suffmu35[wire[i,1],],
                      suffmu35[wire[i,2],])
  lines3d(suffmu35.2, col=suffmu.col, lwd=2)
}
for(i in 1:nrow(suffmu35)) arrow3d(suffmu35[i,], early35mag[i,], 
                                   type = "lines", col = "red", lwd=3)
for(i in 1:nrow(wire)){
  suffmu35.2 <- rbind(suffmu35[wire[i,1],],
                      suffmu35[wire[i,2],])
  lines3d(suffmu35.2, col=suffmu.col, lwd=2)
}
for(i in 1:nrow(suffmu35)) arrow3d(suffmu35[i,], late35mag[i,], 
                                   type = "lines", col = "chartreuse3", lwd=3)

########Day 45 composite#########
# plotting Deficient displacement over sufficient wireframe day 45 
open3d(windowRect=c(20,50,1000,700), box = F, asp=F) # L, T, R, B
spheres3d(suffmu45, col = "black", radius=0.003, alpha=1, main = "")
for(i in 1:nrow(wire)){
  suffmu45.2 <- rbind(suffmu45[wire[i,1],],
                      suffmu45[wire[i,2],])
  lines3d(suffmu45.2, col=suffmu.col, lwd=2)
}
for(i in 1:nrow(suffmu45)) arrow3d(suffmu45[i,], def45mag[i,], 
                                   type = "lines", col = "blue", lwd=3)
for(i in 1:nrow(wire)){
  suffmu45.2 <- rbind(suffmu45[wire[i,1],],
                      suffmu45[wire[i,2],])
  lines3d(suffmu45.2, col=suffmu.col, lwd=2)
}
for(i in 1:nrow(suffmu45)) arrow3d(suffmu45[i,], early45mag[i,], 
                                   type = "lines", col = "red", lwd=3)
for(i in 1:nrow(wire)){
  suffmu45.2 <- rbind(suffmu45[wire[i,1],],
                      suffmu45[wire[i,2],])
  lines3d(suffmu45.2, col=suffmu.col, lwd=2)
}
for(i in 1:nrow(suffmu45)) arrow3d(suffmu45[i,], late45mag[i,], 
                                   type = "lines", col = "chartreuse3", lwd=3)


##########Day 60 composite##########
# plotting Deficient displacement over sufficient wireframe day 60 
open3d(windowRect=c(20,50,1000,700), box = F, asp=F) # L, T, R, B
spheres3d(suffmu60, col = "black", radius=0.003, alpha=1, main = "")
for(i in 1:nrow(wire)){
  suffmu60.2 <- rbind(suffmu60[wire[i,1],],
                      suffmu60[wire[i,2],])
  lines3d(suffmu60.2, col=suffmu.col, lwd=2)
}
for(i in 1:nrow(suffmu60)) arrow3d(suffmu60[i,], def60mag[i,], 
                                   type = "lines", col = "blue", lwd=3)
for(i in 1:nrow(wire)){
  suffmu60.2 <- rbind(suffmu60[wire[i,1],],
                      suffmu60[wire[i,2],])
  lines3d(suffmu60.2, col=suffmu.col, lwd=2)
}
for(i in 1:nrow(suffmu60)) arrow3d(suffmu60[i,], early60mag[i,], 
                                   type = "lines", col = "red", lwd=3)
for(i in 1:nrow(wire)){
  suffmu60.2 <- rbind(suffmu60[wire[i,1],],
                      suffmu60[wire[i,2],])
  lines3d(suffmu60.2, col=suffmu.col, lwd=2)
}
for(i in 1:nrow(suffmu60)) arrow3d(suffmu60[i,], late60mag[i,], 
                                   type = "lines", col = "chartreuse3", lwd=3)







