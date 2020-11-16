###################################################################################
#This is a script to create and plot wireframes for the GH mice, including:
#Whole skull, vault, base, and face
#Workng directory will need to be changed 

###################################################################################
#Packages
library(geomorph)
library(Morpho)
library(rgl)

####################

setwd("/Users/amandaneves/Desktop/GH")

#Read in wire frames
vault_wire <- read.csv("vault_links.csv", header = FALSE)
base_wire <- read.csv("base_links.csv", header = FALSE)
face_wire <- read.csv("face_links.csv", header = FALSE)
cran_wire <- read.csv("cran_links.csv", header = FALSE)

#Load in some files (these correspond with the files used in plot_cran_lms.R)
mesh <- file2mesh(filename = file.path("CranA4760Des.ply"))
lmnames <- read.csv(file = file.path("CranialLMGuide.csv"))
df <- read.csv(file = file.path("A4760LM_fixed.csv"))
lms <- df[,2:4]
lm.count <- 1:nrow(lms)

#Index the lm module names
# lm modules
face.lms <- which(lmnames$Module=="Face"|lmnames$Module2=="Face")
base.lms <- which(lmnames$Module=="Base"|lmnames$Module2=="Base")
vault.lms <- which(lmnames$Module=="Vault"|lmnames$Module2=="Vault")

####################
# The following code plots each wireframe onto the skull: whole skull, vault, base, face
# Links were made using the defined modules in plot_cran_lms.R as a reference 

#### WHOLE SKULL ####
open3d(windowRect=c(20,50,600,400)) # L, T, R, B
plot3d(mesh, col = "lightgrey", alpha=0.3,
       axes = F, box = F, asp=F,
       xlab = "", ylab = "", zlab = "", main = "WHOLE SKULL")
# plot those points in red on the skull
spheres3d(lms, radius = .2, color = "red")
text3d(x = lms, texts =1:nrow(lms), adj=c(0.75,1.5), cex=1.2, font=2)
for(i in 1:nrow(cran_wire)){
  connect.2 <- rbind(lms[cran_wire[i,1],],
                     lms[cran_wire[i,2],])
  lines3d(connect.2, col="red", lwd=2)
}

#### VAULT ####
#Plot vault lms to double-check things
open3d(windowRect=c(20,50,600,400)) # L, T, R, B
plot3d(mesh, col = "lightgrey", alpha=0.3, 
       axes = F, box = F, asp=F,
       xlab = "", ylab = "", zlab = "", main = "VAULT")
# plot vault points in green on the skull
spheres3d(lms[vault.lms,], radius = .2, color = "green")
#other landmarks 
spheres3d(lms[-vault.lms,], radius = .2, color = "black")
# face lm numbers
text3d(x = lms[vault.lms,], texts =lm.count[vault.lms], adj=c(0.75,1.5), cex=1.2, font=2)
# plot wireframe
for(i in 1:nrow(vault_wire)){
  connect.2 <- rbind(lms[vault_wire[i,1],],
                     lms[vault_wire[i,2],])
  lines3d(connect.2, col="green", lwd=2)
}

#### BASE ####
open3d(windowRect=c(20,50,600,400)) # L, T, R, B
plot3d(mesh, col = "lightgrey", alpha=0.3,
       axes = F, box = F, asp=F,
       xlab = "", ylab = "", zlab = "", main = "BASE")
# plot base points in blue on the skull
spheres3d(lms[base.lms,], radius = .2, color = "blue")
# other landmarks
spheres3d(lms[-base.lms,], radius = .2, color = "black")
#Just base lm numbers
text3d(x = lms[base.lms,], texts =lm.count[base.lms], adj=c(0.75,1.5), cex=1.2, font=2)
# plot wireframe
for(i in 1:nrow(base_wire)){
  connect.2 <- rbind(lms[base_wire[i,1],],
                     lms[base_wire[i,2],])
  lines3d(connect.2, col="blue", lwd=2)
}

#### FACE ####
open3d(windowRect=c(20,50,600,400)) # L, T, R, B
plot3d(mesh, col = "lightgrey", alpha=0.3,
       axes = F, box = F, asp=F,
       xlab = "", ylab = "", zlab = "", main = "FACE")
# plot face points in purple on the skull
spheres3d(lms[face.lms,], radius = .2, color = "purple")
#other landmarks 
spheres3d(lms[-face.lms,], radius = .2, color = "black")
# face lm numbers
text3d(x = lms[face.lms,], texts =lm.count[face.lms], adj=c(0.75,1.5), cex=1.2, font=2)
#plot the wireframe
for(i in 1:nrow(face_wire)){
  connect.2 <- rbind(lms[face_wire[i,1],],
                     lms[face_wire[i,2],])
  lines3d(connect.2, col="purple", lwd=2)
}

#### All the modules together ####
open3d(windowRect=c(20,50,600,400)) # L, T, R, B
plot3d(mesh, col = "lightgrey", alpha=0.3,
       axes = F, box = F, asp=F,
       xlab = "", ylab = "", zlab = "", main = "Face (purple), base (blue), and vault (green")
spheres3d(lms[face.lms,], radius = .2, color = "purple")
for(i in 1:nrow(face_wire)){
  connect.2 <- rbind(lms[face_wire[i,1],],
                     lms[face_wire[i,2],])
  lines3d(connect.2, col="purple", lwd=2)
}

spheres3d(lms[base.lms,], radius = .2, color = "blue")
for(i in 1:nrow(base_wire)){
  connect.2 <- rbind(lms[base_wire[i,1],],
                     lms[base_wire[i,2],])
  lines3d(connect.2, col="blue", lwd=2)
}

spheres3d(lms[vault.lms,], radius = .2, color = "green")
for(i in 1:nrow(vault_wire)){
  connect.2 <- rbind(lms[vault_wire[i,1],],
                     lms[vault_wire[i,2],])
  lines3d(connect.2, col="green", lwd=2)
}




