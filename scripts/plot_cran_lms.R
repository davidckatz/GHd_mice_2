####This script was written to help me figure out which cranial landmarks belong to which module. 
###It plots the 51 cranial landmarks onto a mesh of an adult mouse skull. that's it that's all.

################
#Packages
library(geomorph)
library(Morpho)
library(rgl)

################

#Load in data
setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
mesh <- file2mesh(filename = file.path("data", "CranA4760Des.ply"))
df <- read.csv(file = file.path("data", "A4760LM_fixed.csv"))
lmnames <- read.csv(file = file.path("data", "CranialLMGuide.csv"))

#Grab landmark data 
lms <- df[,2:4]

#Code from Day 3 of David's course:
# First use this code to plot mesh, landmarks, and landmark #s
open3d(windowRect=c(20,50,600,400)) # L, T, R, B
plot3d(mesh, col = "lightgrey", alpha=0.3,
       axes = F, box = F, asp=F,
       xlab = "", ylab = "", zlab = "", main = "")
# plot those points in red on the skull
spheres3d(lms, radius = .2, color = "red")
text3d(x = lms, texts =lmnames[,2], adj=c(0.75,1.5), cex=1.2, font=2)

# lm modules
face.lms <- which(lmnames$Module=="Face"|lmnames$Module2=="Face")
base.lms <- which(lmnames$Module=="Base"|lmnames$Module2=="Base")
vault.lms <- which(lmnames$Module=="Vault"|lmnames$Module2=="Vault")

# plot Face lms
open3d(windowRect=c(20,50,600,400)) # L, T, R, B
plot3d(mesh, col = "lightgrey", alpha=0.3,
       axes = F, box = F, asp=F,
       xlab = "", ylab = "", zlab = "", main = "FACE")
# plot those points in red on the skull
spheres3d(lms[face.lms,], radius = .2, color = "red")
text3d(x = lms[face.lms,], texts =lmnames[face.lms,2], adj=c(0.75,1.5), cex=1.2, font=2)

# plot Vault lms
open3d(windowRect=c(20,50,600,400)) # L, T, R, B
plot3d(mesh, col = "lightgrey", alpha=0.3,
       axes = F, box = F, asp=F,
       xlab = "", ylab = "", zlab = "", main = "VAULT")
# plot those points in red on the skull
spheres3d(lms[vault.lms,], radius = .2, color = "green")
text3d(x = lms[vault.lms,], texts =lmnames[vault.lms,2], adj=c(0.75,1.5), cex=1.2, font=2)

# plot Base lms
open3d(windowRect=c(20,50,600,400)) # L, T, R, B
plot3d(mesh, col = "lightgrey", alpha=0.3,
       axes = F, box = F, asp=F,
       xlab = "", ylab = "", zlab = "", main = "BASE")
# plot those points in red on the skull
spheres3d(lms[base.lms,], radius = .2, color = "blue")
text3d(x = lms[base.lms,], texts =lmnames[base.lms,2], adj=c(0.75,1.5), cex=1.2, font=2)

