#set working directory
setwd("~/Documents/University/MScMarineSciences/CSIRO Internship/Ari_s project")

#load data
load("~/Documents/University/MScMarineSciences/CSIRO Internship/Ari_s project/final_environment_with_predictions.RData")

library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
library(texreg)
library(ggspatial)

#Figure 1 (not working with updated version of RStudio)
setEPS()
postscript("figure1.eps")
world <- ne_countries(scale = "medium", returnclass = "sf")
theme_set(theme_bw())
ggplot(data = world) +
  geom_sf(fill = "lemonchiffon1") +
  labs(x ="Longitude", cex = 1.2) + labs(y ="Latitude", cex = 1.2) +
  coord_sf(xlim = c(103, 163), ylim = c(-9, -45), expand = FALSE) +
  coord_sf(xlim = c(103, 163), ylim = c(-9, -45), expand = FALSE) +
  geom_point(data = Data, aes(x = BCH_LONGITUDE_DD, y = BCH_LATITUDE_DD), color = "darkgreen", cex = 1.2)
  annotation_north_arrow(location = "tr", which_north = "true", style = north_arrow_fancy_orienteering)
dev.off()

#Figure 3A
setEPS()
postscript("figure3A.eps")
plot(DataLongPositives$SectionNumber, ylim = c(0, 250), col = "snow3", xlab = "Section Number", ylab = "Frequency", cex.lab =1.6, cex.axis =1.2)
title(main = list("Presence of debris", cex = 1.6))
text(0.5,240,"A", cex = 1.6)
dev.off()

#Figure 3B 
setEPS()
postscript("figure3B.eps")
DataLongPositives$SizeCategory <- factor(DataLongPositives$SizeCategory, ordered = T)
plot(DataLongPositives$SizeCategory, ylim = c(0,350), col = c("black" , "dodgerblue4", "forestgreen", "deepskyblue",  "gold", "darkorange1"), xlab = "Size Category [cm]", ylab = "Frequency", cex.lab = 1.6, xaxt = "n")
title(main = "First item", cex.main = 1.6)
axis(1, at=seq(0.7,6.7,by=1.2), labels = c("0.2-1.2","1.2-2.4","2.4-4.9","4.9-9.6","9.6-19.2","19.2-36.3"), tick = F, cex.axis = 1.1)
text(0.3, 335, "B", cex = 1.6)
dev.off()

#Figure 3C
posi <- as.data.frame(cbind(DataLongPositives$NumericSizeCategory, DataLongPositives$NumericSectionNumber))
names(posi) <- c("SizeCategory", "SectionNumber")
posi$SizeCategory[posi$SectionNumber == 2] <- NA
posi$SizeCategory[posi$SectionNumber == 3] <- NA
posi$SizeCategory[posi$SectionNumber == 4] <- NA
posi$SizeCategory[posi$SectionNumber == 6] <- NA
posi$SizeCategory[posi$SectionNumber == 7] <- NA
posi$SizeCategory[posi$SectionNumber == 8] <- NA
posi$SizeCategory[posi$SectionNumber == 9] <- NA
posi$SectionNumber[posi$SectionNumber == 2] <- NA
posi$SectionNumber[posi$SectionNumber == 3] <- NA
posi$SectionNumber[posi$SectionNumber == 4] <- NA
posi$SectionNumber[posi$SectionNumber == 6] <- NA
posi$SectionNumber[posi$SectionNumber == 7] <- NA
posi$SectionNumber[posi$SectionNumber == 8] <- NA
posi$SectionNumber[posi$SectionNumber == 9] <- NA
posi <- posi[!is.na(posi$SizeCategory),]
posi_tab <- table(posi$SizeCategory,posi$SectionNumber)
ks_1 <- ks.test(posi_tab[,1],posi_tab[,2])
ks_2 <- ks.test(posi_tab[,1], posi_tab[,3])
ks_2
setEPS()
postscript("figure3C.eps")
barplot(posi_tab, col = c("black" , "dodgerblue4", "forestgreen", "deepskyblue",  "gold", "darkorange1"), beside = T, ylim = c(0, 60), xlab = "Section Number", ylab = "Frequency", cex.lab = 1.6, cex.axis = 1.6)
text(1, 57, "C", cex = 1.6)
text(3, 50, "a", cex = 1.2)
text(9.5, 50, "a", cex = 1.2)
text(16.5, 50, "b", cex = 1.2)
title(main = "Size category of the first item", cex.main = 1.6)
dev.off()


#Figure 4 (option 3) --> this one is in the paper
setEPS()
postscript("figure4.eps")
plot(NewData$NumericSectionNumber[NewData$LowerBSize == 0], (NewData$DistMid[NewData$LowerBSize == 0]-NewData$SERelDist[NewData$LowerBSize == 0]), xlab = "Section Number", ylab = "Relative Distance", ylim = c(-0.15,0.15), cex.lab =1.6, cex.axis = 1.6, xaxt = "n", type = "l", lty = 2)
lines(NewData$NumericSectionNumber[NewData$LowerBSize == 0], NewData$DistMid[NewData$LowerBSize == 0]+NewData$SERelDist[NewData$LowerBSize == 0], type = "l", lty = 2)
grid(nx = NA, ny = NULL)
axis(1, at = 1:10, cex.axis =1.6)
#title(main = "Departure from midpoint", cex =1.6)
lines(NewData$NumericSectionNumber[NewData$LowerBSize == 1.2], NewData$DistMid[NewData$LowerBSize == 1.2]-NewData$SERelDist[NewData$LowerBSize == 1.2], type = "l", lty = 2, col = "dodgerblue4" )
lines(NewData$NumericSectionNumber[NewData$LowerBSize == 1.2], NewData$DistMid[NewData$LowerBSize == 1.2]+NewData$SERelDist[NewData$LowerBSize == 1.2], type = "l", lty = 2, col = "dodgerblue4" )
lines(NewData$NumericSectionNumber[NewData$LowerBSize == 2.4], NewData$DistMid[NewData$LowerBSize == 2.4]-NewData$SERelDist[NewData$LowerBSize == 2.4], type = "l", lty = 2, col = "forestgreen" )
lines(NewData$NumericSectionNumber[NewData$LowerBSize == 2.4], NewData$DistMid[NewData$LowerBSize == 2.4]+NewData$SERelDist[NewData$LowerBSize == 2.4], type = "l", lty = 2, col = "forestgreen" )
lines(NewData$NumericSectionNumber[NewData$LowerBSize == 4.9], NewData$DistMid[NewData$LowerBSize == 4.9]-NewData$SERelDist[NewData$LowerBSize == 4.9], type = "l", lty = 2, col = "deepskyblue" )
lines(NewData$NumericSectionNumber[NewData$LowerBSize == 4.9], NewData$DistMid[NewData$LowerBSize == 4.9]+NewData$SERelDist[NewData$LowerBSize == 4.9], type = "l", lty = 2, col = "deepskyblue" )
lines(NewData$NumericSectionNumber[NewData$LowerBSize == 9.6], NewData$DistMid[NewData$LowerBSize == 9.6]-NewData$SERelDist[NewData$LowerBSize == 9.6], type = "l", lty = 2, col = "gold" )
lines(NewData$NumericSectionNumber[NewData$LowerBSize == 9.6], NewData$DistMid[NewData$LowerBSize == 9.6]+NewData$SERelDist[NewData$LowerBSize == 9.6], type = "l", lty = 2, col = "gold" )
lines(NewData$NumericSectionNumber[NewData$LowerBSize == 19.2], NewData$DistMid[NewData$LowerBSize == 19.2]-NewData$SERelDist[NewData$LowerBSize == 19.2], type = "l", lty = 2, col = "darkorange1" )
lines(NewData$NumericSectionNumber[NewData$LowerBSize == 19.2], NewData$DistMid[NewData$LowerBSize == 19.2]+NewData$SERelDist[NewData$LowerBSize == 19.2], type = "l", lty = 2, col = "darkorange1" )
lines(NewData$NumericSectionNumber[NewData$LowerBSize == 0], NewData$DistMid[NewData$LowerBSize == 0],type = "l", lwd=3, col = "black")
lines(NewData$NumericSectionNumber[NewData$LowerBSize == 1.2], NewData$DistMid[NewData$LowerBSize == 1.2], type = "l", lwd = 3, pch = 16, cex = 1.6, col = "dodgerblue4")
lines(NewData$NumericSectionNumber[NewData$LowerBSize == 2.4], NewData$DistMid[NewData$LowerBSize == 2.4], type = "l", lwd = 3, pch = 16, cex = 1.6, col = "forestgreen")
lines(NewData$NumericSectionNumber[NewData$LowerBSize == 4.9], NewData$DistMid[NewData$LowerBSize == 4.9], type = "l", lwd = 3, pch = 16, cex = 1.6, col = "deepskyblue")
lines(NewData$NumericSectionNumber[NewData$LowerBSize == 9.6], NewData$DistMid[NewData$LowerBSize == 9.6], type = "l", lwd = 3, pch = 16, cex = 1.6, col = "gold")
lines(NewData$NumericSectionNumber[NewData$LowerBSize == 19.2], NewData$DistMid[NewData$LowerBSize == 19.2], type = "l", lwd = 3, pch = 16, cex = 1.6, col = "darkorange1")
legend(x = 1, y = 0.15, 
       legend = c("0.1 - 1.2 cm", "1.2 - 2.4 cm", "2.4 - 4.9 cm", "4.9 - 9.6 cm", "9.6 - 19.2 cm", "19.2 - 36.3 cm"),
       col = c("black" , "dodgerblue4", "forestgreen", "deepskyblue",  "gold", "darkorange1"),
       pch = 16,
       pt.cex = 1.2,
       cex = 1.2,
       y.intersp = 0.8)
dev.off()
