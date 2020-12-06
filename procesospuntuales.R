#Cargar Paquetes ####

library(dplyr)
library(spatstat)
library(sp)
library(rgdal)
library(raster)
library(dplyr)
library(maptools)
library(haven)

# Cargar datos ####
data <- read_dta("megabasesecundariaActualizada29-03-19.dta")

mapa <-readOGR("Distritos_R", "Distritosv2")
mapa <- spTransform(mapa, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# Limpieza de datos ####

# Selección de unidades de interés ####
data1 <- data %>% filter(!is.na(X2))
tabla_resumen <- data1 %>% 
  dplyr::select(starts_with("desa"),starts_with("reprot"), cdpr15, cdcan15,cddis15, zona15, rama15, codigo15) %>% 
  dplyr::select(contains("15"))
summary(tabla_resumen)
rm(data)
desa =data1 %>% 
  dplyr::select(starts_with("desert"),starts_with("miit"), zona15, rama15, codigo15,nombre, nombre_ins,X2,Y2) %>% 
  dplyr::select(contains("15"), zona15, rama15, codigo15,nombre, nombre_ins,X2,Y2) %>% 
  dplyr::select(!contains("h"), zona15, rama15, codigo15,nombre, nombre_ins,X2,Y2) %>% 
  dplyr::filter(rama15 == 11)

#### Tabla de descriptivos  #####
resumen =desa %>%  
  dplyr::mutate(desa_total = desa[,1] %>% rowSums(),
                matri_total = desa[,8] %>% rowSums(),
                desa_porc =( desa_total/matri_total)*100)
sd(resumen$desa_porc, na.rm =T)
rm(resumen) #Para no consumir memoria innecesaria
#### Filtrar centros educativos de interés #####
desa =desa %>%  
  dplyr::mutate(desa_total = desa[,1] %>% rowSums(),
                matri_total = desa[,8] %>% rowSums(),
                desa_porc =( desa_total/matri_total)*100) %>% 
  filter(desa_porc > 3) # Centros con desa % en 2015 mayor al 3%
rm(tabla_resumen)

###Convertir a Spatial #### 
sp::coordinates(desa)=~X2+Y2
crs(desa) <- CRS("+proj=longlat +datum=WGS84")
plot(mapa) + plot(desa, add = T)

#### Analisis de procesos puntuales ####

zero <- zerodist(desa)
length(unique(zero[,1]))  ## detecta los duplicados jeje 

### Análisis descriptivo ####

## Cálculo de los índices ####
media_centroX <- mean(desa@coords[,1])
media_centroy <- mean(desa@coords[,2])

sd_centroX <- mean(desa@coords[,1])
sd_centroy <- mean(desa@coords[,2])
standard_distance <- sqrt(sum(((desa@coords[,1]-
                                  media_centroX)^2+(desa@coords[,2]-
                                                      media_centroy)^2))/(nrow(desa)))


plot(desa,pch="+",cex=0.5,main="")
plot(mapa,add=T)
points(media_centroX,media_centroy,col="red",pch=16)
plotrix::draw.circle(media_centroX,media_centroy,radius=standard_distance,border="red",lwd=2)
pdf("plots/Fig1.pdf")
# Eliminar duplicados de los datos ####
desa <- desa %>% remove.duplicates()
# Definir CRS
crs(desa) <- CRS("+proj=longlat +datum=WGS84")
#Definiri la ventana para las observaciones ####
mapa.utm <- spTransform(mapa,CRS("+init=epsg:32630"))
desa.utm <- spTransform(desa,CRS("+init=epsg:32630"))
window <- spatstat::as.owin(mapa.utm)
# Calcular la intensidad ####
desa.ppp <- spatstat::ppp(x=desa.utm@coords[,1],y=desa.utm@coords[,2],window=window)
desa.ppp$n/sum(sapply(slot(mapa.utm, "polygons"), slot, "area"))
desa.ppp$n/sum(sapply(slot(mapa.utm, "polygons"), slot, "area"))
#### plot 1 ####
plot(desa.ppp,pch="+",cex=0.5,main="Exclusión educativa")
plot(spatstat::quadratcount(desa.ppp, nx = 4, ny = 4),add=T,col="blue")


### Mapa intensidad ####

Local.Intensity <- data.frame(Borough=factor(),Number=numeric())
## Provincia
pdf("plots/Fig1a.pdf")
mapa.utm$nom_prov[mapa.utm$nom_prov=="LimÃ³n"] = "Limón"
mapa.utm$nom_prov[mapa.utm$nom_prov=="San JosÃ©"] = "San José"

for(i in unique(mapa.utm$nom_prov)){
  sub.pol <- mapa.utm[mapa.utm$nom_prov==i,]
  sub.ppp <- ppp(x=desa.ppp$x,y=desa.ppp$y,window=as.owin(sub.pol))
  Local.Intensity <-rbind(Local.Intensity,
                          data.frame(Borough=factor(i,levels=unique(mapa.utm$nom_prov)),
                                     Number=sub.ppp$n))
}

colorScale <-plotrix::color.scale(Local.Intensity[order(Local.Intensity[,2]),2],color.spec="rgb",extremes=c("green","red"),alpha=0.8)
par(mfrow=c(1,1))
barplot(Local.Intensity[order(Local.Intensity[,2]),2],names.arg=Local.Intensity[order(Local.Intensity[,2]),1],horiz=T,las=2,space=1,col=colorScale)
dev.off()
## Canton
Local.Intensity <- data.frame(Borough=factor(),Number=numeric())
for(i in unique(mapa.utm$nombre)){
  sub.pol <- mapa.utm[mapa.utm$nombre==i,]
  sub.ppp <- ppp(x=desa.ppp$x,y=desa.ppp$y,window=as.owin(sub.pol))
  Local.Intensity <-rbind(Local.Intensity,
                          data.frame(Borough=factor(i,levels=unique(mapa.utm$nombre)),
                                     Number=sub.ppp$n))
}
Local.Intensity %>% arrange(desc(Number)) %>% head(10)

#### Mapa de calor ####
#Este mapa dura mucho en correr, si da problema se puede hacer la versión que considera 
# solo la provincia de San Jose
par(mfrow=c(2,2))
plot(density.ppp(desa.ppp, sigma =
                   bw.diggle(desa.ppp),edge=T),main=paste("h =",round(bw.diggle(desa.ppp),2)))
plot(density.ppp(desa.ppp, sigma = bw.ppl(desa.ppp),edge=T),main=paste("h
=",round(bw.ppl(desa.ppp),2)))
plot(density.ppp(desa.ppp, sigma =
                   bw.scott(desa.ppp)[2],edge=T),main=paste("h
=",round(bw.scott(desa.ppp)[2],2)))
plot(density.ppp(desa.ppp, sigma =
                   bw.scott(desa.ppp)[1],edge=T),main=paste("h
=",round(bw.scott(desa.ppp)[1],2)))


plot(Gest(desa.ppp),main="Exclusión % Colegial alta")


#### Visualización de una provincia San José ####
projection(desa)=projection(mapa)
overlay <- over(desa,mapa)
desa@data <- cbind(desa@data, overlay)
names(desa@data)[27] <- "nombre_cant"
mapa_sub = subset(mapa, mapa@data$nom_prov == "San JosÃ©")
sj = subset(desa, desa@data$nom_prov == "San JosÃ©")
mapa.utm <- spTransform(mapa_sub,CRS("+init=epsg:32630"))
desa.utm <- spTransform(sj,CRS("+init=epsg:32630"))
window <- spatstat::as.owin(mapa.utm)
# Calcular la intensidad ####
desa.ppp <- spatstat::ppp(x=desa.utm@coords[,1],y=desa.utm@coords[,2],window=window)
desa.ppp$n/sum(sapply(slot(mapa.utm, "polygons"), slot, "area"))
desa.ppp$n/sum(sapply(slot(mapa.utm, "polygons"), slot, "area"))
#### plot 1 ####
plot(desa.ppp,pch="+",cex=0.5,main="Exclusión educativa")
plot(spatstat::quadratcount(desa.ppp, nx = 4, ny = 4),add=T,col="blue")
#### Mapas de calor ####
pdf("plots/Fig2.pdf")
par(mfrow=c(2,2))
plot(density.ppp(desa.ppp, sigma =
                   bw.diggle(desa.ppp),edge=T),main=paste("h =",round(bw.diggle(desa.ppp),2)))
plot(density.ppp(desa.ppp, sigma = bw.ppl(desa.ppp),edge=T),main=paste("h
=",round(bw.ppl(desa.ppp),2)))
plot(density.ppp(desa.ppp, sigma =
                   bw.scott(desa.ppp)[2],edge=T),main=paste("h
=",round(bw.scott(desa.ppp)[2],2)))
plot(density.ppp(desa.ppp, sigma =
                   bw.scott(desa.ppp)[1],edge=T),main=paste("h
=",round(bw.scott(desa.ppp)[1],2)))
dev.off()

