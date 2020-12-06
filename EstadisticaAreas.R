#### Estadística de Áreas ####
#### 
#### library(ggpubr)
library(dplyr)
library(ggplot2)
library(spdep)
library(haven)
library(rgdal)
library(RColorBrewer)
#Cargar datos 
mapa <-readOGR("Distritos_R", "Distritosv2")
mapa <- spTransform(mapa, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
#Eliminar isla Chira porque es el único poligono sin vecino al ser una isla
mapa <- subset(mapa,mapa$nom_distr!="CHIRA")
data <- read_dta("megabasesecundariaActualizada29-03-19.dta")
# Preparar datos 
# Selección de unidades de interés ####
data1 <- data %>% filter(!is.na(X2))
tabla_resumen <- data1 %>% 
  dplyr::select(starts_with("desa"),starts_with("reprot"), cdpr15, cdcan15,cddis15, zona15, rama15, codigo15) %>% 
  dplyr::select(contains("15"))
summary(tabla_resumen)
rm(data, tabla_resumen)
desa =data1 %>% 
  dplyr::select(starts_with("desert"),starts_with("miit"), zona15, rama15, codigo15,nombre, nombre_ins,X2,Y2) %>% 
  dplyr::select(contains("15"), zona15, rama15, codigo15,nombre, nombre_ins,X2,Y2) %>% 
  dplyr::select(!contains("h"), zona15, rama15, codigo15,nombre, nombre_ins,X2,Y2) %>% 
  #dplyr::select(-desa_15) %>% 
  dplyr::filter(rama15 == 11)
desa =desa %>%  
  dplyr::mutate(desa_total = desa[,1] %>% rowSums(),
                matri_total = desa[,8] %>% rowSums(),
                desa_porc =( desa_total/matri_total)*100) %>% 
  filter(desa_porc > 3) # Centros con desa % en 2015 mayor al 3%
names(desa)

###Convertir a Spatial #### 
sp::coordinates(desa)=~X2+Y2
crs(desa) <- CRS("+proj=longlat +datum=WGS84")

#Agregar datos de exclusión escolar al mapa
projection(desa)=projection(mapa)
overlay <- over(desa,mapa)
desa@data <- cbind(desa@data, overlay)
names(desa@data)[27] <- "nombre_cant"
#Crear variables resumen por poligono
tabla <- desa@data %>% 
  dplyr::group_by(cod_dta) %>% 
  dplyr::summarise(prom_desa = mean(desa_porc,na.rm=T ),
                   prom_matri = mean(matri_total, na.rm=T),
                   prom_total =mean(desa_total, na.rm=T))

completar <- left_join(mapa@data, tabla, by = "cod_dta")
completar[is.na(completar)] <- 0
mapa@data <- completar
xy <- coordinates(mapa)

####  Generar el modelo inicial ####

# Por distancia:
plot(w, xy, col='red', lwd=2)
w<- poly2nb(mapa)
plot(Sy1_nb, xy, col='red', lwd=2)



## Calcular los pesos de los vecinos 
Sy0_lw_W <- nb2listw(w)
Sy0_lw_W
summary(sapply(Sy0_lw_W$weights, sum))


#### Pruebas de autocorrelación ####
set.seed(2905)
n <- length(w) -1
rho <- 0.5
autocorr_x <- invIrW(Sy0_lw_W, rho) %*% mapa$prom_desa


### Grafico de autocorrelación ####
oopar <- par(mfrow=c(1,2), mar=c(4,4,3,2)+0.1)
plot(autocorr_x, stats::lag(Sy0_lw_W, autocorr_x),
     xlab="", ylab="",
     main="Autocorrelated random variable", cex.main=0.8, cex.lab=0.8)
lines(lowess(autocorr_x, stats::lag(Sy0_lw_W, autocorr_x)), lty=2, lwd=2)
# Test de moran
moran.test(mapa$prom_desa, listw=Sy0_lw_W)
-1 / (n-1)

#### Estimar índice de moran ####
set.seed(1234)
bperm <- moran.mc(mapa$prom_desa, listw=Sy0_lw_W, nsim=999)
bperm
#Si tengo evidencia para rechazar la hipótesis nula, entonces significa que la
# locación original importa y que por lo tanto si existe autocorrelación espacial.

#### Mapas para clusters #####
oopar <- par(mfrow=c(1,2))
msp <- moran.plot(mapa$prom_desa, listw=nb2listw(w, style="S"), quiet=TRUE)
title("Moran scatterplot")

pdf("plots/Fig3.pdf")
infl <- apply(msp["is_inf"], 1, any)
x <- abs(mapa$prom_desa)
lhx <- cut(x, breaks=c(min(x), mean(x), max(x)), labels=c("L", "H"), include.lowest=TRUE)
wx <- stats::lag(nb2listw(w, style="S"), mapa$prom_desa)
lhwx <- cut(wx, breaks=c(min(wx), mean(wx), max(wx)), labels=c("L", "H"), include.lowest=TRUE)
lhlh <- interaction(lhx, lhwx, infl, drop=TRUE)
cols <- rep(1, length(lhlh))
cols[lhlh == "H.L.TRUE"] <- 2
cols[lhlh == "L.H.TRUE"] <- 3
cols[lhlh == "H.H.TRUE"] <- 4
plot(mapa, col=brewer.pal(4, "Accent")[cols])
legend("topright", legend=c("None", "HL", "LH", "HH"), fill=brewer.pal(4, "Accent"), bty="n", cex=0.8, y.intersp=0.8)
title("Distritos de influencia")
dev.off()

##### Detectar si hay clusters #####
lm1 <- localmoran(mapa$prom_desa, listw=nb2listw(w, style="S"))
r <- sum(mapa$prom_desa)/sum(mapa$prom_matri)
rni <- r*mapa$prom_matri
lw <- nb2listw(w)
sdCR <- (mapa$prom_desa - rni)/sqrt(rni)
wsdCR <- stats::lag(nb2listw(w, style="C"), mapa$prom_desa)
I_CR <- sdCR * wsdCR


#### Mapa de clusters ####
gry <- c(rev(brewer.pal(8, "Reds")[1:6]), brewer.pal(6, "Blues"))
mapa$Standard <- lm1[,1]
mapa$"Constant_risk" <- I_CR
spplot(mapa, c("Standard", "Constant_risk"), at=c(-2.5,-1.4,-0.6,-0.2,0,0.2,0.6,4,7), col.regions=colorRampPalette(gry)(8))
spplot(mapa, c( "Constant_risk"), at=c(-2.5,-1.4,-0.6,-0.2,0,0.2,0.6,4,7), col.regions=colorRampPalette(gry)(8))
