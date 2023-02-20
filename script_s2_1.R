
library(tidyverse)

dexter <- read_csv("/home/marcio/Documentos/Doutorado/Dexter/DataForMarina.csv") %>%
  rename(plot_id = Sort) %>% 
  filter(cluster_membership_corrected == "Amazon Forest" | cluster_membership_corrected == "Cerrado" | cluster_membership_corrected == "SDTF")  %>%
  dplyr::rename(lon = Long10, lat = Lat10, Vegetation_type = cluster_membership_corrected, Leaf_flush = Leaf.flush, MAP = PrecAnn)



s2 <- read_csv("/home/marcio/Documentos/Doutorado/Dexter/sentinel-2/s2_do_Dexter_all.csv") %>%
  rename(evi2 = value) %>% 
  rename(time = date) %>% 
  mutate(time = as.Date(time, format = "%Y-%m-%d")) %>% 
  group_by(plot_id) %>%
  nest

# Para juntar os data frames, devemos observar se usaremos o left_join ou o right_join, pois existe diferença na ordem dos data frames
data <- left_join(dexter, s2, by="plot_id")
#right_join(s2, dexter, by="plot_id") %>% dplyr::select(Leaf_flush) %>% as.data.frame %>% tail

# Apenas um valor por mês. Se há mais de um valor, calcular a média. 
ndvi.mon <- vector("list", length = length(data$plot_id))
for (i in seq_along(data$plot_id)) {
  
  dados <- data %>%
    filter(plot_id==i) %>%
    unnest(cols = c(data)) %>%
    dplyr::select(evi2|time)
  
  time <- dados$time  
  evi2 <- dados$evi2
  
  rm(dados)
  
  df <- data.frame(evi2) %>% xts::xts(order.by = time, "Months")
  
  #df %>% TSstudio::ts_plot()
  
  ep <- xts::endpoints(xts::xts(df), on = "months")
  aux <- xts::period.apply(as.matrix(df$evi2), INDEX = ep,
                           FUN = mean, na.action=na.pass)
  
  ndvi.mon[[i]] <- aux
  
  print(paste0("Passo: ",i))
}

rname <- function(lndvi){
  rownames(lndvi) <- substr(rownames(lndvi), 1, 10)
  rownames(lndvi) <- paste(substr(rownames(lndvi), 1, 8),"01", sep = "")
  rname <- lndvi
}

ndvi.mon.n  <- lapply(ndvi.mon, FUN = rname)


ndvi.mon.n <- ndvi.mon.n %>%
  map(as.data.frame) %>% 
  map(rownames_to_column) %>% 
  map(dplyr::rename, date=rowname)


# ==========
acoords <- SpatialPoints(cbind(dexter$lon, dexter$lat), proj4string = crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

nc.chirps <- ncdf4::nc_open("/home/marina/spa paper/chirps-v2.0.monthly_1981_2021_TSA.nc", readunlim = FALSE)

# nc.chirps <- nc_open("/home/marcio/datasets/chirps-v2.0.monthly_TSA.nc", readunlim = FALSE) # opening chirps netcdf dataset for reading
?nc_open
prec.chirps<- ncdf4::ncvar_get(nc.chirps, "precip") # getting prec values
n <- length(prec.chirps)
pbar <- txtProgressBar()

for (i in 1:dim(prec.chirps)[3]){
  aux <- apply(prec.chirps[ , ,i], 1, rev)
  ChirpsAux <- raster(aux, xmn = -90, xmx = -30, ymn = -20, ymx = 15)
  if (i == 1){
    ChirpsR <- stack(ChirpsAux, nl = 1)
  }else{
    ChirpsR <- addLayer(ChirpsR, ChirpsAux)
  }
  rm(aux, ChirpsAux)
  setTxtProgressBar(pbar, i/n)
}

crs(ChirpsR) <- crs(acoords)
precip_dexter   <- raster::extract(ChirpsR, acoords)

colnames(precip_dexter) <- seq.Date(from = as.Date.factor("1981-01-01"), by="month", length.out = 495) %>% as.character


# As series temporais de cada ponto sao itens de uma lista

precipitation <- vector("list", length = dim(precip_dexter)[1])
for (i in 1:dim(precip_dexter)[1]){
  prec_aux <- data.frame(precipitation = precip_dexter[i,])
  precipitation[i] <- prec_aux %>%
    map(data.frame) %>%
    map(add_column,date = seq.Date(from = as.Date.factor("1981-01-01"), by="month", length.out = 495) %>%
          as.character)
}

precipitation <- precipitation %>% map(dplyr::rename, "rainfall"=.x..i..)

# Join ndvi and precipitation ====
# Juntar os valores de ndvi para cada ponto com as precipitaçoes para os mesmos meses.
library(padr)

evi2c <- ndvi.mon.n %>% map(mutate, date = as.Date(date, fomat="%Y-%m-%d")) %>% map(pad)

prec <- precipitation %>% map(mutate, date = as.Date(date, fomat="%Y-%m-%d"))

tabela <- purrr::map2(prec, evi2c, right_join, 
                                by = "date")


########

#tabela <- joined_prec_ndvi_1

#==================
# Contar os NAs
#==================
# Contar os NAs
n <- length(tabela)
pbar <- txtProgressBar()
soma_de_NAs <- vector("list", length = length(tabela))
for (i in seq_along(tabela)) {
  soma_de_NAs[[i]] <- tabela %>%
    pluck(i) %>%
    dplyr::select(evi2) %>% 
    is.na %>% 
    colSums
  setTxtProgressBar(pbar, i/n)
}

soma_NAs <- soma_de_NAs %>% 
  unlist

data2 <- data %>%
  add_column(soma_NAs) %>% 
  add_column(tabela) %>% 
  filter(soma_NAs<=19)

data2 %>% dim
tabela1 <- data2$tabela

# Usando na.approx

library(zoo)

evi2_aux <- vector("list", length = length(tabela1))
for (i in seq_along(tabela1)) {
  evi2_aux[[i]] <- na.approx(zoo(tabela1 %>%
                                   pluck(i) %>%
                                   dplyr::select(evi2) %>%
                                   as.matrix), tabela1 %>%
                               pluck(i) %>%
                               dplyr::select(date) %>%
                               as.matrix %>%
                               as.Date, rule=2)
}

# Adiciona coluna em cada tabela da lista
dados1 <- vector("list", length = length(tabela1))
for (i in seq_along(tabela1)) {
  dados1[[i]] <- tabela1[[i]] %>% 
    add_column(without_missing = evi2_aux[[i]] %>% as.data.frame)
}

# Acoplamento entre evi2 e chuva
#===============
library(quantmod)


# Usando o NDVI (ndvi_aux) gerado com na.approx

for (i in 1:length(dados1)) {
  evi2_aux <- pluck(dados1,i)[4] %>% as.matrix
  rainfall_aux <- pluck(dados1,i)[1] %>% as.matrix
  try(aux <- cor.test(purrr::map(rainfall_aux, Lag,0) %>% 
                        unlist %>%
                        as.numeric,
                      evi2_aux %>% unlist %>% as.numeric,
                      method = "kendall",
                      na.action = na.pass,
                      exact = F)) # função 'try' ignora os erros
  
  if ( i == 1){
    lag.0 <- c(aux$estimate)
  }else{
    lag.0 <- c(lag.0, aux$estimate)
  }
}


lag.0 %>% head

for (i in 1:length(dados1)) {
  
  try(aux <- cor.test(purrr::map(pluck(dados1,i)[1], Lag,1) %>% 
                        unlist %>%
                        as.numeric,
                      pluck(dados1,i)[4] %>% unlist %>% as.numeric,
                      method = "kendall",
                      na.action = na.pass,
                      exact = F)) # função 'try' ignora os erros
  
  if ( i == 1){
    lag.1 <- c(aux$estimate)
  }else{
    lag.1 <- c(lag.1, aux$estimate)
  }
}

lag.1 %>% head


for (i in 1:length(dados1)) {
  
  try(aux <- cor.test(purrr::map(pluck(dados1,i)[1], Lag,2) %>% 
                        unlist %>%
                        as.numeric,
                      pluck(dados1,i)[4] %>% unlist %>% as.numeric,
                      method = "kendall",
                      na.action = na.pass,
                      exact = F)) # função 'try' ignora os erros
  
  if ( i == 1){
    lag.2 <- c(aux$estimate)
  }else{
    lag.2 <- c(lag.2, aux$estimate)
  }
}

lag.2 %>% head

for (i in 1:length(dados1)) {
  
  try(aux <- cor.test(purrr::map(pluck(dados1,i)[1], Lag,3) %>% 
                        unlist %>%
                        as.numeric,
                      pluck(dados1,i)[4] %>% unlist %>% as.numeric,
                      method = "kendall",
                      na.action = na.pass,
                      exact = F)) # função 'try' ignora os erros
  
  if ( i == 1){
    lag.3 <- c(aux$estimate)
  }else{
    lag.3 <- c(lag.3, aux$estimate)
  }
}

lag.3 %>% head

for (i in 1:length(dados1)) {
  
  try(aux <- cor.test(purrr::map(pluck(dados1,i)[1], Lag,4) %>% 
                        unlist %>%
                        as.numeric,
                      pluck(dados1,i)[4] %>% unlist %>% as.numeric,
                      method = "kendall",
                      na.action = na.pass,
                      exact = F)) # função 'try' ignora os erros
  
  if ( i == 1){
    lag.4 <- c(aux$estimate)
  }else{
    lag.4 <- c(lag.4, aux$estimate)
  }
}

lag.4 %>% head

for (i in 1:length(dados1)) {
  
  try(aux <- cor.test(purrr::map(pluck(dados1,i)[1], Lag,5) %>% 
                        unlist %>%
                        as.numeric,
                      pluck(dados1,i)[4] %>% unlist %>% as.numeric,
                      method = "kendall",
                      na.action = na.pass,
                      exact = F)) # função 'try' ignora os erros
  
  if ( i == 1){
    lag.5 <- c(aux$estimate)
  }else{
    lag.5 <- c(lag.5, aux$estimate)
  }
}

lag.5 %>% head

for (i in 1:length(dados1)) {
  
  try(aux <- cor.test(purrr::map(pluck(dados1,i)[1], Lag,6) %>% 
                        unlist %>%
                        as.numeric,
                      pluck(dados1,i)[4] %>% unlist %>% as.numeric,
                      method = "kendall",
                      na.action = na.pass,
                      exact = F)) # função 'try' ignora os erros
  
  if ( i == 1){
    lag.6 <- c(aux$estimate)
  }else{
    lag.6 <- c(lag.6, aux$estimate)
  }
}

lag.6 %>% head


lags <- cbind(as.data.frame(lag.0),as.data.frame(lag.1),as.data.frame(lag.2),as.data.frame(lag.3),as.data.frame(lag.4),as.data.frame(lag.5),as.data.frame(lag.6))


# 
# mean_ndvi <- vector("list", length = length(tabela1))
# for (i in 1:length(tabela1)) {
#   aux <- purrr::map(pluck(tabela1,i)[3] %>% remove_missing, mean)
#   mean_ndvi[i] <- aux 
# }

mean_ndvi <- vector("list", length = length(tabela1))
for (i in 1:length(tabela1)) {
  aux <- mean(pluck(tabela1,i)[3] %>% remove_missing)
  mean_ndvi[i] <- aux 
}


sd_ndvi <- vector("list", length = length(data2$tabela))
for (i in 1:length(data2$tabela)) {
  aux <- data2$tabela[[i]][3] %>%
    unlist() %>%
    as.numeric %>% 
    sd(na.rm=T)
  sd_ndvi[i] <- aux 
}

data1 <- data2 %>% 
  add_column(dados1) %>%  
  add_column(mean_ndvi = mean_ndvi %>% unlist) %>%
  add_column(sd_ndvi = sd_ndvi %>% unlist) %>% 
  add_column(lags) 

# Maior coupling dentre todos os lags
z <- data1[,29:34]
for (i in 1:dim(z)[1]){
  aux <- z[i,order(-abs(z[i,]))[1]]
  if (i == 1){
    coup_max <- aux
  }else{
    coup_max <- c(coup_max, aux)
  }
}
coup_max <- coup_max %>% unlist

# Lag at maximum coupling
for (i in 1:dim(z)[1]){
  aux <- order(-abs(z[i,]))[1]
  if (i == 1){
    lag <- aux
  }else{
    lag <- c(lag, aux)
  }
}
# data1 <- data1 %>% mutate(Domain = dexter$Domain)
# data1$Domain %>% unique

# Adiciona o maior acoplamento e o lag em que ele ocorreu
data1 <- data1 %>%
  add_column(coupling = coup_max) %>% 
  add_column(lag_max = lag)

setwd("/home/marcio/Documentos/Doutorado/Dexter/sentinel-2/")
# Muda o nome da coluna que identifica o plot
dados_pontos <- read.table("mapbiomas_mudanças_dexter.txt", h = T) %>% 
  remove_missing() %>% 
  rename(plot_id = ID) %>% 
  dplyr::select(plot_id | brasil_coverage_2021 | transition)
dados_pontos%>%head
# Junta transição mapbiomas com o resto dos dados
dados <- left_join(data1, dados_pontos, by = "plot_id")

head(dados)
# Prepara o VPD
vpd <- read_csv("/home/marcio/Documentos/Doutorado/Dexter/sentinel-2/vpd_do_Dexter_all.csv") %>%
  rename(vpd = value) %>% 
  rename(time = date) %>% 
  mutate(time = as.Date(time, format = "%Y-%m-%d")) %>% 
  group_by(plot_id) %>%
  nest

# Junta o VPD no resto dos dados

completo <- left_join(dados, vpd, by = "plot_id") %>% rename(vpd = data.y)

# Calcula a média do vpd
mean_vpd <- vector("list", length = dim(completo)[1])
for (i in 1:dim(completo)[1]) {
  aux <- purrr::map(pluck(completo$vpd,i)[2] %>% remove_missing, mean)
  mean_vpd[i] <- aux 
}

# Calcula o desvio do vpd
sd_vpd <- vector("list", length = dim(completo)[1])
for (i in 1:dim(completo)[1]) {
  aux <- completo$vpd[[i]][2] %>%
    unlist() %>%
    as.numeric %>% 
    sd(na.rm=T)
  sd_vpd[i] <- aux 
}

completo <- completo %>%
  add_column(mean_vpd) %>% 
  add_column(sd_vpd)

a1 <- completo %>% ggplot(aes(x = mean_vpd %>% as.numeric, y = coupling))+
  geom_point()+
  xlab("mean VPD")
a2 <- completo %>% ggplot(aes(x = sd_vpd %>% as.numeric, y = coupling))+
  geom_point()+
  xlab("sd VPD")
a3 <- completo %>% ggplot(aes(x = PrecAnn %>% as.numeric, y = coupling))+
  geom_point()+
  xlab("MAP")
a4 <- completo %>% ggplot(aes(x = CWD %>% as.numeric, y = coupling))+
  geom_point()+
  xlab("CWD")

(a1|a2)/(a3|a4)


b1 <- completo %>% ggplot(aes(x = mean_vpd %>% as.numeric, y = mean_ndvi))+
  geom_point()+
  xlab("mean VPD")+
  ylab("mean EVI2")
b2 <- completo %>% ggplot(aes(x = sd_vpd %>% as.numeric, y = mean_ndvi))+
  geom_point()+
  xlab("sd VPD")+
  ylab("mean EVI2")
b3 <- completo %>% ggplot(aes(x = mean_vpd %>% as.numeric, y = mean_ndvi))+
  geom_point()+
  xlab("mean VPD")+
  ylab("mean EVI2")
b4 <- completo %>% ggplot(aes(x = sd_vpd %>% as.numeric, y = mean_ndvi))+
  geom_point()+
  xlab("sd VPD")+
  ylab("mean EVI2")

(b1|b2)/(b3|b4)


# 14/02/2023
# filtering out what is not natural according to mapbiomas class
dev.off()
# exploring========================================
plot(completo[which(is.na(completo$transition)), ]$lon, completo[which(is.na(completo$transition)), ]$lat)

plot(completo[which(completo$transition == "Changed"), ]$lon,
     completo[which(completo$transition == "Changed"), ]$lat)

boxplot(completo[which(completo$transition == "Changed"), ]$PrecAnn~completo[which(completo$transition == "Changed"), ]$Domain)

boxplot(completo[which(completo$transition == "Not changed"), ]$PrecAnn~completo[which(completo$transition == "Not changed"), ]$Domain)

#==================================================

# by not changed
completoR <- completo[which(completo$transition == "Not changed"),]
#by natural veg
completoRNW <- completoR[which(completoR$brasil_coverage_2021 == 3 | 
                                completoR$brasil_coverage_2021 == 4 | 
                                completoR$brasil_coverage_2021 == 11 | 
                                completoR$brasil_coverage_2021 == 12), ]

completoRN <- completoR[which(completoR$brasil_coverage_2021 == 3 | 
                                completoR$brasil_coverage_2021 == 4 | 
                                completoR$brasil_coverage_2021 == 12), ]

# looking at the EVI2 time series (Sentinel and Landsat)===========
# by domain (Amazonia, Atlantic Forest, Caatinga, Cerrado, Gran Chaco)
dom <- "Amazonia"
dom <- "Cerrado" # not working - problem in leaf flush
dom <- "Caatinga" # not working - problem in leaf flush

exp <- completoRN[which(completoRN$Domain == dom), ]

boxplot(completoRN$mean_ndvi~completoRN$Domain, ylab = "mean EVI2", xlab = "Domain")

i <- 3
plot(exp$dados1[[i]]$evi2)
lines(exp$dados1[[i]]$without_missing$evi2, col = "red")
points(exp$dados1[[i]]$without_missing$evi2, col = "red")

plot(exp$dados1[[i]]$without_missing$evi2, type = "l")

# by leaf flush within the same domain+++++++++++++++++++++++
unique(exp$Leaf_flush)

leaf <- "Evergreen"
leaf <- "Semideciduous"
leaf <- "Deciduous"
exp.leaf <- exp[which(exp$Leaf_flush == leaf), ]
exp.leaf <- test

# ANNUAL CYCLE STORED FOR MEANS
#mon <- matrix(NA, nrow = dim(exp.leaf)[1], ncol = 13)
mon <- matrix(NA, nrow = dim(exp.leaf)[1], ncol = 12)
dev.off()
for (i in 1:dim(exp.leaf)[1]){
  print(i)
  # creating a df to plot the boxplots
  dfaux <- data.frame(exp.leaf$dados1[[i]]$date,
                      exp.leaf$dados1[[i]]$rainfall,
                      exp.leaf$dados1[[i]]$without_missing$evi2,
                      factor(substr(exp.leaf$dados1[[i]]$date, 6, 7)))
  colnames(dfaux) <- c("date", "rain", "evi2", "month")

    # subsetting by month and saving
  for (j in 1:length(levels(dfaux$month))){
    aux <- subset(dfaux, month == levels(dfaux$month)[j])
    mon[i,j] <- mean(aux$evi2)
    #mon[i,j+1] <- mean(aux$evi2)
    
  }
  #mon[i,1] <- exp.leaf$plot_id[i]
  
  # plotting line i
  if (i == 1){
    plot(mon[i,], type = "l", col = "lightgrey",
         ylim = c(0,1), xaxt = "n", xlab = "Months",
         ylab = "EVI2", main = paste0(leaf, " - ", dom))
  }else{
    lines(mon[i,], col = "lightgrey")
  }
}

# plotting the median value of each month
med <- apply(mon, 2, median)
lines(med, lwd = 2, col = "tomato")
axis(side = 1, at = 1:12, labels = c("Jan", "Feb", "Mar",
                                     "Apr", "May", "Jun",
                                     "Jul", "Aug", "Sep",
                                     "Oct", "Nov", "Dec") )
grid()
dev.off()
# boxplots per monthly median
mon.df <- as.data.frame(mon)
colnames(mon.df) <- c("Jan", "Feb", "Mar",
                      "Apr", "May", "Jun",
                      "Jul", "Aug", "Sep",
                      "Oct", "Nov", "Dec")

boxplot(mon.df, ylim = c(0,1), outline = T,
        xlab = "Months", ylab = "EVI2")
# find outliers within the data and in the map!!!!! they are 
# likely to be outside the evergreen zone.

# get mean and Standard deviation
mean = mean(exp.leaf$mean_ndvi)
std = sd(exp.leaf$mean_ndvi)

# get threshold values for outliers
EVI2min = mean-(3*std) %>% data.frame
EVI2max = mean+(3*std) %>% data.frame

# find outlier
exp.leaf$plot_id[which(exp.leaf$mean_ndvi < EVI2min$. | exp.leaf$mean_ndvi > EVI2max$.)]

# remove outlier
exp.leaf$plot_id[which(exp.leaf$mean_ndvi > EVI2min$. & exp.leaf$mean_ndvi < EVI2max$.)]
(dim(exp.leaf)[1])*0.5

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
# PLOTTING ANNUAL CYCLE OF PRECIP AND EVI FOR EACH LOCATION
for (i in 1:dim(exp.leaf)[1]){ #dim(exp.leaf)[1]){
  print(i)
  # creating a df to plot the boxplots
  dfaux <- data.frame(exp.leaf$dados1[[i]]$date,
                      exp.leaf$dados1[[i]]$rainfall,
                      exp.leaf$dados1[[i]]$without_missing$evi2,
                      factor(substr(exp.leaf$dados1[[i]]$date, 6, 7)))
  colnames(dfaux) <- c("date", "rain", "evi2", "month")
  
  # plotting boxplot for rain and evi2
  nf <- layout(matrix(1:2, ncol = 1, byrow = T))
  #layout.show(nf)
  
  par(mar = c(0, 5, 1, 2))
  boxplot(dfaux$rain~dfaux$month, ylim = c(0,max(dfaux$rain)),
          xaxt='n', ylab = "Rainfall (mm/month)", main = paste0(leaf, " - ", dom))
  grid()
  abline(h = 100, col = "tomato", lty = 2, lwd = 2)
  par(mar = c(4, 5, 0, 2))
  boxplot(dfaux$evi2~dfaux$month, ylim = c(0,1), ylab = "EVI2",
          xaxt = "n", xlab = "Months")
  axis(side = 1, at = 1:12, labels = c("Jan", "Feb", "Mar",
                                       "Apr", "May", "Jun",
                                       "Jul", "Aug", "Sep",
                                       "Oct", "Nov", "Dec"))
  grid()
}

#exp.leaf <- exp.leaf[-84,]

dev.off()
exp.leaf %>% dim

plot(ChirpsR, 1); points(acoords, col = "red")
#================

amo <- sample(1:106, 11)


exp.leaf %>% dim
test <- exp.leaf[amo,]

completoRN$sd_ndvi
completoRN$Vegetation.physiognomy %>% unique

a <-  completoRN %>%
  ggplot(aes(y = mean_ndvi, x = Domain, color = factor(Leaf_flush)))+
  geom_boxplot(show.legend = F)+
   ggtitle("a")+
   ylab("mean EVI2")+
  theme(axis.text=element_text(size=12),
        axis.text.x = element_text(size = 12, angle =90, vjust = 0.5, hjust=1),
        text=element_text(size=11),
        axis.title=element_text(size=13),
        title=element_text(size=14),
        legend.text=element_text(size=12))
  

b <- completoRN %>%
  ggplot(aes(y = sd_ndvi, x = Domain, color = factor(Leaf_flush)))+
  geom_boxplot()+
  ggtitle("b")+
  ylab("sd EVI2")+
  labs(color = "Vegetation type")+
  theme(axis.text=element_text(size=12),
        axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust=1),
        text=element_text(size=11),
        axis.title=element_text(size=13),
        title=element_text(size=14),
        legend.text=element_text(size=12))
getwd()
png("mean_e_sd_evi2_by_leaf_flush.png", res = 300, width = 2800, height = 1800)
a|b
dev.off()



?lm
plot(completoRN$mean_ndvi ~ model$fitted.values)
t1 <- completoRN %>%
  dplyr::group_by(Domain) %>% 
  dplyr::summarise(mean_evi2 = mean(mean_ndvi))

t2 <- completoRN %>%
  group_by(Domain) %>% 
  dplyr::summarise(sd_evi2 = mean(sd_ndvi))

completoRN %>% group_by(Domain) %>% 
  dplyr::summarise(min_evi2 = min(mean_ndvi))

left_join(t1, t2, by = "Domain")
