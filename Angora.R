#EVALUACI?N GENETICA PARA CAPRINOS DE ANGORA

# Pendientes
#- INDICE
#- Agregar la parte visual y observaciones al archivo final de DEPs

#### HISTORICO  ####
#2018
#-Arranca nuevamente la Evaluaci?n Gen?tica del Nucleo de Angora

#2019
#-Primera evaluaci?n a Campana

#2020
#-Paso de SAS a R

#COMIENZO EVALUACION

# Parametros --------------------------------------------------------------
#Creo directorios de trAbajo
dir.create("c:/Drive/DATA/ANGORA/POB/2023/")                  #Actualizar
dir.create("c:/Drive/DATA/ANGORA/POB/2023/BLUP1")             #Actualizar
#dir.create("c:/Drive/DATA/ANGORA/POB/2021/BLUP2")             #Actualizar
#dir.create("c:/Drive/DATA/ANGORA/POB/2021/BLUP3")             #Actualizar

procap <- ("c:/Drive/DATA/ANGORA/POB/2023")                   #Actualizar
BLUP1 <- ("c:/Drive/DATA/ANGORA/POB/2023/BLUP1")              #Actualizar
#BLUP2 <- ("c:/Drive/DATA/ANGORA/POB/2021/BLUP2")              #Actualizar
#BLUP3 <- ("c:/Drive/DATA/ANGORA/POB/2021/BLUP3")              #Actualizar

ult.camada <- 2021
est <- 'PIL'
limgen <- 5
per.int <- 2015

library(ggplot2)
datos <- read.csv("c:/Drive/DATA/ANGORA/TEMP/BD Angora.csv", sep=";")

datos$FN<-as.Date(paste(datos$DN,datos$MN,datos$AN,sep='/'),format="%d/%m/%Y")     #Nacimiento
datos$F50<-as.Date(paste(datos$DP50,datos$MP50,datos$AP50,sep='/'),format="%d/%m/%Y") #PreDestete
datos$F100<-as.Date(paste(datos$DP100,datos$MP100,datos$AP100,sep='/'),format="%d/%m/%Y")     #Destete
datos$FE6<-as.Date(paste(datos$DE6,datos$ME6,datos$AE6,sep='/'),format="%d/%m/%Y")     #Esquila6
datos$FE12<-as.Date(paste(datos$DE12,datos$ME12,datos$AE12,sep='/'),format="%d/%m/%Y")     #Esquila12

datos$EDM <- as.numeric(datos$AN-datos$ANM) #Edad de la madre;
datos$ED50 <- as.numeric(datos$F50-datos$FN) #Edad a la se?alada;
datos$ED100 <- as.numeric(datos$F100-datos$FN)  #Edad al destete;
datos$EDE6 <- as.numeric(datos$FE6-datos$FN)  #Edad a la esquila 6;
datos$EDE12 <- as.numeric(datos$FE12-datos$FN)  #Edad a la esquila 12;

datos$ED50 <- ifelse(datos$ED50<20|datos$ED50>80,NA,datos$ED50)
datos$ED100 <- ifelse(datos$ED100<60|datos$ED100>140,NA,datos$ED100)
datos$EDE6 <- ifelse(datos$EDE6>160,NA,datos$EDE6) #Estos son unos animales de Pilca que se desbarrigaron a destiempo



#Descriptivos
vars <- datos[,c("EDM","ED50","ED100","EDE6","EDE12")] #ME quedo con las variables de interes
n <- colSums(!is.na(vars))
mean <- apply(vars, 2, mean, na.rm=TRUE)
median <- apply(vars, 2, median, na.rm=TRUE)
max <- apply(vars, 2, max, na.rm=TRUE)
min <- apply(vars, 2, min, na.rm=TRUE)
tabla1 <- rbind(n,mean,median,max,min)
tabla1

rm(vars,n,mean,median,max,min)

vars <- datos[,c("PCN","PC50","PC100","PCE6","PCE12")] #ME quedo con las variables de interes
n <- colSums(!is.na(vars))
mean <- apply(vars, 2, mean, na.rm=TRUE)
sd <- apply(vars, 2, sd, na.rm=TRUE)
max <- apply(vars, 2, max, na.rm=TRUE)
min <- apply(vars, 2, min, na.rm=TRUE)
tabla2 <- rbind(n,mean,sd,max,min)
tabla2
rm(vars,n,mean,sd,max,min)

vars <- datos[,c("PVS6","PDF6","PVS12","PDF12","KEMP12")] #ME quedo con las variables de interes
n <- colSums(!is.na(vars))
mean <- apply(vars, 2, mean, na.rm=TRUE)
sd <- apply(vars, 2, sd, na.rm=TRUE)
max <- apply(vars, 2, max, na.rm=TRUE)
min <- apply(vars, 2, min, na.rm=TRUE)
tabla2 <- rbind(n,mean,sd,max,min)
tabla2
rm(vars,n,mean,sd,max,min)

#Recodifico variables
datos$TIPO <- ifelse(datos$TIPO==3,2,datos$TIPO)

datos$CATEDM <- ifelse(datos$EDM==2,2,datos$EDM)
datos$CATEDM <- ifelse(datos$CATEDM>5,610,datos$CATEDM)
datos$CATEDM <- ifelse(datos$CATEDM>2 & datos$CATEDM<=5|is.na(datos$CATEDM),35,datos$CATEDM)

datos$GC <- paste(datos$CodCampo,datos$AN,datos$SEXO,sep="")


#Crea el archivo de pedigree
ped <- subset(datos,!is.na(datos$PCE12)|!is.na(datos$PVS12)|!is.na(datos$PDF12))
ped <- ped[,c("IDUI","IDUP","IDUM")]

#Crea los archivos de datos
dat1 <- datos[,c("IDUI","GC","TIPO","CATEDM","EDE12","PCE12","PVS12","PDF12")]
dat1 <- subset(dat1,!is.na(dat1$PCE12)|!is.na(dat1$PVS12)|!is.na(dat1$PDF12))

#Exporta archivo de pedigree y de datos
setwd(BLUP1)
write.table(ped, file = "ped.txt", sep=" ", row.names = F, col.names = F, quote=F, na="0")
write.table(dat1, file = "dat.txt", sep=" ", row.names = F, col.names = F, quote=F, na="0")

#Corre BLUP
#antes copiar renum.par de la corrida del a?o anterior a la carpeta (BLUP1)
shell("renumf90 renum.par")
shell("blupf90 renf90.par")

#Lee BLUP1
VC1 <-read.delim("solutions", skip=1, header=F, sep='')
names(VC1)=c("trait","effect","IDUI.r","solution","se")
VC1 <- subset(VC1,VC1$effect==5)
VC1 <- VC1[,-2]

VC1 <- reshape(VC1, idvar = "IDUI.r", timevar = "trait", direction = "wide")
names(VC1)=c("IDUI.r","vcPCE12","sePCE12","vcPVS12","sePVS12","vcPDF12","sePDF12")

inb1 <-read.delim("renf90.inb", header=F,sep='')
names(inb1)=c("IDUI","inb1","code")
inb1 <- inb1[,-3]

ped.r <-read.delim("renadd05.ped", header=F,sep='')
names(ped.r)=c("IDUI.r","IDUP.r","IDUM.r","v4","V5","V6",'V7',"n.pad","n.mad","IDUI")
ped.r <-ped.r[,c(1,2,3,8:10)]

VC1.j <- merge(ped.r,inb1, by='IDUI')
VC1.j <- merge(VC1,VC1.j, by='IDUI.r')

VC1.j$acPCE12 <- (sqrt(1-(VC1.j$sePCE12*VC1.j$sePCE12)/(4.505*(1+VC1.j$inb/100))))*100
VC1.j$acPVS12 <- (sqrt(1-(VC1.j$sePVS12*VC1.j$sePVS12)/(0.031*(1+VC1.j$inb/100))))*100
VC1.j$acPDF12 <- (sqrt(1-(VC1.j$sePDF12*VC1.j$sePDF12)/(1.566*(1+VC1.j$inb/100))))*100

VC1.j$exa <- apply(VC1.j[ ,c("acPCE12","acPVS12","acPDF12")], 1, mean, na.rm = TRUE)
#Me quedo con las variables de interes para el informe
VC1.j <-VC1.j[,-c(1,3,5,7,9:12)]

#Calcula los deps estandarizados
VC1.j$IDUI <- as.factor(VC1.j$IDUI)
todo.j <- merge(datos, VC1.j, by='IDUI',all.y=T)
todo.j$AN <- as.factor(todo.j$AN)

#A?o Base
ab <- subset(todo.j,todo.j$AN=='2015')

todo.j$vc.PCE12 <- todo.j$vcPCE12 - mean(ab$vcPCE12, na.rm=T)
todo.j$vc.PVS12 <- todo.j$vcPVS12 - mean(ab$vcPVS12, na.rm=T)
todo.j$vc.PDF12 <- todo.j$vcPDF12 - mean(ab$vcPDF12, na.rm=T)

todo.j$depPCE12 <- todo.j$vc.PCE12/2
todo.j$depPVS12 <- todo.j$vc.PVS12/2
todo.j$depPDF12 <- todo.j$vc.PDF12/2

rm(inb1,ped.r)

#configuro decimales y renombro variables

todo.j2 <- todo.j
todo.j2$con <- todo.j2$inb*100

todo.j2$PDF12 <- round(todo.j2$PDF12,digits = 2)
todo.j2$PCE12 <- round(todo.j2$PCE12,digits = 2)
todo.j2$PVS12 <- round(todo.j2$PVS12,digits = 2)

todo.j2$depPDF12 <- round(todo.j2$depPDF12,digits = 2)
todo.j2$depPCE12 <- round(todo.j2$depPCE12,digits = 2)
todo.j2$depPVS12 <- round(todo.j2$depPVS12,digits = 2)
todo.j2$exa <- round(todo.j2$exa,digits = 0)

setwd(procap)
#Listados de PADRES por campo
as.character(todo.j2$AN)
#esto es para el contador de hijos en el periodo de interes
cpdf1 <- subset(todo.j2,as.character(todo.j2$AN)>=per.int)

#aqui tabula padres x el numero de hijos
tx <- as.data.frame.matrix(table(cpdf1$IDUP,cpdf1$CodCampo))
txx <- as.data.frame(table(cpdf1$IDUP,cpdf1$CodCampo))
colnames(txx)[colnames(txx)=="Var1"] <- "IDUP"
colnames(txx)[colnames(txx)=="Var2"] <- "CAMPO.p"
txx <- subset(txx,txx$Freq!=0)
txx <- txx[!duplicated(txx$IDUP),]

#tx$NT <- tx$CAM+tx$DES+tx$MPA+tx$PHZ+tx$PIL 
tx$NT <-tx$PIL #Por ahora no hay otro campo que Pilca desde 2015
tx$IDUP <- rownames(tx)
tx1 <- tx[,c('NT','IDUP')]

#aqui quedan los padres con el numero de campos con progenie
#Por ahora hay un solo campo por lo cual NC=1
tx1$NC <- 1
tx3 <-tx1

#tx2 <- tx
#tx2[tx2!=0] <-1
#tx2$NC <- tx2$CAM+tx2$DES+tx2$MPA+tx2$PHZ+tx2$PIL 
#tx2
#tx2$IDUP <- rownames(tx2)
#tx2 <- tx2[,c('IDUP','NC')]

#tx3 <- merge(tx1,tx2,by='IDUP')
#tx3 <- merge(tx3,txx,by='IDUP')
tx3 <- subset(tx3,tx3$IDUP!=0)

rm(tx,txx,tx1,tx2)

colnames(tx3)[colnames(tx3)=="IDUP"] <- "IDUI"

anp <- cpdf1[,c('IDUP',"ANPad","PADRE")]
anp <- anp[!duplicated(anp$IDUP),]
anp <- anp[!is.na(anp$IDUP),]
names(anp)=c("IDUI","ANPad","PADRE")
anp$IDUI <- as.character(anp$IDUI)

padres <- merge(todo.j2,tx3,by='IDUI',all.y=T)
padres <- merge(padres,anp,by='IDUI',all.y=T)

#padres$CAMPO.p[is.na(padres$CAMPO.p)] <- padres$CodCampo

#padres <- subset(padres,padres$CodCampo==est)
pad.est <- padres[,c('PADRE.y','PADRE.x','ANPad.y','PCE12','PVS12','PDF12','depPCE12','depPVS12','depPDF12','exa','NT','NC','inb1')]
names(pad.est) <- c('ID','PADRE','ANPad','PCE12','PVS12','PDF12','depPCE12','depPVS12','depPDF12','exa','NT','NC','inb1')

pad.est <- pad.est[order(pad.est$ID, decreasing = T),]

write.table(pad.est, file = 'Padres.csv', sep=";", row.names = F, col.names = T, quote=T, na="0")

#Listado de Machos 2D de campo
progM2D.est <- subset(todo.j2,todo.j2$AN==ult.camada & todo.j2$CodCampo==est & !is.na(todo.j2$CodCampo) & todo.j2$SEXO=='M')
progM2D.est <- progM2D.est[,c('ID','PADRE','PCE12','PVS12','PDF12','KEMP12','depPCE12','depPVS12','depPDF12','exa','inb1',"ObsNacDest","ObsEsquila12")]
progM2D.est <- progM2D.est[order(progM2D.est$ID),]

write.table(progM2D.est, file ='Castroncitos 2D.csv', sep=";", row.names = F, col.names = T, quote=T, na="0")

#Listado de Hembras 2D de campo
progH2D.est <- subset(todo.j2,todo.j2$AN==ult.camada & todo.j2$CodCampo==est & !is.na(todo.j2$CodCampo) & todo.j2$SEXO=='H')
progH2D.est <- progH2D.est[,c('ID','PADRE','PCE12','PVS12','PDF12','KEMP12','depPCE12','depPVS12','depPDF12','exa','inb1',"ObsNacDest","ObsEsquila12")]
progH2D.est <- progH2D.est[order(progH2D.est$ID),]

write.table(progH2D.est, file ='Cabrillas 2D.csv', sep=";", row.names = F, col.names = T, quote=T, na="0")

#Listado de Machos 4D de campo
prog4D.est <- subset(todo.j2,todo.j2$AN==ult.camada-1 & todo.j2$CodCampo==est & !is.na(todo.j2$CodCampo) & todo.j2$SEXO=='M' )
prog4D.est <- prog4D.est[,c('ID','PADRE','PCE12','PVS12','PDF12','KEMP12','depPCE12','depPVS12','depPDF12','exa','inb1',"ObsNacDest","ObsEsquila12")]
prog4D.est <- prog4D.est[order(prog4D.est$ID),]

write.table(prog4D.est, file ='Castrones 4D.csv', sep=";", row.names = F, col.names = T, quote=T, na="0")

#Listado de Machos 6D de campo
prog6D.est <- subset(todo.j2,todo.j2$AN==ult.camada-2 & todo.j2$CodCampo==est & !is.na(todo.j2$CodCampo) & todo.j2$SEXO=='M' )
prog6D.est <- prog6D.est[,c('ID','PADRE','PCE12','PVS12','PDF12','KEMP12','depPCE12','depPVS12','depPDF12','exa','inb1',"ObsNacDest","ObsEsquila12")]
prog6D.est <- prog6D.est[order(prog6D.est$ID),]

write.table(prog6D.est, file ='Castrones 6D.csv', sep=";", row.names = F, col.names = T, quote=T, na="0")

#Listado de Machos 8D de campo
prog8D.est <- subset(todo.j2,todo.j2$AN==ult.camada-3 & todo.j2$CodCampo==est & !is.na(todo.j2$CodCampo) & todo.j2$SEXO=='M' )
prog8D.est <- prog8D.est[,c('ID','PADRE','PCE12','PVS12','PDF12','KEMP12','depPCE12','depPVS12','depPDF12','exa','inb1',"ObsNacDest","ObsEsquila12")]
prog8D.est <- prog8D.est[order(prog8D.est$ID),]

write.table(prog8D.est, file ='Castrones 8D.csv', sep=";", row.names = F, col.names = T, quote=T, na="0")

#Tablas
#Tabla 1: Estructura de datos

#datos$SEXO <- droplevels(datos$SEXO)
#datos$TIPO <- droplevels(datos$TIPO)
#datos$AN <- droplevels(datos$AN)


t1 <- table(datos$AN,datos$SEXO)
t2 <- table(datos$AN,datos$TIPO)
table(datos$SEXO)

tabla1 <- as.data.frame(cbind(t1,t2))
tabla1$total <- tabla1$H+tabla1$M
t14 <- (colSums(tabla1))

tabla1 <- rbind(tabla1,t14)
tabla1$AN <- row.names(tabla1)

tabla1$AN <- ifelse(tabla1$AN==nrow(tabla1),'Total',tabla1$AN)

names(tabla1)=c("SEX:h","SEX:m","Tipo:1","Tipo:2",'Total','AN')
tabla1 <- tabla1[,c(ncol(tabla1),ncol(tabla1)-1,1:(ncol(tabla1)-2))] #reordeno las columnas
write.table(tabla1, file = "Tabla 1.csv", sep=";", row.names = F, col.names = T, quote=T, na=".")

rm(t10,t11,t12,t13,t14,tabla1,tendencia,progen,pob1)


#Tablas
#Tabla 1: Estructura de datos PILCA


t1 <- table(datos$AN,datos$SEXO)
t2 <- table(datos$AN,datos$TIPO)
table(datos$SEXO)

tabla1 <- as.data.frame(cbind(t1,t2))
tabla1$total <- tabla1$H+tabla1$M
t14 <- (colSums(tabla1))

tabla1 <- rbind(tabla1,t14)
tabla1$AN <- row.names(tabla1)

tabla1$AN <- ifelse(tabla1$AN==nrow(tabla1),'Total',tabla1$AN)

names(tabla1)=c("SEX:h","SEX:m","Tipo:1","Tipo:2",'Total','AN')
tabla1 <- tabla1[,c(ncol(tabla1),ncol(tabla1)-1,1:(ncol(tabla1)-2))] #reordeno las columnas
write.table(tabla1, file = "Tabla 1.p.csv", sep=";", row.names = F, col.names = T, quote=T, na=".")

rm(t1,t2,t14,tabla1)


#Tabla 3: N? de cr?as por castr?n por a?o
#datos$IDUP <- droplevels(datos$IDUP)

tabla3 <- as.data.frame.matrix(table(cpdf1$PADRE,cpdf1$AN))
tabla3$Total <- rowSums(tabla3)

t15 <- (colSums(tabla3))
tabla3 <- rbind(tabla3,t15)
tabla3$PADRE <- row.names(tabla3)
tabla3[tabla3==0] <-NA
tabla3$PADRE <- ifelse(tabla3$PADRE==nrow(tabla3),'Total',tabla3$PADRE)

tabla3 <- subset(tabla3,tabla3$PADRE!='')
tabla3 <- tabla3[,c(ncol(tabla3),1:(ncol(tabla3)-1))] #reordeno las columnas

write.table(tabla3, file = "Tabla 3.csv", sep=";", row.names = F, col.names = T, quote=T, na=".")

#Tabla 3: N? de cr?as por castr?n por a?o PILCA

cpdf2 <- subset(cpdf1,cpdf1$CodCampo=='PIL')
tabla3 <- as.data.frame.matrix(table(cpdf2$PADRE,cpdf2$AN))
tabla3$Total <- rowSums(tabla3)

t15 <- (colSums(tabla3))
tabla3 <- rbind(tabla3,t15)
tabla3$PADRE <- row.names(tabla3)
tabla3[tabla3==0] <-NA
#tabla3$PADRE <- ifelse(tabla3$PADRE==nrow(tabla3),'Total',tabla3$PADRE)

tabla3 <- subset(tabla3,tabla3$PADRE!='')
tabla3 <- tabla3[,c(ncol(tabla3),1:(ncol(tabla3)-1))] #reordeno las columnas

write.table(tabla3, file = "Tabla 3.p.csv", sep=";", row.names = F, col.names = T, quote=T, na=".")

#Tabla 4. Valores de Cr?a (VC) promedio para el a?o base

vcPCE12 <- mean(ab$vcPCE12, na.rm=T)
vcPVS12 <- mean(ab$vcPVS12, na.rm=T)
vcPDF12 <- mean(ab$vcPDF12, na.rm=T)

m.ab <- cbind(vcPCE12,vcPVS12,vcPDF12)

write.table(m.ab, file = "Tabla 4.csv", sep=";", row.names = F, col.names = T, quote=T, na=".")

#Tendencias gen?ticas
tend <- subset(todo.j2, as.character(todo.j2$AN)>=per.int & todo.j2$CodCampo==est)  #aca defino que a?os incluir en la tend genetica
tend$AN <- droplevels(tend$AN)
t1 <- round(tapply(tend$vc.PCE12,tend$AN, mean, na.rm=T),digits = 2)
t2 <- round(tapply(tend$vc.PVS12,tend$AN, mean, na.rm=T),digits = 2)
t3 <- round(tapply(tend$vc.PDF12,tend$AN, mean, na.rm=T),digits = 2)

tendencia<- as.data.frame(cbind(t1,t2,t3))
tendencia$AN <- rownames(tendencia)
names(tendencia)=c("vcPC12",'vcPVS12',"vcPDF12",'AN')
tendencia <- tendencia[,c(ncol(tendencia),1:(ncol(tendencia)-1))] #reordeno las columnas

write.table(tendencia, file = 'Tendencias.csv', sep=";", row.names = F, col.names = T, quote=T)

tend.p <- reshape(tendencia, idvar = "AN", varying= list(2:4),v.names = 'TRAIT',times = c('vcPC12','vcPVS12','vcPDF12'), direction = "long")

png(file="Tendencias.png", width=1280, heigh=720 )
ggplot(data=tend.p, aes(AN,TRAIT, group=time, color=time)) + geom_path() + geom_point()+
  xlab("A?o de Nacimiento") + 
  ylab("Valores de Cria") + 
  theme_gray()
dev.off()

