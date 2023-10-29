# Universidad de Valparaíso
# Trabajo práctico 3
# IECD323: Estadística Bayesiana
# Alumno: Bastián Barraza Morales
# Profesor: Mauricio Tejo Arriagada

setwd('C:\\Users\\Bastian Barraza M\\OneDrive\\Documentos\\Semestre 8\\Bayesian Statistics\\Trabajo práctico3')

library(ggplot2)
library(dplyr)
library(plotly)
library(gridExtra)
library(pastecs)
library(plot3D)

##################################
######## EJERCISIO 1 #############
##################################

set.seed(2)
n=50
theta = 5
phi1 = 2
SIM50 = data.frame('simulation' = rnorm(n, 
                                theta, 
                                sd = sqrt(phi1)))
stat.desc(SIM50$simulation)
summary(SIM50$simulation)

# X Graph
ggplot(SIM50, aes(x=simulation)) + 
 geom_histogram(aes(y=..density..), 
                colour=1, 
                fill="white",
                bins = 11)+
 geom_density(lwd = 1,
                alpha=.15, 
                fill="#FF6666")  +
 geom_vline(aes(xintercept=mean(simulation)),
            color="blue", size=1)  + 
    xlab("X") + 
    ylab("Probabilidad")+ 
    ggtitle('Simulación desde una distribución N(5, 2)') +
    theme(text = element_text(size = 20)) 
ggsave("P1_XSimulated.png")


##################################
######## EJERCISIO 2 #############
##################################

# hiperparametros
mu0 <- 0
t20 <- 10
s20 <-  25
nu0 <- 10

##################################
######## EJERCISIO 3 #############
##################################

# numero de muestras 
S <- 1000

# matriz para almacenar las muestras
PHI <- matrix(data = NA, nrow = S, ncol = 2)
# ALGORITMO (muestreador de Gibbs)
# 1. inicializar la cadena
#    valor inicial: simular de la previa
#    solo es necesario alguno de los valores
set.seed(1)
phi <- c( rnorm(1, mu0, sd = sqrt(t20)), rgamma(1, nu0/2, s20/2) )
PHI[1,] <- phi

# 2. simular iterativamente de las distribuciones condicionales completas
set.seed(1)
for(s in 2:S) {
        # 2.1 actualizar el valor de theta
        t2n    <- 1/( (t20**(-1)) + n*phi[2] )      
        mun    <- (n*phi[2]*mean(SIM50$simulation)) / (t20**(-1) + n*phi[2])
        phi[1] <- rnorm( 1, mun, sd = sqrt(t2n) )
        # 2.2 actualizar el valor de sigma^2
        nun    <- nu0+n
        s2n    <- s20 + (n-1)*var(SIM50$simulation) + n*(phi[1] - mean(SIM50$simulation))**2
        phi[2] <- rgamma(1, nun/2, s2n/2)
        # 2.3 almacenar
        PHI[s,] <- phi
}

# Plotting theta.
png(file="P3_theta.png",
    width=1000, height=1000)
plot(PHI[2:1000,1],type="l",ylab="Theta",xlab="Iteration")
dev.off()

#Plotting PHI
png(file="P3_phi.png",
    width=1000, height=1000)
plot(PHI[,2],type="l",ylab="Phi",xlab="Iteration")
dev.off()


##################################
######## EJERCISIO 4 #############
##################################

sub_phi = data.frame('post_theta' = PHI[501:1000,1],
                    'post_sigma' = PHI[501:1000,2])
sub_phi

# Theta
mean(sub_phi$post_theta)
var(sub_phi$post_theta)

# PHI
mean(sub_phi$post_sigma)
var(sub_phi$post_sigma)

##################################
######## EJERCISIO 5 #############
##################################

# THETA GRAPH
ggplot(sub_phi, aes(x=post_theta)) + 
 geom_histogram(aes(y=..density..), 
                colour=1, 
                fill="white",
                bins = 15)+
 geom_density(lwd = 1,
                alpha=.15, 
                fill="#FF6666")  +
 geom_vline(aes(xintercept=mean(post_theta)),
            color="blue", size=1)  + 
    xlab("Theta") + 
    ylab("Densidad")+ 
    ggtitle('Simulación de distribución condicional a posteriorir de Theta\n usando las últimas 500 repeticiones.') +
    theme(text = element_text(size = 20))    

ggsave("P5_ThetaSimulated.png")

# PHI GRAPH
ggplot(sub_phi, aes(x=post_sigma)) + 
 geom_histogram(aes(y=..density..), 
                colour=1, 
                fill="white",
                bins = 15)+
 geom_density(lwd = 1,
                alpha=.15, 
                fill="#FF6666")  +
 geom_vline(aes(xintercept=mean(post_sigma)),
            color="blue", size=1)  + 
    xlab("Phi") + 
    ylab("Densidad")+ 
    ggtitle('Simulación de distribución condicional a posteriorir de Phi\n usando las últimas 500 repeticiones.') +
    theme(text = element_text(size = 20))    

ggsave("P5_PhiSimulated.png")


##################################
######## EJERCISIO 6 #############
##################################

theta = cut(sub_phi$post_theta, 25)
sigma = cut(sub_phi$post_sigma, 25)
x = table(theta, sigma)
z = matrix(x, nrow=25, ncol=25, byrow=TRUE)

#Plot distribución conjunta método gibbs
fig <- plot_ly(z = ~z) %>% add_surface(
  contours = list(
    z = list(
      show=TRUE,
      usecolormap = TRUE,
      highlightcolor="#ff0000",
      project=list(z=TRUE)
      )
    )
  )
fig <- fig %>% layout(
    scene = list(camera=list(eye = list(x=1.87, y=0.88, z=-0.64)),
              xaxis = list(ticketmode = 'array',
                          ticktext = seq(4.28, 5.61, 0.05),
                          tickvals = 0:24,
                          range = c(0,24)),

              yaxis = list(ticketmode = 'array',
                          ticktext = round(seq(0.283,0.966, 0.027), 2),
                            tickvals = 0:24,
                            range = c(0,24))
      )
  )

fig


##Generar las muestras mediante aceptacion / rechazo

set.seed(2)
rep = 500
PHI <- matrix(data = NA, nrow = rep, ncol = 2)
phi <- c( rnorm(1, mu0, sd = sqrt(t20)), rgamma(1, nu0/2, s20/2) )
PHI[1,] <- phi

for(s in 2:rep) {
        # 2.1 actualizar el valor de theta
        t2n    <- 1/( (t20**(-1)) + n*phi[2] )      
        mun    <- (n*phi[2]*mean(SIM50$simulation)) / (t20**(-1) + n*phi[2])
        phi[1] <- rnorm( 1, mun, sd = sqrt(t2n) )
        alpha = dnorm(phi[1], mean=5, sd=sqrt(2)) / dnorm(PHI[s-1,1], mean=5, sqrt(2))
        
        if (runif(1) < alpha) { 
          PHI[s, 1] = phi[1]}

        else{
          PHI[s] <- PHI[s - 1, 1]}

        #2.2 actualizar el valor de sigma^2
        
        nun    <- nu0+n
        s2n    <- s20 + (n-1)*var(SIM50$simulation) + n*(phi[1] - mean(SIM50$simulation))**2
        phi[2] <- rgamma(1, nun/2, s2n/2)
        alpha = dgamma(phi[2], shape = nun/2, scale = s2n/2) / dgamma(PHI[s-1,2], shape = nun/2, scale = s2n/2)
        # 2.3 almacenar
        if(runif(1) < alpha){
            PHI[s, 2] <- phi[2]}
        else{
            PHI[s, 2] <- PHI[s - 1, 2]}
}

PHI
hist(PHI[2:500,2])

# Graficar mediante aceptacion / rechazo
theta2 = cut(PHI[2:500,1], 25)
sigma2 = cut(PHI[2:500,2], 25)
x = table(theta2, sigma2)
z = matrix(x, nrow=25, ncol=25, byrow=TRUE)

fig <- plot_ly(z = ~z) %>% add_surface(
  contours = list(
    z = list(
      show=TRUE,
      usecolormap = TRUE,
      highlightcolor="#ff0000",
      project=list(z=TRUE)
      )
    )
  )
fig <- fig %>% layout(
    scene = list(camera=list(eye = list(x=1.87, y=0.88, z=-0.64)),
              xaxis = list(ticketmode = 'array',
                          ticktext = seq(4.28, 5.61, 0.05),
                          tickvals = 0:24,
                          range = c(0,24)),

              yaxis = list(ticketmode = 'array',
                          ticktext = round(seq(0.283,0.966, 0.027), 2),
                            tickvals = 0:24,
                            range = c(0,24))
      )
  )

fig

