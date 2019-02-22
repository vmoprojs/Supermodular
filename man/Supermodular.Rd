\name{Supermodular}
\alias{Supermodular}
\title{Supermodular}
\usage{
Supermodular(ff= NULL,nboot = 10,semilla = NULL,DF = NULL)
}
\description{
Computes values of a dataframe that are supermodular
}


\arguments{
  \item{ff}{String; a string that represents the formula to be evauated in the model}
  \item{nboot}{Numeric; number of bootstrap iterations to obtain the probability distribution of the coefficients.}
  \item{semilla}{Numeric; seed value for the bootstrap in nboot parameter.}
  \item{DF}{data.frame; a dataframe containing the variables in the model.}
}


\author{Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/}
}

\examples{
# set.seed(2539)
# ENAIsam <- ENAI[sample(1:nrow(ENAI),500),]
# save(ENAIsam, file = "ENAIsam.RData")
rm(list = ls())
library(Supermodular)
data(ENAIsam)

# ******** START: ROUTINE INFORMATION
DF <- ENAIsam
rownames(DF) <- DF$id_empresa
ff <- "iventas_prod_mercado~A3 + A5 + gr_extr + exportaprom + creada +
protect_inn + id + otros + lnum_emp_prom + fuentes_externas +
sector_manufactura + sector_servicios + sector_comercio"
nboot <- 10
semilla <- 123

# ******** END: ROUTINE INFORMATION
# ff= ff;nboot = nboot;semilla = semilla;DF = DF
Supermodular(ff= ff,nboot = nboot,semilla = semilla,DF = DF)


}
