# Utility Functions
#
# This is a function named 'row.match'
# which returns the row position that matches
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

row.match <- function(X,x.c)
{
  # X: Matrix
  # x.c: vector to be compared to each row of X
  #Sol: returns the row position that matches
  Sol <- apply(X, 1, (function(x) identical(as.numeric(x), x.c)) )
  return(which(Sol))
}



fm.treat.model <- function(ff,data = DF)
{
  ff <- formula(ff)

  super.vars <- all.vars(ff)[2:3]

  ff.aux <- paste(all.vars(ff)[-1],collapse ="+")
  ff.aux <- paste("(",ff.aux,")","**2")
  ff.aux <- formula(paste(c(all.vars(ff)[1],"~",ff.aux),collapse = ""))
  ff <- ff.aux


  # ST: Model Estimation with al interactions and using step to select the "best"
  bslmregA<-lm(ff,data=DF)
  # Final model
  f.lm <- step(bslmregA,trace = 0);formula(f.lm)
  # END: Model Estimation with al interactions and using step to select the "best"



  #***** ST: step 2: Checking if the final model
  # Checking if the final model contains the supermodular interaction
  fil.1 <- identical(all.vars(formula(f.lm))[2:3],super.vars)

  if(!fil.1)
  {
    stop("Supermodular variables are not in the final model.\nTry another model specification.")
  }

  ff <- formula(f.lm)
  return(ff)
  #***** END: step 2: Checking if the final model

}



superBoot <- function(nboot,semilla,DF,f.lm=ff)
{
  n.init <- nrow(DF)
  muestra <- 1:n.init
  if(!is.null(semilla)) {set.seed(semilla)}
  sol <- NULL
  boot.cont <- 0

  while(boot.cont <nboot)
  {
    sam.boot <- sample(muestra,n.init,replace = TRUE)
    dat.boot <- DF[sam.boot,]
    aux <- coef(lm((f.lm),data =dat.boot ))

    if(!any(is.na(as.numeric(aux))))
    {# It may happen that a "boot" selection gives lots of cero values, so
      # some "coef" be "NA". This "if" condition controls fot his issue
      sol <- rbind(sol,aux)
      boot.cont <- boot.cont+1
      cat("Boot.Iter:", boot.cont,"of ",nboot,"\n")
    }
  }
  return(sol)
}



superBetas <- function(s.control,s.int,nboot,sol)
{
  contador <- 1
  n.tot <- nrow(s.control)*nrow(s.int)*nboot
  res <- array(NA,dim = c(nrow(s.control),nrow(s.int),nboot))
  for(ctrl in 1:nrow(s.control))
  {
    for(itrs in 1:nrow(s.int))
    {
      for(boot in 1:nboot)
      {
        a <- as.numeric(sol[boot,])
        b <- c(1,as.numeric(s.int[itrs,]),as.numeric(s.control[ctrl,]))
        cat("Iteration:",contador,"of",n.tot,"z.hat",a%*%b,"\n")
        contador <- contador+1
        res[ctrl,itrs,boot] <- a%*%b
      }
    }
  }
  return(res)
}


superSochastic <- function(res,s.control,s.int)
{
  z.range <- range(res)
  z.vals <- seq(z.range[1],z.range[2],10)

  SUPER <- NULL
  SUPER.prop <- NULL
  ctrl <- 1; i <- 1; j <- 1
  for(ctrl in 1:nrow(s.control))
  {
    # logic vector that saves the SDom of the comparisons of the target variables
    c.dom.sol <- NULL
    for(i in 1:(nrow(s.int)-1))
    {
      for(j in (i+1):nrow(s.int))
      {
        # print(c(i,j))
        fx <- (ecdf(res[ctrl,i,]))
        fy <- (ecdf(res[ctrl,j,]))

        # Valores x,y in s.int
        x <- as.numeric(s.int[i,])
        y <- as.numeric(s.int[j,])
        # xy is the matrix to find the "arrow" relations
        xy <- matrix(c(x,y),ncol = 2)
        x.min.y <- apply(xy,1,min)
        x.max.y <- apply(xy,1,max)
        # ii,jj are the positions from which z-boot values are extracted
        ii <- row.match(s.int,x.min.y)
        jj <- row.match(s.int,x.max.y)

        fmin <- ecdf(res[ctrl,ii,])
        fmax <- ecdf(res[ctrl,jj,])

        # Stochastic dominance condition:
        c.dom <- (fmax(z.vals)+fmin(z.vals)) >= (fx(z.vals)+fy(z.vals))
        c.dom.sol <- c(c.dom.sol,all(c.dom))
      }
    }
    cat("Iter: ",ctrl,"of: ",nrow(s.control),"\n")
    SUPER <- c(SUPER,all(c.dom.sol))
    SUPER.prop <- c(SUPER.prop,as.numeric(prop.table(table(c.dom.sol))[2]))
  }
  return(SUPER)
}
