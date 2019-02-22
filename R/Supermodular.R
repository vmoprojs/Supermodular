# Supermodular!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
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


Supermodular <- function(ff= NULL,nboot = 10,semilla = NULL,DF = NULL)
{
  if(is.null(ff)) stop("A formula must be passed to funtion to ff argument")
  if(is.null(DF)) stop("A dataframe must be passed to DF argument")

  # ST: Formula treatmet to include interactions:
  ff <- fm.treat.model(ff,data = DF)
  # END: Formula treatmet to include interactions:


  #************************** ST: Step 3: Bootstrap
  sol <- superBoot(nboot,semilla,DF,f.lm=ff)

  #************************** END: Step 3: Bootstrap

  #************************** ST: Step 4: X-MEN
  # Varible definition
  variables <- colnames(sol)
  super.vars <- all.vars(ff)[2:3]
  ddd <- (model.matrix((ff),model.frame(ff,DF)))
  # Remove Duplicates:
  ddd <-ddd[!duplicated(ddd),]
  # *** Creating Spaces
  # (2^17)*(4^3)*(3^1)

  # Position of supermodular vars in model data matrix:
  pos.int <- match(super.vars,colnames(ddd),nomatch = 0)

  if(pos.int[1] ==0 | pos.int[2] ==0)
  {
    # super.sol[[cbn]] <- -1
    stop("Supermodular variables are not in the model data matrix")
  }else
  {
    # ST: Control Space:
    s.control <- ddd[,-c(1,pos.int)]
    if(is.vector(s.control))
    {
      s.control <- matrix(s.control,ncol = 1)
    }
    # END: Control Space

    # ST: Target Space:
    d.aux <- ddd[,c(pos.int)]
    if(is.null(ncol(d.aux)))
    {
      unique.aux <- unique(d.aux)
    }else
    {
      unique.aux <- apply(d.aux,2,unique)
    }
    s.int <- expand.grid(unique.aux[,1],unique.aux[,2])
    # END: Target Space.

    # plot(s.control)
    # plot(s.int)

    # ST: Looping through spaces
    res <- superBetas(s.control,s.int,nboot,sol)
    # END: Looping through spaces

    #************************** ST: Step 4: X-MEN

    #************************** ST: Step 4.1: Stochastic Dominace (Supermodularidad)
    # Application of stochastic dominace (SDom) to the control space
    # Supermodularidad: If there is SDom in all combinations of the target variables
    # Then, there is supermodularity in the control space

    SUPER <- superSochastic(res,s.control,s.int)
    respuesta <- list(values = ddd[SUPER,],ids = rownames(ddd)[SUPER]) # Control values where there is supermodularity
    #************************** ST: Step 4.1: Stochastic Dominace (Supermodularidad)
  }
  #************************** END: Step 4: X-MEN
  return(respuesta)
}
