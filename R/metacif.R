## metacif.R contains functions to perform meta-analysis of summary data with competing risks.
## Copyright (C) 2019 Federico Bonofiglio

## This file is part of comet.

    ## comet is free software: you can redistribute it and/or modify
    ## it under the terms of the GNU General Public License as published by
    ## the Free Software Foundation, either version 3 of the License, or
    ## (at your option) any later version.

    ## comet is distributed in the hope that it will be useful,
    ## but WITHOUT ANY WARRANTY; without even the implied warranty of
    ## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    ## GNU General Public License for more details.

    ## You should have received a copy of the GNU General Public License
    ## along with gcipdr.  If not, see <https://www.gnu.org/licenses/>.


                                        # DESCRIPTION:
# 1) calccinc computes CIF, log ratio of CIFs, and their standard deviation (using delta method derivation).
# 2) metacinc is a routine to meta-analyze these estimates.
# 3) prep.metacr prepares data in the right format.


#' @title CIF computation
#'
#' @description calccinc computes CIF, log ratio of CIFs, and their standard deviation (using delta method derivation).
#'
#' @seealso metacinc

#########################################
###### MAIN COMPUTATIONAL METHODS #######
#########################################

calccinc <- function(events,cause, # competing events in ONE study
                      N, ptime, # riskset in ONE study for both events
                     timepoint # any time sequence
                     ){ ## result: a list

    k <- length(events) # number of causes

  if(length(cause)!=k)
      stop("'cause' and 'events' length differs")
    
  if(length(ptime)!=k)
      stop("'ptime' and 'events' length differs")
    
    t <- length(timepoint) # time sequence

   
  ## first partial derivative in respect to lambda_1
  
  d.lambda1 <- function(lambda1, lambda, timepoint){
      
          res <- (1/lambda1 - 1/lambda + exp(-lambda*timepoint) /
            (1-exp(-lambda*timepoint)) * timepoint)
    res
  }
  ## second partial derivative in respect to lambda_2
  
  d.lambda2 <- function(lambda, timepoint){
    res <- (-1/lambda + exp(-lambda*timepoint) /
            (1-exp(-lambda*timepoint)) * timepoint)
    res
  }
  ##
  lambda1 <- as.numeric( events/ptime ) # cause-spec. hazards, a vector
  lambda  <- as.numeric( sum(lambda1) ) # all-cause hazard
      logvarlam <- as.numeric( 1/(events) ) # ptime cancels out 
  timepoint <- as.numeric(timepoint)

 
    
    ## Cause specific cumulative incidence functions


  cif <- do.call("rbind",
 
     lapply( 1:t , 

           function(j) sapply( 1:k ,

          function(i)
           ( 1-exp(-lambda*timepoint[j]) )*(lambda1[i]/lambda)                   )
                   )

            ) 

            

    ## Asymptotic variance ## 


     varlogcif <- do.call("rbind",
 
     lapply( 1:t , 

           function(j) sapply( 1:k ,

          function(i)
  d.lambda1(lambda1[i], lambda, timepoint[j])^2 * lambda1[i]/ptime[i] +
        d.lambda2(lambda, timepoint[j])^2 * sum(lambda1[-i])/ptime[i]

                              )
                   )

            ) 

    
       

  logcif <- as.vector( log(cif)  ) 
  loglam <- as.vector(
      matrix(rep(log(lambda1), rep(t,k)), ncol=k))

  varlog <- as.vector( varlogcif ) 

    varloglam <- as.vector(
        matrix(rep(logvarlam, rep(t,k)), ncol=k))

  

    RES <- data.frame( LOGRISK= c(logcif, loglam ),
                      VARLR= c(varlog, varloglam),
                      N= rep(rep(N, rep(t,k )),2)  )

    RES$TYPE <- rep( c("CIF","HAZ"), rep(t*k, 2) )
    RES$CAUSE <- rep(rep( cause, rep(t, k)),2)

    
  RES  
  
   } # end calccinc




ccr <- function( event, ptime, cause, N,
                     time, group, study, predt=NULL, #data=NULL,
                pool=FALSE, plot=FALSE, binn=5 ){




  
   SNM <- as.character(unique(study)) #  studies  
   S <- length(SNM) # number of studies
   KNM <- as.character(unique(cause))  # competing events  
   K <- length(KNM) # number of competing events  
   GNM <- as.character(unique(group)) # groups
   G <- length(GNM) # number of groups
       JT <- unique(time) #  timepoints    
      
    JLT <- sort(JT) # sorted time-landmarks 
     if(!is.null(predt)) # with prediction time
       JLT <- sort(c(JT,predt))

   
      if(plot){  # finer segmentation for plotting (with embedded sorted time-landmarks)
          if(is.null(predt))
       JLT <- sort( c(
           seq(0.01,  max(JLT)-0.0569,
               length= 50-length(JT) ),JLT) )
           else
               JLT <- seq( 0.01, max(JLT), length=50)
      }
   
    JLT <- matrix(rep(JLT, S), ncol=S)

   if(!pool){                  # study-specific times
     JLT <- matrix(numeric(S*binn), ncol=S) 
  for(s in 1:S)   
  JLT[,s] <- seq(0.01, JT[s], length.out=binn)    
   }

    colnames(JLT) <- SNM    

   J <- dim(JLT)[1]   # number of timepoints     
    


   DATA <- as.data.frame(
       do.call("rbind",
             mapply(

      function(s)
       do.call("rbind",
               mapply(
     function(g)

           calccinc(event[group==g & study==s],
                                 cause=unique(cause),
                                 N[group==g & study==s],
                 ptime[group==g & study==s], JLT[, s])
           
   , g= GNM ,SIMPLIFY=F   )    )

   , s=SNM  , SIMPLIFY=F

       )    )  )


   
   ######### data creation #########

   

 colnames(DATA) <-
     names(calccinc(1,1,1,1,1))  # colnames extraction
   
    DATA$STUDY <- rep(SNM, rep( G*K*J*2, S ) )
    DATA$GROUP <- rep( rep( GNM, rep( K*J*2, G  ) ) , S)
    DATA$TIME <- as.vector(apply( JLT, 2, function(x)
                         rep( x, K*2*G)  )  )


   DATA
   
    

     } # end ccr


# metacinc DESCRIPTION

# predt: prediction time equal for all studies. If not NULL, overrides argument 'time'.

# pool: If FALSE study specific analyses are reported. If TRUE, a meta-analysis is performed.

# landmark: If 'pool' is TRUE, performs a pooling accordingly to the observed study's times.

# plot: Used with 'pool'=TRUE, if TRUE, it finely segments the time line. The resulting data frame can be fed to plot methods. 

# binn= Used with 'pool'=FALSE, it finely segments the time line.


#' @title Meta-analysis of CIFs (ratios)
#'
#' @description 
#'
#' @details the function has print and plot methods.
#'
#'
#' @seealso calccinc



metacinc <- function( event, N, ptime, cause, 
                     time, group, study, control,
                     predt=NULL, data=NULL, subset=NULL,
                     pool= FALSE, landmark=FALSE,
                     plot=FALSE, binn=5){ # binning time grid
                     

#################################################
######  DATA CATCHING ###########################
################################################# 
    
if (is.null(data))  sys.frame(sys.parent())

 cl <- match.call()  # call for update

 mf <- match.call()  
 mf[[1]] <- as.name("data.frame")
  mf <- eval(mf, data)   # match arguments with data if not NULL


  if(!is.null(mf$subset)){   # subsetting data if subset not NULL
     if(!is.logical(mf$subset))
         stop("Argument 'subset' must be logical.")

         if( (sum(mf$subset) > length(mf$event)) ||
        (length(mf$subset) > length(mf$event)  ) )
      stop("Length of subset is larger than number of trials.")
   else  
         mf <- mf[mf$subset,]  # subsetted data
   }


            event <- mf$event; N <- mf$N; ptime <- mf$ptime
            cause <- mf$cause; time <- mf$time
            group <- mf$group; study <- mf$study

         
                          
####################################################
############## PRELIMINARY WARNINGS ################
####################################################



 if(
                !all(
    aggregate(cause~study+group, FUN=length)[,3]==
                    length( unique(cause) )
  )
  )
     stop("Number of competing events must be equal in all studies.")

            if(
                !all(
    aggregate(group~study+cause, FUN=length)[,3]==
                    length( unique(group))
  )
  )
                stop("Number of groups must be equal in all studies.")


  if(length(grep(control , unique(group) ))==0)
      stop(" Argument 'control' matches with no level from 'group'.")



############################################  
###### ANALYSIS START HERE #################
############################################

  

      

      ###############  KEY STEP: ccr  #########################
      #########################################################
      


  DATA <- ccr( event, ptime, cause, N,
                     time, group, study, predt,
                pool=pool, plot=plot, binn=binn)

  RES <- list(DATA=DATA, binn=binn, call=cl); class(RES) <- "CIFs"
     ##################### pooling option #########################

 if(pool){
     
  GNM <- unique(DATA$GROUP)
  SNM <- unique(DATA$STUDY)
  S <- length(SNM)
  KNM <- unique(DATA$CAUSE)
  K <- length(KNM)
  TNM <- unique(DATA$TIME)
  J <- length(TNM)
  RTP <- c("RR","HR")
  
    ######### CONTRATS: data partition ###########
    
   DC <- subset(DATA, GROUP==control)
   DE <- subset(DATA, GROUP!=control)


  
   contr <- c()
  
 logr <- selr <- np <- c()

  #########  CONTRATS: computation #################
  
  for( g in GNM[GNM!=control]){
   contr <- c(contr,
      paste( g ,"/",control, sep="" )
              )
    logr <- c(logr,
      DE[ DE$GROUP==g, "LOGRISK" ]-DC[ , "LOGRISK"])
    selr <- c(selr,
       sqrt( DE[ DE$GROUP==g, "VARLR" ]+DC[ , "VARLR"] ))
     np <- c(np, DE[ DE$GROUP==g, "N" ]+DC[ , "N"]) # study sample size
  
   
   }# end for 

   
      ######## LOGRATIOS DATASET #############

  
    MD <- data.frame(LOGRATIO= logr, SELR=selr, N=np)
  
  C <- length(contr)


 MD$STUDY <- rep( DC$STUDY, C )
 MD$CAUSE <- rep( DC$CAUSE, C)
 MD$TIME <- rep(DC$TIME, C)
 MD$TYPE <- rep(ifelse(DC$TYPE=="CIF", RTP[1], RTP[2]), C )
 MD$CONTR <- rep(contr , rep(K*J*S*2, C) )

  
      ##########  META - ANALYSIS ###########

  
   ## main function
  
  nested.meta <- function(k,c,i,t){


 md <- subset( MD,
         CAUSE==k &  CONTR==c & TYPE==i & TIME==t ,
               select= 1:4 );  md$fup <- unique(time)
  
  ############ run metagen() ########################

     ma <- metagen( LOGRATIO, SELR, STUDY,  data=md,  
                n.e= md$N )                         
  ############### landmark option ###################
  
  if(landmark & is.null(predt) )
      ma <- metagen( LOGRATIO, SELR, STUDY,  data=md,  
                n.e= md$N, subset= fup>=t )                         
  ###################################################

  #  browser()
 
      ai <- data.frame( k=ma$k, Q=ma$Q, df.Q=ma$df.Q,
                        tau=ma$tau, H.TE=ma$H,
                       H.lower=ma$lower.H,
                       H.upper=ma$upper.H,
                       I2.TE=ma$I2, I2.lower=ma$lower.I2,
                       I2.upper=ma$upper.I2 ) # additional info 


          ## cbind( as.data.frame( ma[c(22:25)] ),
          ##    as.data.frame(summary(ma)[9:10]) ) # additional info 

       res <-       data.frame(  # meta-dataset in loop
  LOGPOOL= c( ma$TE.fixed, ma$TE.random),  # logte
  SELP = c( ma$seTE.fixed, ma$seTE.random),  # logse
  LOGLOW= c( ma$lower.fixed, ma$lower.random),  # logcil
  LOGUP= c( ma$upper.fixed, ma$upper.random),  # logciu
      rbind(ai,ai),                # more infos
       N =  rep(sum(ma$n.e) ,2),   # total participants
       TIME =  rep(t,2), # time
       TYPE =  rep(i,2), # type
      CONTR =  rep(c,2), # contrast
   CAUSE  =  rep(k,2) # cause
                                   )          


   ##         cbind(
   ## unlist( ma[c(6,12)], use=FALSE),  # logte
   ##     unlist( ma[c(7,13)], use=FALSE),  # logse
   ##  unlist( ma[c(8,14)], use=FALSE),  # logcil
   ##  unlist( ma[c(9,15)], use=FALSE),  # logciu
   ##    rbind(ai,ai),                # more infos
   ##       rep(sum(unlist(ma[48],use=FALSE)),2),   # total participants
   ##        rep(t,2), # time
   ##        rep(i,2), # type
   ##        rep(c,2), # contrast
   ##         rep(k,2) # cause
   ##                                 )   

                    res  
             
                         }

    ## results

   MA <- as.data.frame(
       do.call("rbind",

                 mapply(

    function(k)
       mapply(

    function(c)
       mapply(

    function(i)
       mapply(

    function(t)
           
       nested.meta(k,c,i,t)
           
    ,t=TNM     ,SIMPLIFY=F
           ) 
    ,i=RTP  
           )
           
    ,c=contr  
             )
    ,k=KNM  
            )

        )  )
  
                 
# fk <- metagen(1:2,1:2)  # fake for names extraction


  
  ## colnames(MA)[c(1:4,15:19)] <-
  ##     c( "LOGPOOL", "SELP", "LOGLOW",
  ##      "LOGUP", 
  ##   "N", "TIME", "TYPE", "CONTR",
  ##       "CAUSE"
  ##       )
  
    MA$MODEL <- rep(c("FIXED","RANDOM"),  J*2*C*K )
  
  
  RES <- list( MA=MA, MD=MD, DATA=DATA,fup=unique(time),
              plot=plot, call=cl)
  class(RES) <- "RR"
  
   }# end if pool

  
      class(RES)[2] <- "metacif"

  RES
  
    }# end metacinc


  

###################
## print methods ##
###################

#Usage: to print results just type the name of the metacinc object.
# Es. (not run)
# ma1 <- metacinc('your argument here'); ma1 

 print.metacif <- function(object, pred=FALSE, # prediction scale
                           hazr=TRUE, uncombined=TRUE
                           ){

  match.call()

     if (!inherits(object, "metacif")) 
        stop("Argument 'object' must be an object of class \"metacif\"")

       DATA <- object$DATA
    ci <-with( DATA, # log_CI-95%

         ci( LOGRISK, sqrt(VARLR) )
               ) 

             KNM <- unique(DATA$CAUSE) 
             SNM <- unique(DATA$STUDY)
             GNM <- unique(DATA$GROUP)
             RTP <- unique(DATA$TYPE)
             S <- length(SNM); K <- length(KNM)
             G <- length(GNM)

  
     ## original scale (%) ##
   DATA$RISK <- round(exp(DATA$LOGRISK)*100,2) 
     DATA$LOW <- round(exp(ci$lower)*100,3)
     DATA$UP <- round(exp(ci$upper)*100,3)
  if(pred){
  ## prediction scale ##
  DATA$RISK <- round(with(DATA, exp(LOGRISK)*N),0)
    DATA$LOW <- round(exp(ci$lower)*DATA$N,0)
      DATA$UP <- round(exp(ci$upper)*DATA$N,0)
   }
  
  DATA$PRINT <-with(DATA,
                      paste( RISK, "[", LOW, ";", UP, "]", sep=""))
  
  
         if (inherits(object, "CIFs")){  # methods for CIFs
    if(object$binn>10)
        print("Time binning is too fine (binn>10). Printing methods are disabled. Use names() to inspect object.", quote=FALSE)
    else{

        B <- object$binn

        if(!pred)
        cat("\n #####################################\n",
              "### PROPORTION OF EXPECTED CASES ####\n",
              "#####################################\n")
              else
                  cat("\n #####################################\n",
              "##### NUMBER OF EXPECTED CASES ######\n",
              "#####################################\n")
                  
cat( paste("Number of competing events:", K, sep=""), "\n")
 cat( paste("Number of studies:", S, sep="" ), "\n" )
      cat( paste("Number of groups:", G, sep="" ), "\n" )
cat(paste("Time binning:", B, sep=""),"\n")

      cat( " ++++++++++++++++++++++++++++++++++++++ \n")    

          cat("\nLegend: \n")

  cat("IR = Incidence rate (ML-estimate).\n")
     cat("ML= Maximum likelihood.\n")
          cat("CI= 95% Confidence intervals*.\n")
  if(!pred)
        cat("CIFs= Cumulative incidence functions.\n")
        else
            cat("Events= Number of expected events.\n")
      cat( paste("Groups="), as.character(GNM), "\n" )

  ## loop ##
     for(i in KNM){

                 cat(" \n #####################\n",
             "  EVENT:", i, " \n",
             "#####################\n
               \n" )
                 
                  for(s in SNM){
                      
      hz <- subset(DATA, TYPE=="HAZ" & CAUSE==i & STUDY==s)
      cf <- subset(DATA, TYPE=="CIF" & CAUSE==i & STUDY==s)
        TNM <- unique(cf$TIME); J <- length(TNM)
             NP <- cf$N
      
      csh <-unique(subset(hz, select="PRINT" ))
        colnames(csh) <- " "
      rownames(csh) <- GNM      


      tb <- matrix(
          as.matrix(
              subset(cf[order(cf$TIME),],
                           select="PRINT")), ncol=J)
  rownames(tb) <- GNM; colnames(tb) <- paste("time ", TNM, sep="") 
  
      
        cat("\n event: ", i , "\n study: ", s ,
            "\n participants: ", paste(GNM,": ",unique(NP),".  ", sep="") ,"\n \n
        IR ML-estimate (%) [CI]: \n")
            print(csh)
        if(!pred)
        cat( paste("\n CIFs (%) [CI] time-profile:"), "\n \n" )
      else
          cat( paste("\n Events [CI] time-trend:"), "\n \n" )
       print(noquote(tb))
     cat(rep("-  -", B*3 ), "\n")

    } # for i
          } # for s


     cat(" \n
    \nDetails on computational methods:\n")
   cat("- Maximum likelihood estimate for the hazard (IR).\n")
        if(!pred)
        cat("- IR-based Cumulative incidence functions (CIFs).\n")
        else
            cat("- IR-based Prediction of expected cases.\n")
             cat("- Variances: *Delta-method-based.\n")


     
         } # matches else binn

         } # matches if CIFs



          if (inherits(object, "RR")){  # methods for RR
   if(object$plot)
       print("metacif() was in plot mode (plot=TRUE). Printing methods are disabled. Use names() to inspect results. ", quote=FALSE)
              else{
  warning("Cochrane Q test for heterogeneity is performed at different time points. There may be multiple testing problems.", call.=FALSE)

    MA <- object$MA; MD <- object$MD
  
ci <- with(MD, ci(LOGRATIO, SELR))

       ### original scale and rounding ###
  MD$RATIO <- round(exp(MD$LOGRATIO),2)
  MD$LOW <- round(exp(ci$lower),3)
  MD$UP <- round(exp(ci$upper),3)
  
  MA$PEFF <- round(exp(MA$LOGPOOL),2)
    MA$LOW <- round(exp(MA$LOGLOW),3)
    MA$UP <- round(exp(MA$LOGUP),3)
  MA[, 8:14] <- round(MA[, 8:14], 2) 
  
  CTR <- as.character(unique(MA$CONTR))
 RTPM <- as.character(unique(MA$TYPE))
 MDM <- as.character(unique(MA$MODEL))
 TNM <- unique(MD$TIME)
 J <- length(TNM)
  tol <- 0.05
  qst <- with( MA,
      which( abs((PEFF-LOW)-(UP-PEFF))>tol) )
if(length(qst)>0)
  warning(
     paste( " 95% confidence bands for pooled effects are not symmetric. Convergence to normality has failed. Number of events, sample size, or person-time-at-risk are too small. Log standard errors are too big (MA$SELP>0.05)." ) , call.=FALSE)


        p <- with(MA, round( pchisq(Q, df.Q, lower=FALSE), 5) )

  p <- noquote(ifelse(p<1e-04, "< 0.0001", p ))


 ### table editing ###

MA$RR.print <- with(MA, paste(PEFF, "[", LOW, ";", UP, "]", sep="") )
MA$tau2.print <- with(MA, paste( tau , " (", p, ")" ,sep="") )
MA$H.print <- with(MA, paste( H.TE, "[",H.lower, ";", H.upper, "]",sep=""))
MA$I2.print <- with(MA,
       paste( I2.TE*100, "[", I2.lower*100, ";", I2.upper*100, "]",sep="") )
  
MD$PRINT <- with(MD,
           paste(RATIO, "[", LOW, ";", UP, "]", sep="") )
  
  
   ### PRINT RESULTS

 cat(" ############################################\n",
          "####  COMPETING RISKS META-ANALYSIS  #######\n",
          "############################################\n")

  cat("\nLegend: \n")

      cat("CIFs = Cumulative incidence functions.\n")
  cat("RR = Risk Ratio (CIFs Ratio).\n")
    cat("Hazard = Incidence rate (IR).\n")
  cat("HR = Hazard ratio (IR ratio).\n")
  cat("RR* = Pooled RR.\n")
    cat("HR* = Pooled HR.\n")
  cat("CI = 95% Confidence Intervals*.\n")
  cat("k = Number of combined studies.\n")
  cat("N = Number of participants.\n")
  cat("tau^2 = Estimate of between study variance.\n")
  cat("pval = P-value of test Q for heterogeneity.\n")
  cat("H = Q/df.Q ; H=1 means no variation.\n")
  cat("I^2 = % of variation across studies (Higgins and Thompson, 2002).\n")
      cat( paste("Groups="), as.character(GNM), "\n" )
  cat( paste("Comparisons="), as.character(CTR), "\n" )


  
    for(i in KNM){

  dt <- subset(DATA, CAUSE==i)
  md <- subset(MD, CAUSE==i)
  ma <- subset(MA, CAUSE==i)
        
        cat("\n #####################\n",
             "  EVENT:", i, " \n",
            "#####################\n
               \n" )   

        cat("\n
    Hazard (%) [CI] and HR [CI]:\n \n")

  csh <- matrix(
      as.matrix(unique(
          subset(dt,
                 TYPE=="HAZ", select="PRINT"))),ncol=S)
   colnames(csh) <- SNM; rownames(csh) <- GNM

  hr <- matrix(
      as.matrix(unique(
          subset(md,
                 TYPE=="HR", select="PRINT"))),ncol=S)
   colnames(hr) <- SNM; rownames(hr) <- CTR

  print(noquote(rbind(csh,hr)))
  
    if(uncombined){
        
    cat("\n
         Uncombined RR [CI]:\n ")

    for(c in CTR){
        cat("\n Event: ", i," \n Contrast: ", c,"\n" )
            tb <-
              matrix(  as.matrix(
        subset(md[order(md$TIME),],
               CONTR==c & TYPE=="RR" 
               , select="PRINT") ),
                     nrow=S )
           rownames(tb) <- SNM
    colnames(tb) <- paste("time ", TNM, sep="")
  print(noquote(tb))
    }                    }# if uncombo


  
   malabs <- 
        c("RR* [CI]", "Time"," k", "N",
           "tau^2 (pval)"," H [CI]"," I^2 (%) [CI]")
  
            cols <- c(24, 16, 5, 15, 25:27)
  if(!hazr)
      RTPM <- RTPM[1]
  for(c in CTR){
                            cat("\n \n Combined results \n",
                      "######################\n",
                       "Event: ", i,"\n Contrast : ", c,"\n",
                      "######################\n" )
  for(r in RTPM)    
      for(m in MDM){
  if(r=="HR")
          malabs[1] <- "HR* [CI]"
          if(m=="RANDOM"){
              cols <- cols[1:4]
          malabs <- malabs[1:4] }

      cat(" \n Event: ", i," \n",
      "Pooled", r ,"(",tolower(m), "effects model):\n
         \n" )    
  
  rf <- subset(ma,
               CONTR== c & TYPE== r & MODEL== m,
               select=cols )
   colnames(rf) <- malabs

  
      print(rf)
               } # end for MDM

                        }# end for CTR
      } # end for KNM


   # last Details

     cat(" \n
    \nDetails on meta-analytical method:\n")
  cat("- Maximum likelihood estimate for the hazard (IR).\n")
        cat("- IR-based Cumulative incidence function (CIFs).\n")
   cat("- Time-dependent competing risks meta-analysis.\n")
cat("- Inverse variance method.\n")
cat("- DerSimonian-Laird estimator for tau^2.\n")
   cat("- Estimator for RR: CIFs ratio.\n")
   cat("- Variances: *Delta-method-based.\n
          \n")



       } #matches else plot
          } # matches if RR
     

  } # end print.metacif


####################
# plotting methods #
####################

# dependency: ggplot2


plot.metacif <- function(object, contr=1,
                          ev_lab=NULL, arm_lab=NULL,
                          x_lab=NULL, y_lab=NULL,title=NULL,
                          pred=FALSE, hazr=TRUE, # prediction scale
                          fixed=TRUE,text=FALSE,
                           logtrans=FALSE, ci=TRUE,# conv=NULL,
                          x_extra= 0.0, y_extra= 0.0){ # xlim, ylim 
                          
     match.call()
     
         if (!inherits(object, "metacif")) 
        stop("Argument 'object' must be an object of class \"metacif\"")


     if (inherits(object, "CIFs")){ # if CIFs


  DATA <- subset(object$DATA, TYPE=="CIF")

      ci <-with( DATA, # log_CI-95%

         ci( LOGRISK, sqrt(VARLR) )
               ) 

     ## original scale (%) ##
   DATA$RISK <- exp(DATA$LOGRISK) 
     DATA$LOW <- exp(ci$lower)
     DATA$UP <- exp(ci$upper)
         yl <- ylab("Cumulative Event Probability")
  tl <- "Risk Profile: proportion of expected cases"
   if(pred){
  ## prediction scale ##
  DATA$RISK <- round(with(DATA, exp(LOGRISK)*N),0)
    DATA$LOW <- round(exp(ci$lower)*DATA$N,0)
      DATA$UP <- round(exp(ci$upper)*DATA$N,0)
           yl <- ylab("Events")
  tl <- "Risk Profile: number of expected cases"

    }

  #labels
        if(!is.null(ev_lab)){
       DATA$CAUSE <- as.factor(DATA$CAUSE)
       levels(DATA$CAUSE) <- ev_lab
   }
  
       DATA$GROUP <- as.factor(DATA$GROUP)
        if(!is.null(arm_lab)){
       levels(DATA$GROUP) <- arm_lab
   }
  
       if(!is.null(title))  
           tl <- title
       # layers
  mp <- aes( TIME, RISK, linetype= GROUP ) 
  fg <- facet_grid( CAUSE ~ STUDY )
  
  geom <- geom_line( aes(linetype= GROUP ), 
                    size=1)  

           
  lg <- scale_linetype_discrete(name="Arm")
  pl <- ggplot(DATA, mp )
       
  bckgr <- theme_bw()
  tl <- ggtitle(tl)
       th <- theme(plot.title=element_text(face="bold"))

                        mint <- min(pl$data$TIME) # time-range
 maxt <- max(pl$data$TIME)

       xlim <- NULL
       if(x_extra!=0.0)
 xlim <- xlim( mint , maxt + x_extra)

       tmtr <- NULL
       if( logtrans )
     tmtr <- scale_x_continuous(trans="log10", labels= identity)

  
       xl <- xlab("Time")   
           if(!is.null(x_lab))
               xl <- xlab(x_lab)
           if(!is.null(y_lab))
               yl <- ylab(y_lab)
       
     #  print(
           
 out <-  pl + xlim + fg + geom + bckgr + lg + yl + xl + tl +th + tmtr
           
 
           
     #      )

   }  # end if CIFs



     if (inherits(object, "RR")){  # if RR

         fup <- object$fup
         call <- object$call
         DATA <- object$MA
         CTR <- unique(DATA$CONTR)
         DATA$PEFF <- exp(DATA$LOGPOOL)
    DATA$LOW <- exp(DATA$LOGLOW)
    DATA$UP <- exp(DATA$LOGUP)
    
         DFR <- subset(DATA, CONTR==CTR[contr]) # select contrast 
      if(!fixed)                                # select model
           DFR <- subset(DFR, MODEL=="RANDOM" )
          else
         DFR <- subset(DFR, MODEL=="FIXED")    

         if(!hazr)
             DFR <- subset(DFR, TYPE=="RR")

                 ## relabelling ##
         DFR$TYPE <-with(DFR,
    ifelse(TYPE=="RR", "Pooled CIF Ratio",  # CUMULATIVE RISK
                      "Pooled HAZARD Ratio"   ) # INCIDENCE RATE
         )
DFR$MODEL <- with(DFR,
    ifelse(MODEL=="FIXED", "Fixed-effect","Random-effect" )
         )
         

     #labels
        if(!is.null(ev_lab)){
       DFR$CAUSE <- as.factor(DFR$CAUSE)
       levels(DFR$CAUSE) <- ev_lab
   }

              
         
      tl <- paste("Pooled effect: ",CTR[contr] ,".\n",  # automatic title
                  as.character(unique(DFR$MODEL)),
                  " model\n", sep="") 

     if(!is.null(title))
          tl <- title
  # ++++++++++++++   PLOT LAYERS +++++++++++++++++++++++  #                    
  tl <- ggtitle(tl)
  th <- theme(plot.title=element_text( face="bold"))
         
  
    mp <- aes(TIME, PEFF )
              

       
                  pl <- ggplot(DFR , mp)
      panel <- facet_grid(CAUSE~TYPE)

         geom1 <- geom_line( size=1)
         

  bckgr <- theme_bw()
        
 one <- geom_hline( yintercept=1, linetype=2 , colour="black")


         
                        mint <- min(pl$data$TIME)
 maxt <- max(pl$data$TIME)
 miny <- min(pl$data$LOW)
 maxy <- max(pl$data$UP)

       xlim <- ylim <- NULL

        if(y_extra!=0.0)

            ylim <- ylim( miny- y_extra,  maxy+ y_extra)

        if(x_extra!=0.0)
 xlim <- xlim( mint , maxt + x_extra)



    txtu <- txtl <- geom2 <- NULL; sz <- 0.0
                      
                       if(ci){
         geom2 <- geom_ribbon(aes( ymin=LOW, ymax=UP ),
                              alpha=0.2)
         sz <- 3.2

   txtl <- geom_text(data=subset(pl$data, TIME==max(TIME) ),
                   aes( TIME, LOW,  label= round( LOW, 3 ) ),
                    fontface= 1 , size= sz ,  
                    vjust= 1.4, hjust= 0.7
                     )

         txtu <- geom_text(data=subset(pl$data, TIME==max(TIME) ),
                   aes( TIME, UP,  label= round( UP, 3 )),
                    fontface= 1 , size= sz ,  
                    vjust= 1.4, hjust= 0.7
                     )


         
     }
         
 

         txtc <- geom_text(data=subset(pl$data, TIME==max(TIME) ),
                   aes( TIME, PEFF,  label=round( PEFF, 3 )),
                    fontface= 1 , size=3.2 ,  
                    vjust= 1.4, hjust= 0.7
                     )
         
          pltck <- lgnd2 <- lgnd <- txt2 <- txt1 <-  geom3 <- NULL



        if(text & !is.null(call$landmark) )
            if(call$landmark){


                
   ticks <-with(pl$data,
                data.frame(x=sort(fup),y=min(LOW)-y_extra,
                       k=unique(k),
                       N=unique(N) ) )
   S <- dim(ticks)[1]

   
    colfunc <- colorRampPalette(c("red", "green"))
   
   geom3 <- geom_point( aes(x,y, colour=factor(k) ), data=ticks
                      , shape=17  #, colour= colfunc(  S )
                       ) 

   txt1 <- geom_text(data=ticks, aes( x, y,  label= k),
                    fontface= 3 , size=3 ,  
                    vjust= -0.5, hjust= 0
                     )

  txt2 <- geom_text(data=ticks, aes( min(x), y,  label= "N. Studies:"),
                    fontface= 1 , size=3.6 ,  #changed font 6 to 3
                    vjust= -0.55, hjust= 1.3  # ch -0.65, 1.1
                     )
   mt <- round(median(c(0, min(ticks$x))),1)
   
   pltck <- scale_x_continuous( breaks=c(0,mt, ticks$x),
             labels=c(0, mt, round(ticks$x[1],1),
                 rep( " ", S-2), ticks$x[S]))
   

  lgnd2 <- scale_colour_manual(name="N. Participants",
   values= rev(colfunc( S)),
                               breaks=pl$data$k,
                        labels=as.character(pl$data$N))

   
         }

    
         
        lgtr <- NULL
        if(logtrans)

            lgtr <- scale_y_continuous(trans="log"
                       # something changed in ggplot2. add new breaks arg.
                      ## , breaks=trans_breaks('log',
                        ##      function(x) round(exp(x),2) )
                                       
                                         )

         yl <- ylab("Pooled Ratio")
        xl <- xlab("Time") 

        if(!is.null(x_lab))
               xl <- xlab(x_lab)
           if(!is.null(y_lab))
               yl <- ylab(y_lab)



   #     print(

  out <- pl + xlim + ylim + panel + one + geom1 + geom2 + bckgr + xl + yl + tl + th + geom3  + lgnd + txtl + txtc + txtu + txt1 + txt2 + lgnd2 + lgtr + pltck
   #        ) 

     

   }  # end if RR
   

  out
     
   } # end plot.metacif



##################################
#### data preparation methods #####
##################################


prep.metacr <- function( event, ptime, N, # list of data.frames**
                 studlab, fup, data=NULL){
                 # **list name : group
                 # **matrix column: cause

 if (is.null(data))  sys.frame(sys.parent())

  mf <- match.call()

  mf[[1]] <- as.name("list")

 mf <- eval(mf, data)

 event <- mf$event; ptime <- mf$ptime
 studlab <- mf$studlab; fup <- mf$fup
 N <- mf$N


 obj <- list(event, ptime, N)

if(is.null(names(event)) || is.null(names(ptime)) ||
   is.null(names(N)))
 stop("'event', 'ptime', 'N': list elements must have a name.")

 if(!is.list(event) || !is.list(ptime) || !is.list(N))
     stop("'event', 'ptime', 'N': must be a 'list'.")
 if(!all(
  unlist(lapply(obj, function(x) lapply(x,
         function(x) is.data.frame(x))))) )
 stop("'event', 'ptime', 'N': list element must be a data.frame.")
    
 
 GNM <- names(event)
 G <- length(GNM)
 KNM <- colnames(event[[1]])
 K <- length(KNM)
    
    if(!all(
  unlist(lapply(obj, 
         function(x) names(x)==GNM))) )
 stop(paste("'event', 'ptime', 'N': character vector",
            GNM, "differs among list elements. Groups must be the same across arguments."))
       
 if(!all(
  unlist(lapply(obj, function(x) lapply(x,
         function(x) colnames(x)==KNM)))) )
 stop(paste("'event', 'ptime', 'N': character vector",
            KNM, "differs among list elements. Causes must be the same across groups."))

    if(!all(
  unlist(lapply(obj, function(x) lapply(x,
         function(x) dim(x)[1]==dim(event[[1]])[1])))) )
    stop("'event', 'ptime', 'N': dim(x)[1] differs among list elements. Number of records must be the same across arguments.")
       if(!all(
  unlist(lapply(obj, function(x) lapply(x,
         function(x) dim(x)[2]==dim(event[[1]])[2])))) )
    stop("'event', 'ptime', 'N': dim(x)[2] differs among list elements. Number of causes must be the same across arguments.")

          S <- length(studlab)

if(S!=dim(event[[1]])[1])
    stop("'studlab': lenght differs from 'event' record.")
          if(S!=length(fup))
    stop("'fup': lenght differs from number of studies.")


   event <- unlist(as.matrix(event), use=FALSE)
      ptime <- unlist(as.matrix(ptime), use=FALSE)
              N <- unlist(as.matrix(N), use=FALSE)
DATA <- data.frame(event=event, ptime=ptime,
                   N=N, study=studlab, time=fup)
          
           DATA$group <- rep( GNM, rep( S*K ,G) )
 DATA$cause <- rep(rep(KNM, rep(S, K)),G)


 DATA
 

  } # end prep.metacr







