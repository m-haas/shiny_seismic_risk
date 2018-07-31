#############################################################################################
# DESERVE: Seismic Risk Calculator - Local Scale
#
# Israeli Towns
#     - Arad
#     - Megiddo
# Jordan Towns
#     - Kerak
#     - Madaba
# Palestinian Towns
#     - Nablus
#
#############################################################################################

library(shiny)
library(shinyIncubator)
library(RPostgreSQL)

shinyServer(function(input, output,session) {

  ###########################################################################
  ###########################################################################
  ##                                                                       ##
  ##                  A)     Calculation Routines                          ##
  ##                                                                       ##
  ###########################################################################
  ###########################################################################  
  
  getInt <- reactive({
  ###########################################################
  # Create grid from user input                             #
  ###########################################################
  # Sets up a grid depending on the user input:
  # input$lon, input$lat, input$radius, input$grdStep
   getGrd <- function(){ 
    
    # Convert inputs from km to deg
    # Earth radius [km]
    R <- 6371
    # Longitude extent [?]
    rx <- input$radius/(R*cos(input$lat*pi/180))*180/pi
    # Latitude extent [?]
    ry <- input$radius/R*180/pi
    # Grid step longitudinal [?]
    dx <- input$grdStep/(R*cos(input$lat*pi/180))*180/pi
    # Grid step latitudinal [?]
    dy <- input$grdStep/R*180/pi
    
    #number of nodes (size)
    nx <- 2*round(rx/dx)+1
    ny <- 2*round(ry/dy)+1
    size <- nx*ny
    
    #Initialize x,y with upper left corner coordinates
    x  <- input$lon - rx
    y  <- input$lat + ry
    
    #create mesh (grd[,1]: longitudes ; grd[,2]: latitudes)
    #first entry
    grd <- c(x,y)
    # loop through whole array
    for (i in 2:size){
      # coordinates to write in grd
      # check if max longitude reached (1 row of array completed)  
      if ((i-1)%%nx == 0){
        x <- input$lon - rx #set longitude back to start value
        y <- y - dy #decrease y by grid step
      }else x <- x + dx #increase x by grid step
      #attach to grd
      grd <- rbind(grd,c(x,y))
      
    }#end for  
        
    return(grd)
  }#end function "getGrd"
  
  ###########################################################
  # Approximate endpoints of line                           #
  # Wells&Coppersmith 1994 FLR/FWR(Surface Rupture)         # 
  ###########################################################
  getEnds <- function(){
    # Estimating surface rupture length and width
    rl <- -2.42 + 0.58 * input$Mw
    rw <- -1.61 + 0.41 * input$Mw
    
    flt <- c(10^rl,10^rw)
    
    lat1 <- input$lat
    lon1 <- input$lon
    len <- flt[1]/2
    theta <- input$strike*pi/180
    for (d in seq(0.05,4,0.05)){
      
      lat2 <- lat1+d*sin(theta)
      lon2 <- lon1+d*cos(theta)
      if (exists("est")) est <- rbind(est,d) else est <- d
      if (exists("dst")){
        tmp <- spDistsN1(as.matrix(cbind(lon1,lat1)),c(lon2,lat2),longlat=TRUE)
        dst <- rbind(dst,tmp)
      }else{
        dst <- spDistsN1(as.matrix(cbind(lon1,lat1)),c(lon2,lat2),longlat=TRUE)
      }
    }
    dsel <- which(dst >= len)[1]
    p1 <- c(lon1+est[dsel]*cos(theta) , lat1+est[dsel]*sin(theta)) 
    p2 <- c(lon1-est[dsel]*cos(theta) , lat1-est[dsel]*sin(theta))
    
    return(rbind(p1,p2))
  }
  
  ###########################################################
  # Calculate Distances                                     #
  ###########################################################
  
  getDist <- function(grd){
    #grd <- getGrd()
    if(input$sourceType=="pt"){
    dst <- spDistsN1(grd,c(input$lon,input$lat),longlat=TRUE)
    }else{
      ends <- getEnds()  
      
      # rupture as vector approximating Earth as sphere
      R <- 6371
      
      dst <- matrix(0,length(grd[,1]),1)
      for (i in 1:length(grd[,1])){
        # function convert to cartesian
        xyz <- function(lon,lat){
          rlon <- lon*pi/180
          rlat <- lat*pi/180
          x <- cos(rlat)*cos(rlon)
          y <- cos(rlat)*sin(rlon)
          z <- sin(rlat)
          return(cbind(x,y,z))
        } 
        # function crossproduct
        xprod <- function(u,v){
          res <- c(u[2]*v[3]-u[3]*v[2],u[3]*v[1]-u[1]*v[3],u[1]*v[2]-u[2]*v[1])
          return(res)
        }
        #function norm of vector
        vnorm <- function(v){
          res <- sqrt(v%*%v)
          return(res)
        }
        #function xyz to spherical
        xyz2sph <- function(v){
          d <- sqrt(v[1]^2+v[2]^2+v[3]^2)
          lon <- atan2(v[2],v[1]);
          lat <- asin(v[3]/d);
          return(cbind(lon,lat))
        }
        #function xy to polar
        xy2rphi <- function(v){
          r <- sqrt(v[1]^2+v[2]^2)
          phi <- atan2(v[2],v[1])
          return(cbind(r,phi))
        }
        
        # convert given lon lat to cartesian
        ll <- rbind(ends[1,],ends[2,],grd[i,])
        cc <- matrix(0,3,3)
        for (j in 1:3){
          cc[j,] <- xyz(ll[j,1],ll[j,2])
        }
        
        # calculation (great circle estimate)
        v0 <- cc[1,] #start line
        v1 <- cc[2,] #end line
        v2 <- cc[3,] #point of interest
        plane     <- xprod(v0,v1)
        circle    <- xprod(v2,plane)
        intersect <- xprod(plane,circle)
        intersect <- intersect/vnorm(intersect)*R
        
        #back to lon lat
        pt <- xyz2sph(intersect)*180/pi
        # get polar coordinates on great circle
        c1 <- xy2rphi(v0)
        c2 <- xy2rphi(v1)
        c3 <- xy2rphi(intersect)
        #check if intersect point within arcsegment and get distance
        if(c3[2] <= max(c1[2],c2[2]) && c3[2] >= min(c1[2],c2[2]) ){
          dst[i] <- spDistsN1(as.matrix(cbind(grd[i,1],grd[i,2])),pt,longlat=TRUE)  
        }else{
          dStart <- spDistsN1(as.matrix(cbind(grd[i,1],grd[i,2])),ends[1,],longlat=TRUE)
          dEnd <- spDistsN1(as.matrix(cbind(grd[i,1],grd[i,2])),ends[2,],longlat=TRUE)   
          dst[i] <- min(dStart,dEnd)
        }        
      } #end for
    }#end else sourceType
    return(dst)
  } #end function getDist
  
  ###########################################################
  # Intensity Prediction Equations                          #
  ###########################################################
  
    grd <- getGrd()
    dst <- getDist(grd)
    h <- input$z
    #case 1: Point & EMCA
    if((input$sourceType == "pt") && (input$IPE == "0")){
      
      xEP <- c(1.0074538,-2.0045088,3.2980663,2.6920855,4.2344195e-04)
      sigmaEP <- 0.69073468
      int <- xEP[1]*input$Mw+xEP[2]*log10(h)+xEP[3]-xEP[4]*0.5*log10((dst/h)^2+1)-xEP[5]*(sqrt(dst^2+h^2)-h)
    
    }else if((input$sourceType == "li") && (input$IPE == "0")){
    #case 2: Line & EMCA
      xEX <- c(0.88335012,-1.7719603,3.5466419,2.5910521,9.8232418e-05)
      sigmaEX <- 0.69808861 
      int <- xEX[1]*input$Mw+xEX[2]*log10(h)+xEX[3]-xEX[4]*0.5*log10((dst/h)^2+1)-xEX[5]*(sqrt(dst^2+h^2)-h)      
      
    }else if((input$sourceType == "pt") && (input$IPE == "0")){
    #case 3: Point & Global (Allen Wald Worden 2012 JOSE)
      c0 <-  2.085
      c1 <-  1.428
      c2 <- -1.402
      c4 <-  0.078
      m1 <- -0.209
      m2 <-  2.042
      
      #hypocentral distance
      R_hyp <- sqrt(dst^2+h^2)
      
      Rm <- m1 + m2 * exp(input$Mw-5)
      if(R_hyp <= 50){
        int <- c0 + c1*input$Mw + c2 * log(sqrt(R_hyp^2 + Rm^2))
      }else{
        int <- c0 + c1*input$Mw + c2 * log(sqrt(R_hyp^2 + Rm^2)) + c4 * log(R_hyp/50)
      }
      
      s1 <- 0.82
      s2 <- 0.37
      s3 <- 22.9
      
      sigma <- s1 + s2/(1+(R_hyp/s3)^2)
      
    }else{
    #case 4: Line & Global (Allen Wald Worden 2012 JOSE)
      c0 <-  3.950
      c1 <-  0.913
      c2 <- -1.107
      c3 <-  0.813   
      
      R_rup <- sqrt(dst^2+h^2)
      Rm <- 1 + c3 * exp(input$Mw-5)
      int <- c0 + c1 * input$Mw + c2 * log(sqrt(R_rup^2+Rm^2))    
    
    }
    #create raster
    n <- sqrt(length(grd[,1]))
    rast <- rasterFromXYZ(cbind(grd,int))
    
    #save in file
    if (OS==0) path2 <- "/Intensities/" else path2 <- "\\Intensities\\"
    writeRaster(rast,filename=paste(path,path2,"intensities.tif",sep=""),
                format="GTiff",overwrite=TRUE)
    return(rast)
  })#end function getInt
  
  ################################################################################################################
  # LOSS COMPUTATION
  ################################################################################################################
  
  # estimate the loss of a set of settlements
  # in term of expected fatalities, based on a
  # shakemap in macroseismic intensity provided 
  # as raster.
  # saves a shapefile with the settlements exposed to 
  # mmi > 5. for each settlements the bounds for the fatalities
  # are provided, in absolute and relative terms.
  # dbin, dbout: list(db_conn,layer,overwrite) to specify input and output layers in the dbs
    
  compute_loss_points<-function(shakemap,nighttime=FALSE,dbin,dbout)
  {
    # Marc: Clean urbanstructures_loss table if it is not specifically queried (type 2)
    sql <- "delete from urbanstructures_loss;" 
    drv <- dbDriver("PostgreSQL")
    con <- dbConnect(drv, host='lhotse21.gfz-potsdam.de', user= "postgres", password="postgres", dbname="losstool")
    post <- dbGetQuery(con, statement = sql)
    dbDisconnect(con)
    dbUnloadDriver(drv)	     
  
    #define vulnerability - draft
    vuln_emca <- c('b','b','b','d','d','e','e','c','d','a','c','d','e')
    mmi<-c('I','II','III','IV','V','VI','VII','VIII','IX','X')
    occupancy <- c(28, 8, 140, 60,200,200,200,200,140,5,28,5,5)
    nighttime_coeff <- 0.5
    daytime_coeff <- 0.3
    nclass<-13
    # off: offset, depends on the number of columns of input data frame
    off<-4
    
    # daytime or nighttime?
    coeff<-ifelse(nighttime,nighttime_coeff,daytime_coeff)
    
    #load exposure
    # load coords of settlements as spatial points	
    data<-readOGR(dsn=dbin$dbconn,layer=dbin$layer,verbose=FALSE)
    if (is.na(proj4string(data)))
    {
      print('No projection on input file, assigning WGS84')
      proj4string(data)<-"+proj=longlat +datum=WGS84 +no_defs"
    }
    # filter out places with no population information	
    places<-subset(data,data$population>0 & data$population < 1e6)  
    
    # evaluate intensity for each settlement using the shakemap
    intens<-round(extract(shakemap,places@coords,method='bilinear'))
   
     
    
    #ind<-which(intens>5)
    #plot(shakemap)
    #points(places[ind,],col='red')
    
    #### Michael: check added to prevent unneeded error messages
    #check if there is any loss at all --> avoid NA
    check <- intens
    check[which(is.na(intens))] <-0
    
    if(max(check)==0||max(check)<=5){
      # Return safely from calculation if no loss
      msg <- "With this configuration the model expects no loss!"
      print(msg)
      # Marc: clean places_loss table in database
      sql <- "delete from places_loss;" 
      drv <- dbDriver("PostgreSQL")
      con <- dbConnect(drv, host='lhotse21.gfz-potsdam.de', user= "postgres", password="postgres", dbname="losstool")
      post <- dbGetQuery(con, statement = sql)
      dbDisconnect(con)
      dbUnloadDriver(drv)	     
      return(list(totloss_min=0,totloss_med=0,totloss_max=0,msg))
      }else{
      #keep only intensities >5
      ind <- which(intens >5)
      }
      
    loss<-matrix(0,ncol=6,nrow=length(ind))
    tot_loss_min<-0
    tot_loss_med<-0
    tot_loss_max<-0
    i<-1
    for (place in ind) # loop on settlements
    {
      n_fat_min<-0
      n_fat_med<-0
      n_fat_max<-0
      for (type in 1:nclass) # loop on building types
      {	
        #get dpm for most likely vulnerability
        dpm<-ems98_DPM(vuln_emca[type],3)
        
        #compute prob of collapse
        pd4<-dpm[mmi[intens[place]],'d4']
        pd5<-dpm[mmi[intens[place]],'d5']
        
        # expected min (5%) num  
        b_min<-(places@data)[place,off+type]	
        # expected median 
        b_med<-(places@data)[place,off+nclass+type]	
        # expected max (95%) num  
        b_max<-(places@data)[place,off+2*nclass+type]	
        
        # number of expected collapses
        # uses damage grade 5 plus 25% of damage grade 4 (see Spence 2002)
        pd<-(pd5+0.25*pd4)
        n_coll_min<-pd*b_min
        n_coll_med<-pd*b_med
        n_coll_max<-pd*b_max
        
        # add on the expected fatalities
        n_fat_min<-n_fat_min+(occupancy[type]*coeff*n_coll_min)
        n_fat_med<-n_fat_med+(occupancy[type]*coeff*n_coll_med)
        n_fat_max<-n_fat_max+(occupancy[type]*coeff*n_coll_max)
      }
      
      # absolute loss
      if( !(is.na(n_fat_min)) && !(is.na(n_fat_min)) && !(is.na(n_fat_min))){
      loss[i,1]<-round(n_fat_min)
      loss[i,2]<-round(n_fat_med)
      loss[i,3]<-round(n_fat_max)
      
      # total loss
      tot_loss_min<-tot_loss_min+loss[i,1]
      tot_loss_med<-tot_loss_med+loss[i,2]
      tot_loss_max<-tot_loss_max+loss[i,3]
      
      # relative loss: fatalities as percentage of population
      pop<-(places@data)$population[place]
      loss[i,4]<-n_fat_min/pop
      loss[i,5]<-n_fat_med/pop
      loss[i,6]<-n_fat_max/pop
      }
      i<-i+1
      
    }
    
    # write loss estimates to the db
    df2<-data.frame(osm_id=places@data$osm_id[ind],name=places@data$name[ind],population=places@data$population[ind],loss=loss)
    df3<-SpatialPointsDataFrame(places[ind,],df2)
    
    
    #db_conn<-"PG:dbname=test_db host='lhotse21.gfz-potsdam.de' user='postgres' password='postgres'"
    writeOGR(df3,dbout$dbconn,dbout$layer, "PostgreSQL",overwrite_layer=dbout$overwrite)
    msg <- 'Loss for places estimated successfully'
    #sum of loss_min and loss_max values
    msg2 <- sprintf("%s%i%s%i",'Number of total fatalities estimated between: ',sum(loss[,1]),' and ',sum(loss[,3]))
    print(msg)
    print(msg2)
    return (list(totloss_min=tot_loss_min,totloss_med=tot_loss_med,totloss_max=tot_loss_max,msg))
  }
  
  #compute detailed loss for surveyed towns
  # in term of expected fatalities, based on a
  # shakemap in macroseismic intensity provided 
  # as raster.
  # saves a shapefile with the settlements exposed to 
  # mmi > 5. for each settlements the bounds for the fatalities
  # are provided, in absolute and relative terms.
  # dbin, dbout: list(db_conn,layer,overwrite) to specify input and output layers in the dbs
  
  compute_loss_geocells<-function(shakemap,nighttime=FALSE,dbin,dbout)
  {
    # Marc: Clean places_loss table if it is not specifically queried (type 1)
    sql <- "delete from places_loss;" 
    drv <- dbDriver("PostgreSQL")
    con <- dbConnect(drv, host='lhotse21.gfz-potsdam.de', user= "postgres", password="postgres", dbname="losstool")
    post <- dbGetQuery(con, statement = sql)
    dbDisconnect(con)
    dbUnloadDriver(drv)	     
  
    nighttime_coeff <- 0.5
    daytime_coeff <- 0.3
    coeff<-ifelse(nighttime,nighttime_coeff,daytime_coeff)
    # load geocells from db
    tmpdata<-readOGR(dsn=dbin$dbconn,layer=dbin$layer,verbose=FALSE)  #"urbanstructures")
    
    # project into WGS84
    tmpdata<-spTransform(tmpdata,CRS('+proj=longlat +datum=WGS84 +no_defs'))
    # filter out non-built-up and industrial/commercial areas
    geocells<-subset(tmpdata,tmpdata@data$structure >0)
    
    #compute centroids of geocells
    geocells_centroids<-coordinates(geocells)
    
    # evaluate intensity for each settlement using the shakemap
    intens<-round(extract(shakemap,geocells_centroids,method='bilinear'))
    
    #### Michael: check added to prevent unneeded error messages
    #check if there is any loss at all --> avoid NA
    check <- intens
    check[which(is.na(intens))] <-0
    
    if(max(check)==0||max(check)<=5){
      # Return safely from calculation if no loss
      msg <- 'With this configuration the model expects no loss!'
      print(msg)
      # Marc: clean urbanstructures_loss table in database
      sql <- "delete from urbanstructures_loss;" 
      drv <- dbDriver("PostgreSQL")
      con <- dbConnect(drv, host='lhotse21.gfz-potsdam.de', user= "postgres", password="postgres", dbname="losstool")
      post <- dbGetQuery(con, statement = sql)
      dbDisconnect(con)
      dbUnloadDriver(drv)
      return(list(tot_fatalities=0,msg))       
    }else{
      #keep only intensities >5
      ind <- which(intens >5)
    }
    
    #keep only intensities >5
    ind<-which(intens>5)
    
    # get vulnerability composition for all strata
    vuln_strata<-vuln_strata()
    
    # select uncertainty level. type=3 refers to most likely fragility values
    type<-3
    # init DPMs
    dpm_a<-ems98_DPM('a',type)
    dpm_b<-ems98_DPM('b',type)
    dpm_c<-ems98_DPM('c',type)
    dpm_d<-ems98_DPM('d',type)
    dpm_e<-ems98_DPM('e',type)
    dpm_f<-ems98_DPM('f',type)
    
    #service function to process each geocell
    process_geocells<-function(cell)
    {
      stratum<-geocells[cell,]@data$structure
      nb<-round(geocells[cell,]@data$nbuildings_mean)
      ncel<-0
      ncol<-dpm_a[intens[cell]-4,'d5']*vuln_strata[stratum,1]
      ncol<-ncol+dpm_b[intens[cell]-4,'d5']*vuln_strata[stratum,2]
      ncol<-ncol+dpm_c[intens[cell]-4,'d5']*vuln_strata[stratum,3]
      ncol<-ncol+dpm_d[intens[cell]-4,'d5']*vuln_strata[stratum,4]
      ncol<-ncol+dpm_e[intens[cell]-4,'d5']*vuln_strata[stratum,5]
      ncol<-ncol+dpm_f[intens[cell]-4,'d5']*vuln_strata[stratum,6]
      
      ncol<-round(ncol*nb)
      #average occupancy
      avg_occ <- 10
      #estimated fatalities, using modifier coefficient
      nfat<-ncol*avg_occ*coeff
    }
    
    #for each geocell estimate fatalities
    est_fat<-(sapply(ind,FUN=process_geocells))
    
    if(sum(est_fat)<1){
    msg <- 'With this configuration the model expects no loss!'
    print(msg)
    return(list(tot_fatalities=0,msg))
    }
    
    df2<-data.frame(geocells[ind,]@data,est_fat)
    df3<-SpatialPolygonsDataFrame(geocells[ind,],df2)
    
    writeOGR(df3,dbout$dbconn,dbout$layer, "PostgreSQL",overwrite_layer=dbout$overwrite)
    
    #uncomment to save them as a shapefile
    #path<-'/home/max/Documents/GFZ_sync/workspace/yurta/Loss'
    #writeOGR(df3,path,"geocells_loss", driver="ESRI Shapefile")
    msg <- 'Loss for geocells estimated successfully'
    msg2 <- sprintf("%s%i",'Number of total fatalities estimated: ',sum(est_fat))
    print(msg)
    print(msg2)
    return(list(tot_fatalities=sum(est_fat),msg))
  }
  
  ###########################################################################################################
  # VULNERABILITY COMPUTATION
  ###########################################################################################################
  
  # get vulnerability for each geocell stratum
  vuln_strata<-function()
  {
    # load average vulnerability distributions
    # from postgresql db
    sql <- "select structure as s, unnest(post_vuln) as vuln from avg_posteriors;" 
    
    drv <- dbDriver("PostgreSQL")
    con <- dbConnect(drv, host='lhotse21.gfz-potsdam.de', user= "postgres", password="postgres", dbname="losstool")
    post <- dbGetQuery(con, statement = sql)
    
    # Closes the connection
    dbDisconnect(con)
    # Frees all the resources on the driver
    dbUnloadDriver(drv)	
    
    # return avg vulnerability distribution of all strata 
    #post$vuln[post$s==stratum]
    t(sapply(1:21,FUN=function(x){post$vuln[post$s==x]}))
  }
  
  # Damage Probability Matrix (DPM) for a specific vulnerability class
  #type: 1: least possible Vulnerability
  #      2: least plausible "
  #      3: most likely "
  #      4: greatest plausible "
  #      5: greatest possible "
  ems98_DPM<-function(ems,type)
  {
    V<-ems98_C2V(ems,type)
    #init dpm
    mus<-sapply(seq(5,10),function(x) rep(ems98_DG(x,V),6))
    ind<-rep(seq(0,5),6)
    d<-sapply((1:36),function(x) P_DS(ind[x],mus[x]))
    dpm<-matrix(d,nrow=6,ncol=6,byrow=TRUE,
                dimnames=list(c('V','VI','VII','VIII','IX','X'),c('d0','d1','d2','d3','d4','d5')))
    dpm
  }
  
  # probability of discrete damage states by integration of prob. density
  P_DS<-function(dg,mu_d,t=8)
  {
    integrate(function(x)pbeta(x,mu_d,t),dg,dg+1)$value
  }
  
  # beta function used to define damage state probability density
  pbeta<-function(x,mu_d,t=8)
  {
    a<-0
    b<-6
    r<-t*(0.007*mu_d^3-0.0525*mu_d^2+0.2875*mu_d)
    p1<-gamma(t)/(gamma(r)*gamma(t-r))
    p2<-p1*((x-a)^(r-1)*(b-x)^(t-r-1))/((b-a)^(t-1))
    p2
  }
  
  # vulnerability index for a specific EMS98 vulnerability class
  #type: 1: least possible
  #      2: least plausible
  #      3: most likely
  #      4: greatest plausible
  #      5: greatest possible
  ems98_C2V<-function(ems,type=3)
  {
    f<-c(0.22,0.14,0.1,0.06,-0.08)
    e<-c(0.38,0.3,0.26,0.22,0.14)
    d<-c(0.54,0.46,0.4,0.38,0.3)
    c<-c(0.7,0.62,0.58,0.54,0.46)
    b<-c(0.86,0.78,0.74,0.7,0.62)
    a<-c(1.02,0.94,0.9,0.86,0.78)
    vmat<-rbind(a,b,c,d,e,f)
    vmat[ems,type][[1]]
  }
  
  # damage grade for a specific MMI and vulnerability index (and ductility)
  ems98_DG<-function(I,V,Q=2.3)
  {
    res <- 2.5*(1+tanh((I+6.25*V-13.1)/Q))
    res
  }
  
  
  ###########################################################################################################
  # Landslide (based on Annamarias work)
  ###########################################################################################################
  
  getLS <- function(shakemap){    
    #read in Annamarias susceptibility map computed from topography, geo-tectonics and land_use
    inScpt <- raster("/var/www/losstool/apps/landslides/LSI3_gc.tif")
    #check for overlap
    #crop to extent of szenario (shakemap)
    suppressWarnings(tryCatch(inScptAOI <- intersect(inScpt,shakemap),error=function(e) 'Landslide Susceptibility: Modeled area outside Kyrgyzstan'))
    
    ##calculate a weight layer for the intensities
    #fun1 <- function(x){x[x>0 && x<8]<-x*(-1.62)}
    #weight <- calc(shakemap,fun1)
    #
    #maxInt <- cellStats(shakemap,stat='max')
    #if (maxInt > 8){
    #fun2 <- function(x){x[x>8 && x < 8.09999] <- x*(-0.754)}
    #weight <- calc(weight,fun2) 
    #  if (maxInt > 8.09999){ 
    #  fun3 <- function(x){x[x>8.09999]<-x*1.467}
    #  weight <- calc(weight,fun3)
    #  }
    #}
    ##multiply weights with susceptibility
    #landslide <- overlay(inScptAOI,weight,fun=function(x,y){return(x*y)})
    ##apply threshold for landslide
    #landslide <- calc(landslide, fun=function(x){ x[x < 5.659001] <- NA; return(x)} )
    
    #####################
    # Keep only where int larger 8
    
    if(exists("inScptAOI")&&(cellStats(shakemap,stat='max')>8)){
      inScptAOI <- resample(inScptAOI,shakemap,method="bilinear")
      largerInt8 <- calc(shakemap, fun=function(x){ x[x <= 8] <- NA ; return(x)})
      if(exists("largerInt8")){
        landslide <- mask(inScptAOI,largerInt8)
      }
    }else{
      landslide <- raster(shakemap) #empty raster
    }
    ######################
    
    #write file
    writeRaster(landslide,filename="/var/www/losstool/apps/landslides/landslide.tif",
                format="GTiff",overwrite=TRUE)
    
  } #end getLS
 
  ###########################################################################################################
  # Function starting the actual calculations 								    #
  ###########################################################################################################

  getLoss <- reactive({
    # get shakemap 
    shakemap <- getInt()
    # get landslide binary map
    landslide <- getLS(shakemap)
    dbconn2<-"PG:dbname=losstool host='lhotse21.gfz-potsdam.de' user='postgres' password='postgres'"
    type <- as.integer(input$type)
    if (type==1)
    {
      dbin<-list(dbconn=dbconn2,layer='places',overwrite=FALSE)
      dbout<-list(dbconn=dbconn2,layer='places_loss',overwrite=TRUE)	
      result <- compute_loss_points(shakemap,nighttime=FALSE,dbin,dbout)
    } 
    else if (type==2)
    {
      # test geocell losses
      dbin<-list(dbconn=dbconn2,layer='urbanstructures',overwrite=FALSE)
      dbout<-list(dbconn=dbconn2,layer='urbanstructures_loss',overwrite=TRUE)
      result <- compute_loss_geocells(shakemap,nighttime=FALSE,dbin,dbout)
    }
  #}
  return(msg<-tail(result,1))
  }) #end getLoss
  
  
  #########################################################################
  #########################################################################
  ##                                                                     ##  
  ##                    B)    Plotting Routines                          ##        
  ##                                                                     ##
  #########################################################################
  #########################################################################
  
  ################################################################
  # Plot intensities                                             #
  ################################################################
  
  output$msg <- renderPrint({
  
    #calculate loss
    msg <-getLoss()
    withProgress(session, min=1, max=10, {
      setProgress(message = 'Calculation in progress', detail = 'This may take a while...')
      for (i in 1:10) {
        setProgress(value = i)
        Sys.sleep(0.5)
      }
    })
  }) #end renderPlot
}) #end server.R
