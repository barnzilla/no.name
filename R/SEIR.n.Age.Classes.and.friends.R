#' SEIR model computations
#' @export
SEIR.n.Age.Classes = function(file.name, sheets.names, differential.eqns.func=NULL, post.processing.func=NULL ,post.processing.companion.kit=NULL, also.get.flows=NULL)
{
  #sheets.names is list(parms.1d,  parms.2d,  initial.conditions,  model.flow,  auxiliary.vars)

  input.info.verbatim = sheets.names
  input.info.verbatim$file.name = file.name
  #==========================================================================
  #   Initial conditions and parameters
  #==========================================================================
  #browser()

  tmp_time_stuff   = sheets.names$parms.1d            # 1 dimensional parameters
  tmp_time_stuff_m = sheets.names$parms.2d            # 2 dimensional parameters
  input_stuff      = sheets.names$initial.conditions  # initial values/conditions

  if(!is.data.frame(tmp_time_stuff))
    tmp_time_stuff   <- tbl_to_df( readxl::read_excel(file.name, sheet = tmp_time_stuff) )  # 1 dimensional parameters

  if(!is.data.frame(tmp_time_stuff_m))
    tmp_time_stuff_m <- tbl_to_df( readxl::read_excel(file.name, sheet = tmp_time_stuff_m) ) # 2 dimensional parameters

  if(!is.data.frame(input_stuff))
    input_stuff  <- tbl_to_df( readxl::read_excel(file.name, sheet = input_stuff) ) # initial values/conditions

  raw.init.conditions = input_stuff # keep snapshot for output.  input_stuff may be altered later

  nagegrp <- length(unique(tmp_time_stuff$agegrp))        # number of age groups
  nagegrp = ncol(raw.init.conditions) - 1             # number of age groups
  agegrp.suffix = ""
  if(nagegrp > 1)
    agegrp.suffix = 1:nagegrp

  raw.compartments =   raw.init.conditions$NAME # e.g.  S D L
  raw.compartments.age = c ( t( outer( raw.compartments,agegrp.suffix,paste0) ) ) # e.g. S1 S2 S3 S4 S5 D1 D2 ...


  nrow_   <- dim(tmp_time_stuff)[1]/nagegrp

  time_stuff <- dplyr::arrange(tmp_time_stuff, tmin, agegrp) # sort by tmin agegrp
  time_stuff <- time_stuff %>%
    mutate(isim = rep(1:nrow_, each=nagegrp))

  time_stuff_m <- dplyr::arrange(tmp_time_stuff_m, tmin, cagegrp, ragegrp) # sort by tmin cagegrp  ragegrp
  time_stuff_m <- time_stuff_m %>%
    mutate(isim = rep(1:nrow_, each = nagegrp*nagegrp))

  #===================================================================
  # Build function eval.post.processing.func (if not provided)
  #===================================================================

  code.body.df = sheets.names$post.processing # data.frame or string
  if(!is.data.frame(code.body.df))
    code.body.df = tbl_to_df( readxl::read_excel(file.name, sheet = code.body.df ) )

  code.body.df$code[is.na(code.body.df$code)] = ""

  eval.post.processing.func = post.processing.func
  if(is.null(eval.post.processing.func))
  {
    code.body.char = code.body.df$code
    code.head.char = c("function(list.solution.etcetera){",
                       "# list.solution.etcetera contains $solution and possibly $solution.inflows and/or $solution.outflows" ,
                       "# Probably could use    with(list.solution.etcetera)   instead of copying" ,
                       "solution          = list.solution.etcetera$solution",
                       "solution.inflows  = list.solution.etcetera$solution.inflows",
                       "solution.outflows = list.solution.etcetera$solution.outflows",
                       "parms.1d          = list.solution.etcetera$input.info$parms.1d",
                       "parms.2d          = list.solution.etcetera$input.info$parms.2d",
                       "companion.kit     = list.solution.etcetera$post.processing.companion.kit",
                       "sommaire          = data.frame(osqvhfipumzcdjblkwrgnexaty=4)" ,
                       "## BEGIN code from excel"  )
    code.tail.char = c("## END code from excel",
                       "sommaire$osqvhfipumzcdjblkwrgnexaty=c()",
                       "list.solution.etcetera$solution = solution # only solution may be modified",
                       "list.solution.etcetera$sommaire = sommaire # add on new kid on the block",
                       "list.solution.etcetera",
                       "#list(sommaire=sommaire,solution.inflows=solution.inflows,solution.outflows=solution.outflows)" )

    eval.post.processing.func = paste(c(code.head.char , code.body.char , code.tail.char , "}"),collapse="\n" )
    eval.post.processing.func = eval(parse(text= eval.post.processing.func )) # from text to function
  }

  #===================================================================
  # Build function eval.differential.eqns.func (if not provided)
  #===================================================================
  model_flows_tmp = auxiliary.vars = "Not used since was argument differential.eqns.func was provided"
  eval.differential.eqns.func = differential.eqns.func
  flows.requested = !is.null(also.get.flows) && any(also.get.flows %in% c("inflows","outflows"))
  if(flows.requested || is.null(eval.differential.eqns.func))
  {
    model_flows_tmp = sheets.names$model.flow
    auxiliary.vars  = sheets.names$auxiliary.vars

    if(!is.data.frame(model_flows_tmp))
      model_flows_tmp  <- tbl_to_df( readxl::read_excel(file.name, sheet = model_flows_tmp) )  # arrows in flowchart

    if(!is.data.frame(auxiliary.vars))
      auxiliary.vars   <- tbl_to_df( readxl::read_excel(file.name, sheet = auxiliary.vars ) )  # auxiliary variables in equations

    model_flows_tmp = subset( model_flows_tmp,!is.na(expression) & activation ==1)
    model_flows = model_flows_tmp     # keep model_flows_tmp for output. Work on model_flows from this point on


    # BEGIN update model_flows and input_stuff as per requested flows (if any)

    init.cond.numeric.vars = setdiff(colnames(input_stuff),"NAME")
    silly.inflow.model = silly.outflow.model = list()

   #if(!is.null(flows.prefix$inflows))
    if(any(also.get.flows == "inflows"))
      silly.inflow.model = get.silly.model.chunk(model_flows, input_stuff, init.cond.numeric.vars, "inflow" ,  "SlackBoxForInflows", "inflows.")

   #if(!is.null(flows.prefix$outflows))
    if(any(also.get.flows == "outflows"))
      silly.outflow.model= get.silly.model.chunk(model_flows, input_stuff, init.cond.numeric.vars, "outflow", "outflows.", "SlackBoxForOutflows")

    model_flows = rbind(model_flows, silly.inflow.model$model, silly.outflow.model$model)
    input_stuff = rbind(input_stuff, silly.inflow.model$init , silly.outflow.model$init )

    # END   update model_flows and input_stuff as per requested flows (if any)


    # Get init_list.  Same content as input_stuff but in a list format
    init_list <- list()
    for(k in input_stuff$NAME)
      init_list[[k]] <- as.matrix( subset(input_stuff, NAME == k)[,init.cond.numeric.vars] )


    # BEGIN Build character vector differential.eqns.char
    model_flows$multiplier = model_flows$From

    lazy.range=range(model_flows$lazy)
    if(is.na(lazy.range[1]))
      stop("\nstyle (lazy or not) must be specified for all transitions")
    if(diff(lazy.range)>0 || !(lazy.range[1] %in% c(0,1)) )  # may be relaxed later
      stop("\nstyle (lazy or not) should be all zeros or all ones")
    if(flows.requested && lazy.range[1]==1)
      stop("\nLazy style should not be used when inflows and/or outflows are requested.")  # may be relaxed later

    if(lazy.range[1] == 0)  # not lazy
      model_flows$multiplier = NA

    model_flows = model_flows %>%
      mutate(fromstar = ifelse(is.na(multiplier),"",paste0(From,"*(")) ,
             closing.parenthesis = ifelse(is.na(multiplier),"",")") ,
             expression = paste0(fromstar,expression,closing.parenthesis) )

    boxes = unique(c(model_flows$From,model_flows$To))
    outflows=rep("#####$$$$$$$$$$$!@@@@@@@@",length(boxes))
    names(outflows) = boxes
    inflows=outflows

    for(k in seq(nrow(model_flows)))
    {
      this.expression = model_flows$expression[k]
      outflows[model_flows$From[k]] = paste(outflows[model_flows$From[k]],"+", this.expression)
      inflows [model_flows$To[k]  ] = paste( inflows[model_flows$To[k]  ],"+", this.expression)
    }
    differential.eqns.char = paste("(flow.multiplier$inflow)*(",inflows,") + (flow.multiplier$outflow)*(",outflows,")")
    differential.eqns.char = gsub(  "#####$$$$$$$$$$$!@@@@@@@@ +","" ,differential.eqns.char ,fixed=TRUE)
    differential.eqns.char = gsub("( #####$$$$$$$$$$$!@@@@@@@@ )","0",differential.eqns.char ,fixed=TRUE)
    differential.eqns.char = paste0("derivative.",boxes,"=",differential.eqns.char)

    #if( any(grepl("#####$$$$$$$$$$$!@@@@@@@",differential.eqns.char)) )  browser()


    rm(model_flows,boxes,inflows,outflows,this.expression) # clean up intermediate objects no longer required

    # END Build character vector differential.eqns.char


    # BEGIN Build function eval.differential.eqns.func

    auxiliary.vars.assignment = paste(auxiliary.vars$variable,"=",auxiliary.vars$expression)
    compartments =   input_stuff$NAME # e.g.  S D L
    compartments.age = c ( t( outer( compartments,seq(nagegrp),paste0) ) ) # e.g. S1 S2 S3 S4 S5 D1 D2 ...
    out.char      = paste0("out=c(", paste0("derivative.", paste(compartments    , collapse=", derivative.")),")")
    names.out.char= paste0("out=c(", paste0("derivative.", paste(compartments.age, collapse=",derivative." )),")")
    names.out.char = gsub(",","','",names.out.char ,fixed=TRUE) # replace , by ','
    names.out.char = gsub("(","('" ,names.out.char ,fixed=TRUE) # replace ( by ('
    names.out.char = gsub(")","')" ,names.out.char ,fixed=TRUE) # replace ) by ')
    names.out.char = gsub("out=","names(out)=",names.out.char ,fixed=TRUE)

    #function.char = c(auxiliary.vars.assignment,differential.eqns)
    within.with.part.char = paste(c(auxiliary.vars.assignment,differential.eqns.char,out.char,names.out.char,"list(out)"),collapse="\n") # First tried ";\n" but realised that "\n" is sufficient

    function.char = "function(list.inits, list.parms,flow.multiplier=list(inflow=1,outflow=-1)){"
    with.char =   " with(as.list(c(list.inits, list.parms)),{"


    eval.differential.eqns.func = paste(c(function.char,with.char,within.with.part.char,"})" ,"}"),collapse="\n" ) # code of eval.differential.eqns.func
    # cat(eval.differential.eqns.func) # quick look at how it look like
    eval.differential.eqns.func = eval(parse(text =  eval.differential.eqns.func ))

    #eval.differential.eqns.func = eval(parse(text= "function(X) {\n cos(X)\n}" ))
    #eval.differential.eqns.func = eval(parse(text= "function(X) { cos(X)}" ))
    #eval.differential.eqns.func (pi)  # gives -1 as it should

    rm(auxiliary.vars.assignment, out.char, within.with.part.char, function.char, with.char) # clean up
    # END Build function eval.differential.eqns.func
  }


  #==========================================================================
  #  Main routine
  #==========================================================================
  # The SEIR model with N age classes
  #
  #' SEIR model computations within N age groups
  #' @export
  SEIR.n.Age.Classes.within.loop <- function( time=NULL, age.parms = NULL, age.age.parms = NULL,list.inits = NULL, not.parms=  c("tmin", "tmax", "agegrp", "cagegrp", "ragegrp", "isim"))
  {
    nage = nrow(age.parms)

    if (is.null(age.parms))
      stop("undefined 'age.parms'")

    if (is.null(age.age.parms))
      stop("undefined 'age.age.parms'")

    if (is.null(time))
      stop("undefined 'time'")

    if (is.null(list.inits))
      stop("undefined 'list.inits'")


    list.parms <- list()
    for(k in setdiff(names(age.parms), not.parms )){
      list.parms[[k]] <- age.parms[,k]
    }
    for(k in setdiff(names(age.age.parms), not.parms ))
    {
      temp<- array(NA, c(nage,nage))
      temp[cbind(age.age.parms$cagegrp, age.age.parms$ragegrp)] <- age.age.parms[,k]
      list.parms[[k]] <- temp
      if(any(is.na(temp)))
        stop(paste0(k," matrix has some missing entries"))
    }

    ### to write the system of differential equations
    calculate_derivatives <- function(time, vec.inits, list.parms, names.inits) {


      iota <- seq(length(vec.inits)/length(names.inits))
      list.inits <- list()
      for(k in names.inits){
        if (length(iota) > 1) {
          list.inits[[k]] <- vec.inits[paste0(k, iota)]
        }else{list.inits[[k]] <- vec.inits[k] }
      }

      eval.differential.eqns.func(list.inits, list.parms)


    } #end of function calculate_derivatives

    ###---------------------------------------------------------------------------------
    ### Solver for Ordinary Differential Equations (ODE), Switching Automatically
    ### Between Stiff and Non-stiff Methods
    ###---------------------------------------------------------------------------------

    output <-  lsoda(y = unlist(list.inits),
                     times = time,
                     func =  calculate_derivatives,
                     parms = list.parms,
                     names.inits = names(list.inits))

    return(output)
  }  # END  of function SEIR.n.Age.Classes.within.loop

  ################################################################
  #        To run the example
  ################################################################

  excluded_names <- c("tmin", "tmax","agegrp","cagegrp","ragegrp","isim")

  nSim <- max(time_stuff$isim)
  listOut <- list()
  previous.tmax <- 0
  out<-NULL
  df.out = c()

  for(i in seq(1, nSim, 1))
  {
    parameter.by.age     <- subset(time_stuff  , isim == i)
    parameter.by.age.age <- subset(time_stuff_m, isim == i)

    tmin <- unique(c(parameter.by.age$tmin, parameter.by.age.age$tmin))
    tmax <- unique(c(parameter.by.age$tmax, parameter.by.age.age$tmax))

    if(length(tmin)>1 || length(tmax)>1 || tmin>=tmax )
      stop(paste0("Unexpected pattern in tmin, tmax for interval ", i))

    tt <- seq(0, tmax - tmin, by = 1)

    if(tmin != previous.tmax)
      stop(paste(interval.label , "\n  Interval lower bound not equal to previous interval upper bound"))

    previous.tmax <- tmax
    out <- SEIR.n.Age.Classes.within.loop( time=tt,
                                           age.parms = parameter.by.age,
                                           age.age.parms = parameter.by.age.age,
                                           list.inits = init_list)

    out <- as.data.frame(out)
    # ode/lsoda Output diagnostic #######################
    #diagn <- diagnostics.deSolve(out)

    out$time <- seq(tmin,tmax,1)
    out_for_init <- out %>%
      slice(nrow(out)) %>%
      pivot_longer(-time)
    init <- out_for_init$value
    names(init) <- out_for_init$name

 #   rowns <- names(select(out,-c(time)))
 #   out <- out %>%
 #     mutate(N_tot = rowSums(.[rowns]))  # Total number of individuals


    #updating the initial values
    for(k in 1:length(init_list)){
      init_list[[k]][1:nagegrp] <- init[seq(nagegrp*(k-1)+1,nagegrp*k)]
    }

   if(i < nSim)
      out = out[-nrow(out),]
    # Add outputs to the list
    df.out = rbind(df.out,out)
    listOut[[i]] <- out
  } #end for(i in seq(1, nSim, 1))

  resultat = list( differential.eqns.func = eval.differential.eqns.func ,
                   eval.post.processing.func = eval.post.processing.func,
                   post.processing.companion.kit = post.processing.companion.kit)

  #BEGIN add on $solution.inflows, $solution.outflows and $solution to resultat
  solution = df.out
  if(any(also.get.flows == "inflows"))
  {
    var.names = c ( t( outer( paste0("inflows.",raw.compartments),agegrp.suffix,paste0) ) )
    resultat$solution.inflows  = solution[,c("time",var.names)]
    colnames( resultat$solution.inflows) = c("time",raw.compartments.age)
  }
  if(any(also.get.flows == "outflows"))
  {
    var.names = c ( t( outer( paste0("outflows.",raw.compartments),agegrp.suffix,paste0) ) )
    resultat$solution.outflows  = solution[,c("time",var.names)]
    colnames( resultat$solution.outflows) = c("time",raw.compartments.age)
  }
  resultat$solution = solution[,c("time",raw.compartments.age)]
  #END add on $solution.inflows, $solution.outflows and $solution to resultat

  resultat$input.info = list(parms.1d=tmp_time_stuff, parms.2d=tmp_time_stuff_m, initial.conditions=raw.init.conditions,
                             auxiliary.vars=auxiliary.vars, model.flow=model_flows_tmp, post.processing=code.body.df)
  resultat$input.info.verbatim = input.info.verbatim

  eval.post.processing.func(resultat)
} #end of SEIR.n.Age.Classes function

#' Get SEIR model inputs
#' @export
get.silly.model.chunk = function(model_flows, init.cond, init.cond.numeric.vars, which.flow, NewFrom, NewTo)
{
  #browser()
  inflow.orphans = outflow.orphans = c()
  silly_model = model_flows
  if(which.flow %in% c("inflow","inflows"))
  {
    silly_model$To   = paste0(NewTo,silly_model$To)
    silly_model$From = NewFrom
    init.prefix      = NewTo
    init.slack.name  = NewFrom
    inflow.orphans   = setdiff(init.cond$NAME,unique(model_flows$To))
  }
  if(which.flow %in% c("outflow","outflows"))
  {
    silly_model$From = paste0(NewFrom,silly_model$From)
    silly_model$To   = NewTo
    init.prefix      = NewFrom
    init.slack.name  = NewTo
    outflow.orphans  = setdiff(init.cond$NAME,unique(model_flows$From))
  }

  silly_init  = init.cond
  silly_init$NAME  = paste0(init.prefix,silly_init$NAME)
  silly_init[,init.cond.numeric.vars] = 0
  init.slack = silly_init[1,]
  init.slack$NAME = init.slack.name
  silly_init = rbind(silly_init,init.slack)

  #BEGIN handle the orphans
  model.inflow.orphans = model.outflow.orphans = c()
  if(length(inflow.orphans) > 0)
  {
    model.inflow.orphans = silly_model[rep(1,length(inflow.orphans)),]
    model.inflow.orphans$To = paste0(NewTo,inflow.orphans)
    model.inflow.orphans$expression = "0"
  }
  if(length(outflow.orphans) > 0)
  {
    model.outflow.orphans = silly_model[rep(1,length(outflow.orphans)),]
    model.outflow.orphans$From = paste0(NewFrom,outflow.orphans)
    model.outflow.orphans$expression = "0"
  }
  #END handle the orphans

  list(init=silly_init , model = rbind(silly_model,model.inflow.orphans,model.outflow.orphans) )
}

#' Get other SEIR model inputs
#' @export
smooch.parms.df.into.list = function(df.parms.1d,df.parms.2d,these.are.not.parms)
{ # takes 1d and 2d parms in data.frame format and makes a list.
  # for 2d parms, they are in skinny format in input data frame, but are square matrices in the output list.
  list.parm.1d = as.list(df.parms.1d [,setdiff(colnames(df.parms.1d),these.are.not.parms)]) # list rather than data.frame
  list.parm.2d = as.list(df.parms.2d [,setdiff(colnames(df.parms.2d),these.are.not.parms)]) # list rather than data.frame

  # BEGIN fix list.parm.2d which is not quite what we want.  We want list of matrices, not vectors.
  list.parm.2d.tmp = list.parm.2d
  list.parm.2d = list()
  matrix.template = array(NA,c(1,1)*max(list.parm.2d.tmp$ragegr))
  matrix.indices = cbind(list.parm.2d.tmp$ragegrp,list.parm.2d.tmp$cagegrp)
  for(k in setdiff(names(list.parm.2d.tmp),c("ragegrp" , "cagegrp")) )
  {
    list.parm.2d[[k]] = matrix.template
    list.parm.2d[[k]] [matrix.indices] = list.parm.2d.tmp[[k]]
  }
  # END fix list.parm.2d which is not quite what we want.  We want list of matrices.
  as.list(c(list.parm.1d,list.parm.2d))
}


#' Ever been here function
#' @export
ever.been.here = function(SEIR.object)
{
  n_agegrp = ncol(SEIR.object$input.info$initial.conditions) - 1 # number of age groups
  agegrp.suffix = ""
  if(n_agegrp > 1)
    agegrp.suffix = 1:n_agegrp

 # SEIR.object = SEIR.object.5age.betavec.not.lazy
  range(diff(SEIR.object$solution$time)) # check that times are in increments of 1 (feature implicitly used below)
  differential.eqns.func=SEIR.object$differential.eqns.func

  these.are.not.parms = c("tmin",     "tmax" ,  "agegrp", "isim" ) # actaully want to keep "ragegrp",  "cagegrp"
  inflows  = c()
  outflows = c()
  for(tt in SEIR.object$solution$time) # [1:4]
  {
    cat("\nicitte ",tt)
    #BEGIN get list.parm.1d.2d  (all parameters stored in a list  )
    df.parm.1d = SEIR.object$input.info$parms.1d
    df.parm.2d = SEIR.object$input.info$parms.2d
    df.parm.1d = subset(df.parm.1d ,overlap.length(tt,tt+1,df.parm.1d$tmin,df.parm.1d$tmax) > 0 )
    df.parm.2d = subset(df.parm.2d ,overlap.length(tt,tt+1,df.parm.2d$tmin,df.parm.2d$tmax) > 0 )
    list.parm.1d.2d = smooch.parms.df.into.list(df.parm.1d,df.parm.2d,these.are.not.parms)
    #END get list.parm.1d.2d  (all parameters stored in a list  )

   # browser()
    #BEGIN list.compartments  ... will contain vectors D, R, S, ... at time tt
    list.compartments = list() # will contain vectors D, R, S, ....
    mat1row.compartments = data.matrix(subset(SEIR.object$solution,time == tt))
    for(k in SEIR.object$input.info$initial.conditions$NAME)
      list.compartments[[k]] = mat1row.compartments[1,paste0(k,agegrp.suffix)]
    #END list.compartments  ... will contain vectors D, R, S, ... at time tt
    # browser()
    inflows  = rbind(inflows  , differential.eqns.func(list.compartments,list.parm.1d.2d,flow.multiplier =list(inflow = 1,outflow=0))[[1]] )
    outflows = rbind(outflows , differential.eqns.func(list.compartments,list.parm.1d.2d,flow.multiplier =list(inflow = 0,outflow=1))[[1]] )
  }
  browser()
  colnames(inflows ) = gsub("derivative.","",colnames(inflows ))
  colnames(outflows) = gsub("derivative.","",colnames(outflows))
  solution.lean = SEIR.object$solution[,setdiff(names(SEIR.object$solution),c("time","N_tot")) ]


  ever.been.here.from.inflows = rbind(solution.lean[1,],inflows[-nrow(inflows),])
  ever.been.here.from.inflows = as.data.frame( apply(ever.been.here.from.inflows,2,cumsum)  )
  ever.been.here.from.inflows$time = SEIR.object$solution$time

  ever.been.here.from.outflows = rbind(0*solution.lean[1,],outflows[-nrow(outflows),])
  ever.been.here.from.outflows = as.data.frame( apply(ever.been.here.from.outflows,2,cumsum)  )
  ever.been.here.from.outflows = ever.been.here.from.outflows[,colnames(solution.lean)] + solution.lean
  ever.been.here.from.outflows$time = SEIR.object$solution$time

  list(ever.been.here.from.inflows=ever.been.here.from.inflows, ever.been.here.from.outflows=ever.been.here.from.outflows)
}

if(FALSE)
{
  # Compare 1 age group results with 5 age groups
  results.1age  = "provide some object"
  results.5ages = "provide some object"

  mat.range = c()
  vars.to.check = setdiff(colnames(results.1age$solution),c( "time","N_tot" ))
  for(this.box in vars.to.check )
  {
    somme = apply(results.5ages$solution[,paste0(this.box,1:5)],1,sum)
    actual.diff = results.1age$solution[,this.box] - somme
    rel.diff = actual.diff /( 1e-9 +  pmax(results.1age$solution[,this.box] , somme) )
    this.range = c(range(rel.diff ) , range(actual.diff))
    mat.range = rbind(mat.range,this.range)
    cat("\n",this.box,"\t", this.range)
  }
  rownames(mat.range) = vars.to.check
  matplot(seq(nrow(mat.range)),mat.range[,1:2],main='Relative differences')
  matplot(seq(nrow(mat.range)),mat.range[,3:4],main='Actual   differences')
}

#' Save SEIR model to Excel file
#' @export
SaveModelInExcel = function(input.info.list,file_name,map.names)
{  # map.names is named vector saying map.names["parms.1d"  ] =  "Parameters by Age"
 # map.names = unlist(sheet_names) # allows to map "parms.1d" into "Parameters by Age" for instance
  names(input.info.list) = map.names[names(input.info.list)] # $parms.1d is now $'Parameters by Age'
  openxlsx::write.xlsx(input.info.list, file_name, colWidths = c(NA, "auto", "auto"))

}


#' Parameter sweep function
#' @export
try.various.parms.values = function(SEIR.object,parm.cloud.grid.specs,only.show.parms.to.try=FALSE)
{
  #parm.cloud.grid.specs is a list that should contain the following 7 things
  # *  $hypercube.lower.bounds , $hypercube.upper.bounds, $hypercube.apex.mode
  # *  $n.repeat.within.hypercube
  # *  $LatinHypercubeSampling
  # *  $racine
  # *  $tmin.alter.scope
  # *  $backend.transformation
  # *  $reference.alteration

  lower.bound      = parm.cloud.grid.specs$hypercube.lower.bounds
  upper.bound      = parm.cloud.grid.specs$hypercube.upper.bounds
  apex             = parm.cloud.grid.specs$hypercube.apex.mode  # may be NULL
  n.repeat         = parm.cloud.grid.specs$n.repeat.within.hypercube
  tmin.alter.scope = parm.cloud.grid.specs$tmin.alter.scope
# racine           = parm.cloud.grid.specs$racine
# backend.transformation = parm.cloud.grid.specs$backend.transformation
# reference.alteration      = parm.cloud.grid.specs$reference.alteration

  set.seed(parm.cloud.grid.specs$racine)

  operation_list = list(overwrite = function(current,new)         {0*current+new        } ,
                        add       = function(current,increment  ) {  current+increment  } ,
                        multiply  = function(current,mult_factor) {  current*mult_factor} )
  operation_func = operation_list[[parm.cloud.grid.specs$reference.alteration]]
  operation.label = c(overwrite=".overwrite",add=".add",multiply=".multiplier")[parm.cloud.grid.specs$reference.alteration]


  # Compute all possible candidate by multiplier combinations
  #multiplier_combos <- expand.grid(candidates, KEEP.OUT.ATTRS = FALSE)
  lower.bound.expanded = expand.grid(lower.bound, KEEP.OUT.ATTRS = FALSE)
  upper.bound.expanded = expand.grid(upper.bound, KEEP.OUT.ATTRS = FALSE)
  if(any(dim(lower.bound.expanded) != dim(upper.bound.expanded) ))
    stop('dimension mismatch')
  if(any(colnames(lower.bound.expanded) != colnames(upper.bound.expanded) ))
    stop('column name mismatch')

  lower.bound.expanded = lower.bound.expanded[rep(seq(nrow(lower.bound.expanded)) , n.repeat ),,drop=F]
  upper.bound.expanded = upper.bound.expanded[rep(seq(nrow(upper.bound.expanded)) , n.repeat ),,drop=F]
  lower.bound.expanded = data.matrix(lower.bound.expanded)
  upper.bound.expanded = data.matrix(upper.bound.expanded)
  if(!is.null(apex))
  {
    apex.expanded = expand.grid(apex, KEEP.OUT.ATTRS = FALSE)
    apex.expanded = apex.expanded[rep(seq(nrow(apex.expanded)) , n.repeat ),]
    apex.expanded = data.matrix(apex.expanded)
  }

  if(parm.cloud.grid.specs$LatinHypercubeSampling)
    random.vect= c( lhs::randomLHS( nrow(upper.bound.expanded),ncol(upper.bound.expanded) ) )
  else
    random.vect= runif(prod(dim(upper.bound.expanded)))

  # parms.to.try = (upper.bound.expanded - lower.bound.expanded) *  random.vect
  # parms.to.try = lower.bound.expanded + parms.to.try
  if(is.null(apex))
    parms.to.try =    stats::qunif    (random.vect,lower.bound.expanded,upper.bound.expanded)
  else
    parms.to.try = triangle::qtriangle(random.vect,lower.bound.expanded,upper.bound.expanded,apex.expanded)

  parms.to.try = 0*lower.bound.expanded + parms.to.try
  parms.to.try = parm.cloud.grid.specs$backend.transformation(parms.to.try)

  parms.to.try.verbose = as.data.frame(parms.to.try)
  colnames(parms.to.try.verbose) = paste0(colnames(parms.to.try),operation.label)
  rownames(parms.to.try.verbose) = c()

  if(only.show.parms.to.try)
    return (list(parms.to.try=parms.to.try.verbose))


  # Recover some info from baseline/template (i.e. SEIR.object)

  sheet_names_for_sweep = SEIR.object$input.info           # Can use either of those two lines ... in theory (not tested)
  sheet_names_for_sweep = SEIR.object$input.info.verbatim  # Can use either of those two lines ... in theory (not tested)

  baseline.parms.1d = SEIR.object$input.info$parms.1d  # data frame of 1d parameters
  baseline.parms.2d = SEIR.object$input.info$parms.2d  # data frame of 2d parameters

  post.process.kit  = SEIR.object$post.processing.companion.kit

  flows.of.interest = gsub("solution.","",intersect(names(SEIR.object),c("solution.inflows","solution.outflows")))


  list.sweep = list() #        store results in list  ...
  df.sweep = c()      # ... or store results in data.frame
  list.sweep.ever.from.inflows = list.sweep.ever.from.outflows = list.sweep
  df.sweep.ever.from.inflows =   df.sweep.ever.from.outflows =   df.sweep  # later maybe
  outcomes.summary.df = c()
  ever.been.here.info= list(ever.been.here.from.inflows=0         , ever.been.here.from.outflows=0         ) # numeric   --> do not bypass
  ever.been.here.info= list(ever.been.here.from.inflows="Not done", ever.been.here.from.outflows="Not done") # character -->        bypass

  for(i in 1:nrow(lower.bound.expanded)) {
    row <-  parms.to.try[i,]
    names(row) = colnames(parms.to.try)
    #this.label <- paste0(names(row),       ".multiplier= ", row, collapse = " , ")
    this.label <- paste0(names(row), operation.label, "= ", row, collapse = " , ")


    cat("\n Doing",i,"of",nrow(parms.to.try),"simulations :",this.label)

    parms.1d = baseline.parms.1d  # data frame of 1d parameters to be altered
    parms.2d = baseline.parms.2d  # data frame of 2d parameters to be altered

    # Modify the parameter values at the specified tmin.alter.scope
    for(parameter in names(row)) {
      if(parameter %in% names(parms.1d)) {
        parms.1d[parms.1d$tmin %in% tmin.alter.scope, parameter] = operation_func(subset(parms.1d,tmin %in% tmin.alter.scope)[[parameter]] , row[[parameter]] )
      } else if(parameter %in% names(parms.2d)) {
        parms.2d[parms.2d$tmin %in% tmin.alter.scope, parameter] = operation_func(subset(parms.2d,tmin %in% tmin.alter.scope)[[parameter]] , row[[parameter]] )
      } else {
        stop(paste(parameter,"was not found in any parameter sheet"))# Parameter was not found in any parameter sheet
      }
    }

    sheet_names_for_sweep$parms.1d = parms.1d # altered data.frame goes in sheet_names_for_sweep
    sheet_names_for_sweep$parms.2d = parms.2d # altered data.frame goes in sheet_names_for_sweep

    this.result = SEIR.n.Age.Classes(file_name,sheet_names_for_sweep,also.get.flows=flows.of.interest,post.processing.companion.kit = post.process.kit)

   # # Add on other time series like cumulative incidence
   # this.result$solution = Add.Other.Outcomes(this.result$solution,parms.1d,"third argument currently inoperant")

    # Add on univariate stuff like maxI or maxI.time to outcomes.summary.df
    summary.template = data.frame(etiquette = this.label)
    for(parameter in names(row))
      summary.template[[paste0(parameter, operation.label)]] <- row[[parameter]]

    summary.chunk = this.result$sommaire
    summary.chunk = cbind(summary.template,summary.chunk)
    # summary.chunk$etiquette = this.label # not sure if this is useful to keep
    # for(parameter in names(row))
    #  summary.chunk[[paste0(parameter, operation.label)]] <- row[[parameter]]

    outcomes.summary.df = rbind(outcomes.summary.df,summary.chunk)
    #print(names(outcomes.summary.df))

    # Get ever.been.here info
    if(!is.character(ever.been.here.info$ever.been.here.from.inflows))
      ever.been.here.info = ever.been.here(this.result) # this.result has been altered above but should not cause harm

    #Update list.sweep and df.sweep
    list.sweep[[this.label]] = this.result[c("solution","input.info")] # this.result$solution
    list.sweep.ever.from.inflows [[this.label]] = ever.been.here.info$ever.been.here.from.inflows
    list.sweep.ever.from.outflows[[this.label]] = ever.been.here.info$ever.been.here.from.outflows

    # Add this.label and the parameter multpliers to the this.result$solution data frame
    this.result$solution$etiquette <- this.label
    for(parameter in names(row)) {
      this.result$solution[[paste0(parameter, operation.label)]] <- row[[parameter]]
    }

    df.sweep = rbind(df.sweep,this.result$solution)
  }
  rownames(outcomes.summary.df) = c()

  list(parm.cloud.grid.specs = parm.cloud.grid.specs,
       parms.to.try = parms.to.try.verbose,
       outcomes.summary.df = outcomes.summary.df,
       df.sweep = df.sweep,
       # list.sweep.ever.from.inflows =list.sweep.ever.from.inflows,
       # list.sweep.ever.from.outflows=list.sweep.ever.from.outflows,
       list.sweep=list.sweep)
}

#' Compute correlation between parameters and an outcome variable
#' @export
Assess.covariate.importance = function (don,X,Y,method)
{
  #  if(method == "ANOVA SS type ???")
  #  if(method == "Stepwise")
  #  if(method == "Regression tree")
  #  if(method == "Et cetera...")

  if(method %in% c( "negative-log-p-value","t-test","pearson-partial-correlation-fast") )
  {
    formule = eval(parse(text= paste(Y,"~",paste(X,collapse="+")) ))
    reg = summary(lm(formule,data=don))
    if(method =="negative-log-p-value" )
      result = -log(reg$coefficients[-1,"Pr(>|t|)"] )
    else
    { # "t-test" or "pearson-partial-correlation-fast"
      result = reg$coefficients[-1,"t value"] # result for "t-test" ... "PPCF" needs more work below
      if(method == "pearson-partial-correlation-fast")
        result = result / sqrt(reg$df[2] + result**2) #  "pearson-partial-correlation-fast" aka "PPCF"
    }
  }

  if(grepl("-partial-correlation-slow",method ) )
  {
    cor.method = sub("-partial-correlation-slow","",method)
    result = c()
    for(this.X in X)
    {
      other.covariates = paste(setdiff(X,this.X),collapse="+")
      formule.Y      = eval(parse(text= paste(Y     ,"~",other.covariates) ))
      formule.this.X = eval(parse(text= paste(this.X,"~",other.covariates) ))
      result = c(result,cor(lm(formule.Y,data=don)$residuals,lm(formule.this.X,data=don)$residuals,method=cor.method))
    }
    names(result) = X
  }

  result
}
