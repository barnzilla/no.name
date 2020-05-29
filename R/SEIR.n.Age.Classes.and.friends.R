#' Main SEIR model function
#' @export
#' @importFrom adaptivetau ssa.adaptivetau
#' @importFrom deSolve lsoda
#' @importFrom dplyr arrange mutate slice
#' @importFrom magrittr %>%
#' @importFrom readxl read_excel
#' @importFrom tidyr pivot_longer
#' @param file.name A string
#' @param sheets.names A string vector
#' @param differential.eqns.func A string
#' @param just.get.functions A logical
#' @param also.get.flows A string
#' @param post.processing.func A string
#' @param post.processing.companion.kit A string
#' @param agegrp.glue A string
#' @param CTMC.random.seeds A vector of integers

SEIR.n.Age.Classes = function(file.name, sheets.names, differential.eqns.func=NULL, just.get.functions=FALSE, also.get.flows=NULL, 
                                       post.processing.func=NULL ,post.processing.companion.kit=NULL, agegrp.glue="", CTMC.random.seeds=NULL) # agegrp.glue=".ag"
{
  expected.sheets.names = c("parms.0d","parms.1d","parms.2d","initial.conditions","model.flow","auxiliary.vars")
  
  if( !all( expected.sheets.names %in% names(sheets.names) ) )
    stop( paste("\n sheets.names must provide the following:",paste(expected.sheets.names,collapse=", ")) )
  
  input.info.verbatim = sheets.names
  input.info.verbatim$file.name = file.name
  input.info.verbatim$agegrp.glue = agegrp.glue
  input.info.verbatim$CTMC.random.seeds = CTMC.random.seeds
  #==========================================================================
  #   Initial conditions and parameters
  #==========================================================================
  
  tmp_parm_value_0d = sheets.names$parms.0d            # 0 dimensional parameters
  tmp_parm_value_1d = sheets.names$parms.1d            # 1 dimensional parameters
  tmp_parm_value_2d = sheets.names$parms.2d            # 2 dimensional parameters
  input_stuff       = sheets.names$initial.conditions  # initial values/conditions
    
  
  if(!is.data.frame(tmp_parm_value_0d))
    tmp_parm_value_0d <- from.tbl.to.df( readxl::read_excel(file.name, sheet = tmp_parm_value_0d) )  # 1 dimensional parameters
  
  if(!is.data.frame(tmp_parm_value_1d))
    tmp_parm_value_1d <- from.tbl.to.df( readxl::read_excel(file.name, sheet = tmp_parm_value_1d) )  # 1 dimensional parameters
  
  if(!is.data.frame(tmp_parm_value_2d))  
    tmp_parm_value_2d <- from.tbl.to.df( readxl::read_excel(file.name, sheet = tmp_parm_value_2d) ) # 2 dimensional parameters
  
  if(!is.data.frame(input_stuff)) 
    input_stuff  <- from.tbl.to.df( readxl::read_excel(file.name, sheet = input_stuff) ) # initial values/conditions
  
  raw.init.conditions = input_stuff # keep snapshot for output.  input_stuff may be altered later
  
  nagegrp <- length(unique(tmp_parm_value_1d$agegrp)) # number of age groups
  nagegrp = ncol(raw.init.conditions) - 1             # number of age groups
  
  #BEGIN get raw.compartments.age and pretty.compartments.age
  agegrp.suffix.pretty = agegrp.suffix = ""
  if(nagegrp > 1)
  {
    agegrp.suffix        = 1:nagegrp
    agegrp.suffix.pretty = paste0(agegrp.glue,1:nagegrp) # if agegrp.glue = ".ag" this generates ".ag1" ".ag2" ".ag3" 
  }
    
     raw.compartments =   raw.init.conditions$NAME # e.g.  S D L
     raw.compartments.age = c ( t( outer( raw.compartments,agegrp.suffix       ,paste0) ) ) # e.g. S1    S2    S3    S4    S5    D1    D2    ...
  pretty.compartments.age = c ( t( outer( raw.compartments,agegrp.suffix.pretty,paste0) ) ) # e.g. S.ag1 S.ag2 S.ag3 S.ag4 S.ag5 D.ag1 D.ag2 ...
  
  if(length(unique(raw.compartments.age)) != length(raw.compartments.age)  )
    stop("\n Ambuiguity from compartment names and ages")
  # Ambiguity could arise if there are 2 compartments named S and S1 (say) and there are > 10 age groups
  # Then S11 could be 11th age group of S or first age group of S1
  #END get raw.compartments.age and pretty.compartments.age

  nrow_   <- dim(tmp_parm_value_1d)[1]/nagegrp
  
  parm_value_0d <- dplyr::arrange(tmp_parm_value_0d, tmin        ) # sort by tmin 
  parm_value_0d$isim = 1:nrow_  # (1:nrow(parm_value_0d))
  
  parm_value_1d <- dplyr::arrange(tmp_parm_value_1d, tmin, agegrp) # sort by tmin agegrp
  parm_value_1d <- parm_value_1d %>%
    mutate(isim = rep(1:nrow_, each=nagegrp)) 
  
  parm_value_2d <- dplyr::arrange(tmp_parm_value_2d, tmin, cagegrp, ragegrp) # sort by tmin cagegrp  ragegrp
  parm_value_2d <- parm_value_2d %>%
    mutate(isim = rep(1:nrow_, each = nagegrp*nagegrp))
 
  #===================================================================
  # Build function eval.post.processing.func (if not provided) 
  #===================================================================
 
  code.body.df = sheets.names$post.processing # data.frame or string
  if(!is.data.frame(code.body.df))
    code.body.df = from.tbl.to.df( readxl::read_excel(file.name, sheet = code.body.df ) )
 
  code.body.df$code[is.na(code.body.df$code)] = ""
  
  eval.post.processing.func = post.processing.func
  if(is.null(eval.post.processing.func))
  {
    code.body.char = code.body.df$code
    code.head.char = c("function(list.solution.etcetera){",
                       "# list.solution.etcetera contains $solution and possibly $solution.inflows and/or $solution.outflows" , 
                       "# Probably could use    with(list.solution.etcetera)   instead of copying" , 
                       "solution            = list.solution.etcetera$solution",
                       "solution.inflows    = list.solution.etcetera$solution.inflows",
                       "solution.outflows   = list.solution.etcetera$solution.outflows",
                       "CTMC.df             = list.solution.etcetera$CTMC.out$CTMC.df",
                       "parms.0d            = list.solution.etcetera$input.info$parms.0d",
                       "parms.1d            = list.solution.etcetera$input.info$parms.1d",
                       "parms.2d            = list.solution.etcetera$input.info$parms.2d",
                       "initial.conditions  = list.solution.etcetera$input.info$initial.conditions",
                       "companion.kit       = list.solution.etcetera$functions$post.processing.companion.kit",
                       "input.info.verbatim = list.solution.etcetera$input.info.verbatim",
                       "sommaire            = data.frame(osqvhfipumzcdjblkwrgnexaty=4)" ,
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
    
    if(!is.data.frame(auxiliary.vars))
      auxiliary.vars   <- from.tbl.to.df( readxl::read_excel(file.name, sheet = auxiliary.vars ) )  # auxiliary variables in equations
    
    auxiliary.vars = subset(auxiliary.vars, !is.na(code)   )  # auxiliary.vars$code will be top portion of with(...,{ ... })  part
   
    
    if(!is.data.frame(model_flows_tmp))
      model_flows_tmp  <- from.tbl.to.df( readxl::read_excel(file.name, sheet = model_flows_tmp) )  # arrows in flowchart
    
    model_flows_tmp = subset( model_flows_tmp,!is.na(expression) & activation ==1)
    model_flows = model_flows_tmp     # keep model_flows_tmp for output. Work on model_flows from this point on

 #   transitions.for.adaptivetau.early = paste("list(" , paste(paste('c(',model_flows$From,'=-1,',model_flows$To,'=+1)'),collapse=' , ') , ')')
 #   transitions.for.adaptivetau.early = eval(parse(text=transitions.for.adaptivetau.early))
 #   names(transitions.for.adaptivetau.early) = paste(model_flows$From,'-->',model_flows$To)
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

    # BEGIN Build character vector differential.eqns.char  (and ratefunc.for.adaptivetau.chunks + transitions.for.adaptivetau )
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
    ratefunc.for.adaptivetau.chunks = c()
    transitions.for.adaptivetau = list()
    
    for(k in seq(nrow(model_flows)))
    {
      this.expression = model_flows$expression[k]
      outflows[model_flows$From[k]] = paste(outflows[model_flows$From[k]],"+", this.expression)
      inflows [model_flows$To[k]  ] = paste( inflows[model_flows$To[k]  ],"+", this.expression)
      
      ratefunc.for.adaptivetau.chunks = c(ratefunc.for.adaptivetau.chunks,this.expression)
      transitions.for.adaptivetau.chunk = c(-1,+1)
      names(transitions.for.adaptivetau.chunk) = model_flows[k,c("From","To")]
      transitions.for.adaptivetau[[k]] = transitions.for.adaptivetau.chunk
    }
    names(transitions.for.adaptivetau) = paste(model_flows$From,'-->',model_flows$To)
    names(ratefunc.for.adaptivetau.chunks) = names(transitions.for.adaptivetau)
    
    differential.eqns.char = paste("(flow.multiplier$inflow)*(",inflows,") + (flow.multiplier$outflow)*(",outflows,")")
    differential.eqns.char = gsub(  "#####$$$$$$$$$$$!@@@@@@@@ +","" ,differential.eqns.char ,fixed=TRUE)
    differential.eqns.char = gsub("( #####$$$$$$$$$$$!@@@@@@@@ )","0",differential.eqns.char ,fixed=TRUE)
    differential.eqns.char = paste0("derivative.",boxes,"=",differential.eqns.char)
    
    #if( any(grepl("#####$$$$$$$$$$$!@@@@@@@",differential.eqns.char)) )  browser()

    
    rm(model_flows,boxes,inflows,outflows,this.expression) # clean up intermediate objects no longer required
    
    # END Build character vector differential.eqns.char  (and ratefunc.for.adaptivetau.chunks + transitions.for.adaptivetau )
    

    # BEGIN Build function eval.differential.eqns.func
    
     ODE.func.header.char = c("function(list.boxes, list.parms, time.now, flow.multiplier=list(inflow=1,outflow=-1)){" ,"#browser() #not a good place for it")
    CTMC.func.header.char = c("function(list.boxes, list.parms, time.now) {"                                           ,"#browser() #not a good place for it")  # , "print(time.now)"
    
    # BEGIN get ODE.within.with.part.char and CTMC.within.with.part.char
    compartments =   input_stuff$NAME # e.g.  S D L   (note it is not be the same as raw.compartments if also.get.flows != NULL)
    compartments.age = c ( t( outer( compartments,seq(nagegrp),paste0) ) ) # e.g. S1 S2 S3 S4 S5 D1 D2 ...
    out.char      = paste0("out=c(", paste0("derivative.", paste(compartments    , collapse=", derivative.")),")")
    names.out.char= paste0("out=c(", paste0("derivative.", paste(compartments.age, collapse=",derivative." )),")") 
   #names.out.char = gsub(",","','"  ,names.out.char ,fixed=TRUE) # replace , by ','     may generate a line that is too long causing problems
    names.out.char = gsub(",","',\n'",names.out.char ,fixed=TRUE) # replace , by ',\n'   invoids the problematic single line that is too long
    names.out.char = gsub("(","('"   ,names.out.char ,fixed=TRUE) # replace ( by ('
    names.out.char = gsub(")","')"   ,names.out.char ,fixed=TRUE) # replace ) by ')
    names.out.char = gsub("out=","names(out)=",names.out.char ,fixed=TRUE)
    
     ODE.within.with.part.char = paste(c(differential.eqns.char,"#browser()",out.char,names.out.char,"#print(out)","list(out)"),collapse="\n") # First tried ";\n" but realised that "\n" is sufficient
    
    CTMC.within.with.part.char = paste("c(", paste(ratefunc.for.adaptivetau.chunks, collapse=" ,\n "),")")
    # END get ODE.within.with.part.char and CTMC.within.with.part.char
    
    ODE.with.char = c("relevant.parm=sum(time.now >= list.parms$tmin.vect)",
                       "names(time.now     )='time.now'",
                       "# list.parms$parms.1d = list.parms$parms.1d[[relevant.parm]] # for CTMC.with.char    # DELETE",
                       "# list.parms$parms.2d = list.parms$parms.2d[[relevant.parm]] # for CTMC.with.char    # DELETE",
                       "## list.parms = as.list(c(list.parms['tmin.vect'], list.parms$parms.1d[[relevant.parm]] , list.parms$parms.2d[[relevant.parm]])) # for CTMC.with.char",
                       "list.parms=as.list(c( time.now, list.parms ))" ,
                       "with(as.list(c(list.boxes, list.parms)),{",                                 
                       "#browser() # a good place for it" )
   # ODE.with.char = gsub("list.parms$parms.1d[[relevant.parm]] , list.parms$parms.2d[[relevant.parm]]","list.parms",CTMC.with.char,fixed=TRUE)
    CTMC.with.char = gsub("##","",ODE.with.char,fixed=TRUE)
  
    
    # BEGIN Create functions eval.differential.eqns.func and ratefunc.for.adaptivetau
    eval.differential.eqns.func = paste(c( ODE.func.header.char,  ODE.with.char, auxiliary.vars$code,  ODE.within.with.part.char,"})" ,"}"),collapse="\n" ) # code of eval.differential.eqns.func
    ratefunc.for.adaptivetau    = paste(c(CTMC.func.header.char, CTMC.with.char, auxiliary.vars$code, CTMC.within.with.part.char,"})" ,"}"),collapse="\n" ) # code of ratefunc.for.adaptivetau
    
    # cat(eval.differential.eqns.func) # quick look at how it look like
    
    eval.differential.eqns.func = eval(parse(text = eval.differential.eqns.func ))
    ratefunc.for.adaptivetau    = eval(parse(text = ratefunc.for.adaptivetau    ))  
    # END Create functions eval.differential.eqns.func and ratefunc.for.adaptivetau
    
    
    rm(CTMC.within.with.part.char, CTMC.func.header.char, CTMC.with.char,
        ODE.within.with.part.char,  ODE.func.header.char,  ODE.with.char, out.char) # clean up 
    # END Build function eval.differential.eqns.func
  }
 
  resultat = list(functions=list( differential.eqns.func    = eval.differential.eqns.func , 
                                  CTMC.eqns.func            = ratefunc.for.adaptivetau,
                                  post.processing.func      = eval.post.processing.func,
                                  post.processing.companion.kit = post.processing.companion.kit )) 
  
  if(just.get.functions)
    return ( resultat$functions )
  
  #==========================================================================
  #  Main routine
  #==========================================================================
  # The SEIR model with N age classes
  #
  SEIR.n.Age.Classes.within.loop <- function( time=NULL, tmin.vect=NULL, ageless.parms = NULL, age.parms = NULL, age.age.parms = NULL,list.inits = NULL, not.parms=  c("tmin", "tmax", "agegrp", "cagegrp", "ragegrp", "isim"))
  {
    if (is.null(ageless.parms))
      stop("undefined 'ageless.parms'")
    
    if (is.null(age.parms))
      stop("undefined 'age.parms'")    
    
    if (is.null(age.age.parms))
      stop("undefined 'age.age.parms'")
    
    if (is.null(time))
      stop("undefined 'time'")
    
    if (is.null(list.inits))
      stop("undefined 'list.inits'")
    
    if (is.null(tmin.vect))
      stop("undefined 'tmin.vect'")
    
    nage = nrow(age.parms)
    
    list.parms.for.lsoda <- list(tmin.vect=tmin.vect)
    
    for(k in setdiff(names(ageless.parms), not.parms ))
      list.parms.for.lsoda[[k]] <- ageless.parms[,k]
    
    for(k in setdiff(names(age.parms), not.parms ))
      list.parms.for.lsoda[[k]] <- age.parms[,k]
    
    for(k in setdiff(names(age.age.parms), not.parms ))
    {
      temp<- array(NA, c(nage,nage))
      temp[cbind(age.age.parms$ragegrp, age.age.parms$cagegrp)] <- age.age.parms[,k]
      list.parms.for.lsoda[[k]] <- temp
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
      
      eval.differential.eqns.func(list.inits, list.parms, time)
      
      
    } #end of function calculate_derivatives
    
    ###---------------------------------------------------------------------------------
    ### Solver for Ordinary Differential Equations (ODE), Switching Automatically
    ### Between Stiff and Non-stiff Methods
    ###---------------------------------------------------------------------------------
    #browser()
    list.parms  = list.parms.for.lsoda
    output <-  lsoda(y = unlist(list.inits),
                     times = time,
                     func =  calculate_derivatives,
                     parms = list.parms,
                     names.inits = names(list.inits))
    
    return(output)
  }  # END  of function SEIR.n.Age.Classes.within.loop
  
  lean.ssa.adaptivetau <- function(sim_number, init.conditions, params,racines,etiquettes) 
  {
    if(!is.na(racines[sim_number]))
      set.seed(racines[sim_number])
    
    #browser()
   #simu.all = ssa.exact ( init.values = init.conditions,  # need to remove argument tl.params
    simu.all = ssa.adaptivetau(init.values = init.conditions, 
                           transitions = transitions.for.adaptivetau, 
                           rateFunc = ratefunc.for.adaptivetau, 
                           params = params, 
                           tf = tmax,
                           tl.params = list(epsilon=0.01))
    
    temps = sort(c(simu.all[,"time"],seq(tmax-1)))
    int.whereabouts  = which(temps %in% 0:tmax)   # int.whereabouts  such that temps   [int.whereabouts] = 0:tmax
    int.whereabouts2 = int.whereabouts - (0:tmax) # int.whereabouts2 such that simu.all[int.whereabouts2,"time"][k+1] is largest time <= k
    simu.lean = as.data.frame(simu.all[int.whereabouts2,]) # snapshots at t= 0, 1, 2, ... , tmax
    simu.lean$time = 0:tmax
    simu.lean$what.is.it = etiquettes[sim_number] 
    
    list(simu.all=simu.all,simu.lean=simu.lean)
  }
  ################################################################
  #        To run the example
  ################################################################
  
  
  excluded_names <- c("tmin", "tmax","agegrp","cagegrp","ragegrp","isim")
  
  nTimeSegments <- max(parm_value_1d$isim)
  previous.tmax <- 0
  df.out = c()
  CTMC.parms.info = list(tmin.vect=sort(unique(parm_value_1d$tmin)),parms.0d=list(),parms.1d=list(),parms.2d=list())
  updated_init_list = init_list
  
  for(segment in seq(1, nTimeSegments, 1))
  {

    this.parameter.by.nothing <- subset(parm_value_0d, isim == segment)
    this.parameter.by.age     <- subset(parm_value_1d, isim == segment)
    this.parameter.by.age.age <- subset(parm_value_2d, isim == segment)
    CTMC.parms.info$parms.0d[[segment]] = this.parameter.by.nothing
    CTMC.parms.info$parms.1d[[segment]] = this.parameter.by.age
    CTMC.parms.info$parms.2d[[segment]] = this.parameter.by.age.age
    
    tmin <- unique(c(this.parameter.by.nothing$tmin, this.parameter.by.age$tmin, this.parameter.by.age.age$tmin))
    tmax <- unique(c(this.parameter.by.nothing$tmax, this.parameter.by.age$tmax, this.parameter.by.age.age$tmax))
    
    
    if(length(tmin)>1 || length(tmax)>1 || tmin>=tmax )
      stop(paste0("Unexpected pattern in tmin, tmax for interval ", segment))
    
  # tt <- seq(0   , tmax - tmin, by = 1)
    tt <- seq(tmin, tmax       , by = 1)
    
    if(tmin != previous.tmax)
      stop(paste(interval.label , "\n  Interval lower bound not equal to previous interval upper bound"))
    
 
    
    previous.tmax <- tmax
    out <- SEIR.n.Age.Classes.within.loop( time=tt,
                                           tmin.vect = CTMC.parms.info$tmin.vect,
                                           ageless.parms = this.parameter.by.nothing,
                                           age.parms     = this.parameter.by.age,
                                           age.age.parms = this.parameter.by.age.age,
                                           list.inits = updated_init_list)
    
    out <- as.data.frame(out)
    # ode/lsoda Output diagnostic #######################
    #diagn <- diagnostics.deSolve(out)
    
    out$time <- seq(tmin,tmax,1) 
    out_for_init <- out %>%
      slice(nrow(out)) %>%    # select last row
      pivot_longer(-time)     # fat to skinny
    init <- out_for_init$value
    names(init) <- out_for_init$name     
    
 #   rowns <- names(select(out,-c(time)))
 #   out <- out %>%
 #     mutate(N_tot = rowSums(.[rowns]))  # Total number of individuals 
    
    
    #updating the initial values  
    for(k in 1:length(updated_init_list)){
      updated_init_list[[k]][1:nagegrp] <- init[seq(nagegrp*(k-1)+1,nagegrp*k)] 
    }
    
   if(segment < nTimeSegments)
      out = out[-nrow(out),]
    # Add outputs to the list
    df.out = rbind(df.out,out)

  } #end for(segment in seq(1, nTimeSegments, 1))
  

  CTMC.out=NULL
  if(nagegrp == 1 && !is.null(CTMC.random.seeds))
  {
    iota.CTMC = seq(length(CTMC.random.seeds))
    CTMC.etiquettes = paste("CTMC sim",iota.CTMC)
    CTMC.list = lapply(iota.CTMC, FUN = lean.ssa.adaptivetau, unlist(init_list), CTMC.parms.info, CTMC.random.seeds, CTMC.etiquettes)
    names(CTMC.list) = CTMC.etiquettes 
    CTMC.df = c()
    for(k in CTMC.etiquettes)
    {
      CTMC.df = rbind(CTMC.df,CTMC.list[[k]]$simu.lean)
      CTMC.list[[k]] = CTMC.list[[k]]$simu.all
    }
    CTMC.out=list(CTMC.list=CTMC.list,CTMC.df=CTMC.df)
  }
  
  #browser()
  if(!is.null(CTMC.out))
  {
    resultat$CTMC.out = CTMC.out
   #resultat$functions$CTMC.eqns.func = ratefunc.for.adaptivetau # always provide or not ?
  }
  
  #BEGIN add on $solution.inflows, $solution.outflows and $solution to resultat
  solution = df.out
  if(any(also.get.flows == "inflows"))
  {
   # var.names = c ( t( outer( paste0("inflows.",raw.compartments),agegrp.suffix,paste0) ) )
    var.names = paste0("inflows.",raw.compartments.age)
    resultat$solution.inflows  = solution[,c("time",var.names)]
    colnames( resultat$solution.inflows) = c("time",pretty.compartments.age)
  }
  if(any(also.get.flows == "outflows"))
  {
    #var.names = c ( t( outer( paste0("outflows.",raw.compartments),agegrp.suffix,paste0) ) )
    var.names = paste0("outflows.",raw.compartments.age)
    resultat$solution.outflows  = solution[,c("time",var.names)]
    resultat$solution.outflows[,-1] = - resultat$solution.outflows[,-1] 
    colnames( resultat$solution.outflows) = c("time",pretty.compartments.age)
  }  
  resultat$solution = solution[,c("time",raw.compartments.age)]  
  colnames( resultat$solution) = c("time",pretty.compartments.age)
  #END add on $solution.inflows, $solution.outflows and $solution to resultat
  
  resultat$input.info = list(parms.0d=tmp_parm_value_0d, parms.1d=tmp_parm_value_1d, parms.2d=tmp_parm_value_2d, initial.conditions=raw.init.conditions,
                             auxiliary.vars=auxiliary.vars, model.flow=model_flows_tmp, post.processing=code.body.df)
  resultat$input.info.verbatim = input.info.verbatim
  
  eval.post.processing.func(resultat)
} #end of SEIR.n.Age.Classes function

#' Get model chunk
#' @export
#' @param model_flows A string
#' @param init.cond A string
#' @param init.cond.numeric.vars A vector string
#' @param which.flow A string
#' @param NewFrom A string
#' @param NewTo A string

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
  n.age.grp = length(init.cond.numeric.vars)
  zeros.char = paste0("c(",paste(rep("0",n.age.grp),collapse=","),")")
  model.inflow.orphans = model.outflow.orphans = c()
  if(length(inflow.orphans) > 0)
  {
    model.inflow.orphans = silly_model[rep(1,length(inflow.orphans)),]
    model.inflow.orphans$To = paste0(NewTo,inflow.orphans)
    model.inflow.orphans$expression = zeros.char
  }
  if(length(outflow.orphans) > 0)
  {
    model.outflow.orphans = silly_model[rep(1,length(outflow.orphans)),]
    model.outflow.orphans$From = paste0(NewFrom,outflow.orphans)
    model.outflow.orphans$expression = zeros.char
  }
  #END handle the orphans
  
  list(init=silly_init , model = rbind(silly_model,model.inflow.orphans,model.outflow.orphans) )
}

#' Put model parameters into a list
#' @export
#' @param df.parms.1d A data frame
#' @param df.parms.2d A data frame
#' @param these.are.not.parms A data frame

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

#' Another SEIR model function
#' @export
#' @param SEIR.object A SEIR model object

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

#' Save SEIR model as Excel file
#' @export
#' @importFrom openxlsx write.xlsx
#' @param input.info.list A data frame
#' @param file_name A string
#' @param map.names A vector string

SaveModelInExcel = function(input.info.list,file_name,map.names)
{  # map.names is named vector saying map.names["parms.1d"  ] =  "Parameters by Age"
 # map.names = unlist(sheet_names) # allows to map "parms.1d" into "Parameters by Age" for instance
  names(input.info.list) = map.names[names(input.info.list)] # $parms.1d is now $'Parameters by Age'
  openxlsx::write.xlsx(input.info.list, file_name, colWidths = c(NA, "auto", "auto"))

}

#' Perform parameter sweeping
#' @export
#' @importFrom lhs randomLHS
#' @importFrom stats qunif runif
#' @importFrom triangle qtriangle
#' @param SEIR.object A SEIR model object
#' @param parm.cloud.grid.specs A data frame
#' @param only.show.parms.to.try A logical

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

  baseline.parms.0d = SEIR.object$input.info$parms.0d  # data frame of 0d parameters (i.e. not by age) 
  baseline.parms.1d = SEIR.object$input.info$parms.1d  # data frame of 1d parameters (i.e. by age)
  baseline.parms.2d = SEIR.object$input.info$parms.2d  # data frame of 2d parameters (i.e. by age and age)
  
  post.process.kit  = SEIR.object$functions$post.processing.companion.kit
  agegrp.glue       = SEIR.object$input.info.verbatim$agegrp.glue
  CTMC.racines.alea = SEIR.object$input.info.verbatim$CTMC.random.seeds
  
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
    
    parms.0d = baseline.parms.0d  # data frame of 0d parameters to be altered
    parms.1d = baseline.parms.1d  # data frame of 1d parameters to be altered
    parms.2d = baseline.parms.2d  # data frame of 2d parameters to be altered
    
    # Modify the parameter values at the specified tmin.alter.scope
    for(parameter in names(row)) {
      if(parameter %in% names(parms.0d)) 
        parms.0d[parms.0d$tmin %in% tmin.alter.scope, parameter] = operation_func(subset(parms.0d,tmin %in% tmin.alter.scope)[[parameter]] , row[[parameter]] )
      else if(parameter %in% names(parms.1d)) 
        parms.1d[parms.1d$tmin %in% tmin.alter.scope, parameter] = operation_func(subset(parms.1d,tmin %in% tmin.alter.scope)[[parameter]] , row[[parameter]] )
      else if(parameter %in% names(parms.2d)) 
        parms.2d[parms.2d$tmin %in% tmin.alter.scope, parameter] = operation_func(subset(parms.2d,tmin %in% tmin.alter.scope)[[parameter]] , row[[parameter]] )
      else 
        stop(paste(parameter,"was not found in any parameter sheet"))# Parameter was not found in any parameter sheet
      
    }
    
    sheet_names_for_sweep$parms.0d = parms.0d # altered data.frame goes in sheet_names_for_sweep
    sheet_names_for_sweep$parms.1d = parms.1d # altered data.frame goes in sheet_names_for_sweep
    sheet_names_for_sweep$parms.2d = parms.2d # altered data.frame goes in sheet_names_for_sweep
    
    this.result = SEIR.n.Age.Classes(file_name,sheet_names_for_sweep,also.get.flows=flows.of.interest,post.processing.companion.kit = post.process.kit, agegrp.glue=agegrp.glue, CTMC.random.seeds=CTMC.racines.alea) 
    
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
  outcomes.summary.df$etiquette = c() # drop etiquette.  Not really useful to keep.
  
  list(parm.cloud.grid.specs = parm.cloud.grid.specs,
       parms.to.try = parms.to.try.verbose, 
       outcomes.summary.df = outcomes.summary.df,
       df.sweep = df.sweep,
       # list.sweep.ever.from.inflows =list.sweep.ever.from.inflows,
       # list.sweep.ever.from.outflows=list.sweep.ever.from.outflows,
       list.sweep=list.sweep) 
}

#' Assess covariate importance
#' @export
#' @importFrom stats cor lm time
#' @param don A data frame
#' @param X A numerical vector
#' @param Y A numerical vector
#' @param method A string

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
