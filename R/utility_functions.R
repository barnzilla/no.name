# This a subset of \\hird1/Saslibe$/SAS_user_defined_macros/Utilities.R (aka S:/SAS_user_defined_macros/Utilities.R)

#' Fetch a package
#' @export
#' @importFrom utils install.packages
#' @param nom_du_package A string
#' @param try.this A string
#' @param verbose A logical

fetch.package =function(nom_du_package,try.this=c("nothing","load","install"),verbose=FALSE)
{
  if(verbose)
    cat(paste("\n try.this =",try.this[1]))
  
  if(try.this[1] =="load")
    try( suppressPackageStartupMessages(library(nom_du_package,character.only=TRUE,quietly=TRUE)) )#,warn.conflicts=FALSE)
  
  if(try.this[1] =="install")    # This part needs more work
  {
    #   browser()
    cat("\ntTrying to install ",nom_du_package)
    #install.packages(nom_du_package)
    install.packages( pkgs = nom_du_package  )   #,repos = "file:////fld6filer/packagerepo-depotprogiciel/miniCRAN"
    found = !grepl("no package called",try(find.package(nom_du_package)) )
    if(found)
    {
      suppressPackageStartupMessages(library(nom_du_package,character.only=TRUE,quietly=TRUE))#,warn.conflicts=FALSE)
      cat("\nInstalled and loaded ",nom_du_package)
    }
    else 
      cat("\nCould not install ",nom_du_package ) # stop(paste("\nCould not install",nom_du_package) )      
  }
  
  # found=require(nom_du_package,character.only=TRUE,quietly=TRUE,warn.conflicts=FALSE)
  found=  paste("package:",nom_du_package,sep="") %in% search()
  
  if(!found && length(try.this)>1)
    found = fetch.package(nom_du_package,try.this[-1])  # recursive call
  
  
  found
}

#' Convert a tibble to a data frame
#' @export
#' @importFrom dplyr is.tbl
#' @param data.frame.ish A tibble
#' @param rename A vector or data frame

from.tbl.to.df  = function(data.frame.ish,rename=NULL) # rename is vector or data.frame with at least $var.names.in
{ 
  if (!fetch.package("dplyr"))
    stop("package dplyr not found")  # needed for dplyr::is.tbl
  
  #browser()
  #BEGIN make rename = data.frame(var.names.in=..., var.names.out=..., colClass=...)
  if(is.null(rename))
    rename = names(data.frame.ish) 
  
  if(!is.data.frame(rename))
    rename = data.frame(var.names.in=rename,var.names.out=rename,stringsAsFactors = FALSE)
  
  if(!any(names(rename)=="var.names.out")) 
    rename$var.names.out=rename$var.names.in
  
  if(!any(names(rename)=="colClass"))
    rename$colClass='unspecified'
  #END   make rename = data.frame(var.names.in=..., var.names.out=...)
  
  this.is.tbl = dplyr::is.tbl(data.frame.ish)
  
  simple.df=data.frame(iota...seq=seq(nrow(data.frame.ish)))
  for(k in 1:nrow(rename))
  {
    if(this.is.tbl)
    {
      if(rename$colClass[k] == "character"  ) simple.df[,rename$var.names.out[k]] = as.character( data.frame.ish[,rename$var.names.in[k]][[1]] )
      if(rename$colClass[k] == "numeric"    ) simple.df[,rename$var.names.out[k]] = as.numeric  ( data.frame.ish[,rename$var.names.in[k]][[1]] ) 
      if(rename$colClass[k] == "unspecified") simple.df[,rename$var.names.out[k]] =             ( data.frame.ish[,rename$var.names.in[k]][[1]] )      
    }
    else
      simple.df[,rename$var.names.out[k]] =  data.frame.ish[,rename$var.names.in[k]]  
  }
  
  simple.df$iota...seq = NULL   
  simple.df
}

#' Overlap length function
#' @export
#' @param L1 A vector
#' @param U1 A vector
#' @param L2 A vector
#' @param U2 A vector

overlap.length = function(L1,U1,L2,U2)
{
  pmax(pmin(U1,U2) - pmax(L1,L2),0)
}


#verbose.save = function(object.name,path.with.trailing.slash="",prefix.suffix=c(prefix="This file contains an R object called ",suffix=".SavedFromR"),time.stamp="")

#' Verbose save function
#' @export
#' @param object.name A string
#' @param path.with.trailing.slash A string
#' @param prefix.suffix A string
#' @param time.stamp A time object

verbose.save = function(object.name,path.with.trailing.slash="",prefix.suffix=c(prefix="This file contains an R object called ",suffix=".SavedFromR"),time.stamp=gsub(":","-",Sys.time()))
{
  if(time.stamp != "")
    time.stamp = paste0(" (", time.stamp,")")
  
  code = paste0(prefix.suffix["prefix"],object.name,time.stamp,prefix.suffix["suffix"])
  code = paste0("save(",object.name,",file='",path.with.trailing.slash,code,"')")
  eval(parse(text=code))
}
