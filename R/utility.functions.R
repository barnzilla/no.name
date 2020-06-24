#' Return the complete path to a demo model Excel workbook
#'
#' @export
#'
#' @param x a character element representing a file name ("demo.1.xls", "demo.3.xls", "demo.4.xls").
#'
#' @examples
#' path.to.demo.model <- get.path("demo.1.xls")
#' 
#' @return a character element representing the complete path to the file. If the file cannot be found, an empty element will be returned.

get.path <- function(x) {
  system.file("extdata", x, package = "sweepr")
}