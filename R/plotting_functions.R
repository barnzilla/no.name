  #' Render a scatter plot
#' @export
#' @importFrom ggplot2 aes coord_flip element_blank element_text geom_bar geom_point geom_text ggplot scale_x_continuous scale_y_continuous theme theme_minimal xlab ylab
#' @importFrom plotly ggplotly
#' @param x A vector
#' @param y A vector
#' @param x_label_text A string
#' @param y_label_text A string
#' @param geom_point_size A numeric
#' @param element_text_size A numeric

get_scatter_plot <- function(x, y, x_label_text = deparse(substitute(x)), y_label_text = deparse(substitute(y)), geom_point_size = 2, element_text_size = 12) {
  if(is.null(x) | is.null(y)) {
    return()
  } else {
    df <- data.frame(x, y)
    ggplotly(
		ggplot(data = df, aes(x, y)) +
		geom_point(color = "#428bca", size = geom_point_size) +
		xlab(x_label_text) +
		ylab(y_label_text) +
		scale_y_continuous(labels = comma) +
		scale_x_continuous(labels = comma) +
		theme_minimal() +
		theme(
			plot.title = element_text(size = element_text_size),
			axis.title.x = element_text(size = element_text_size),
			axis.title.y = element_text(size = element_text_size),
			legend.text = element_blank(),
			legend.title = element_blank(),
			legend.position = "none"
		)
	)
  }
}

#' Render a tornado plot
#' @export
#' @importFrom forcats fct_reorder
#' @importFrom ggplot2 aes coord_flip element_blank element_text geom_bar geom_point geom_text ggplot scale_x_continuous scale_y_continuous theme theme_minimal xlab ylab
#' @importFrom plotly ggplotly
#' @param outcome_variable A string
#' @param parameters The parms.tried.df data frame
#' @param outcomes The outcomes.summary.df data frame
#' @param method A string
#' @param bin_width A numeric
#' @param element_text_size A numeric
#' @param order_by_absolute_value A logical
#' @param add_label A logical

get_tornado_plot <- function(outcome_variable, parameters = parms.tried.df, outcomes = outcomes.summary.df, method = "kendall-partial-correlation-slow", bin_width = 0.5, element_text_size = 12, order_by_absolute_value = FALSE, add_label = FALSE) {
  if(is.null(outcome_variable) | is.null(parameters) | is.null(outcomes) | is.null(method) | is.null(outcome_variable)) {
    return()
  } else {
    what.matters = Assess.covariate.importance(outcomes,names(parameters), outcome_variable, method)
    correlations <- tibble(variable = names(what.matters), coefficient = what.matters)
    correlations$variable <- factor(correlations$variable)
    if(isTRUE(order_by_absolute_value)) {
      correlations$variable <- fct_reorder(correlations$variable, abs(correlations$coefficient), .desc = FALSE)
    } else {
      correlations$variable <- fct_reorder(correlations$variable, correlations$coefficient, .desc = FALSE)
    }
    if(isTRUE(add_label)) {
      label_content <- round(correlations$coefficient, 3)
    } else {
      label_content <- ""
    }
    ggplotly(ggplot(correlations, aes(x = variable, y = coefficient)) +
		geom_bar(color = "#428bca", fill = "#428bca", stat = "identity", width = bin_width, aes()) +
		geom_text(label = label_content, size = 3.5, hjust = -3) +
		coord_flip() +
		theme_minimal() +
		ylab(paste0("Strength of correlation with ", outcome_variable)) +
		theme(
			plot.title = element_text(size = element_text_size),
			axis.title.y = element_blank(),
			axis.title.x = element_text(size = element_text_size),
			legend.text = element_blank(),
			legend.title = element_blank(),
			legend.position = "none"
		),
		tooltip = "text"
    )
  }
}
