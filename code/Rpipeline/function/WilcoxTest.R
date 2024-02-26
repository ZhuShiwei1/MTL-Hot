#'  Wilcoxon test
#' 
#' @param  x 
#' @param  y 
#' @param types The names of the x and y vectors
#' @param paired Whether x and y are paired samples
#' @param  alternative Alternate assumption, default is "two-sided"

WilcoxTest <- function(x, y, types = NULL, paired = FALSE, alternative = "two.sided"){
	library(coin)
	if(paired){
		if(all(x == y)){ 
			p.value = 1
			z.statistic = 0
		}else{
			z <- wilcoxsign_test(x ~ y, alternative = alternative) 
			p.value <- pvalue(z) 
			z.statistic <-  statistic(z) 
		}
	}else{

		test.data <- data.frame(value = c(x, y), type = factor(c(rep(types[1], length(x)), rep(types[2], length(y))), levels = types))
		z <- wilcox_test(value ~ type, data = test.data, alternative = alternative)
		p.value <- pvalue(z) 
		z.statistic <-  statistic(z) 
	}
	result <- data.frame(p.value = p.value, z.statistic = z.statistic)
	return(result)
}
