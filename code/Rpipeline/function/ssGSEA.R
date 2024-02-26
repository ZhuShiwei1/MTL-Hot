
###  
#' @return Returns a matrix, behavior signature, listed as samples, with the value of ssGSEA score for each signature in each sample.
##################2022-5-6#############

ssGSEA_Score <- function(expr,  	
								sig.list,         
								ssgsea.norm=TRUE,
								z.score=FALSE){   
                                 
                                                                      
	library(GSVA)
	
	# ssGSEA score
	ssGSEA.score <- gsva(expr = expr, 
			gset.idx.list = sig.list, 
			method="ssgsea", 
			ssgsea.norm = ssgsea.norm) 

	if(z.score){
	  ssGSEA.score <- t(apply(ssGSEA.score, 1, scale))
	  colnames(ssGSEA.score) <- colnames(expr)
	}
	
	return(ssGSEA.score)
}

