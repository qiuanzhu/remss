#' Hazard consistency measurement.
#' 
#' \code{hz_cons_measure} returns Hazard consistency for the grouping scheme.  
#' 
#' @param main_list main_list includes survival indicator variable, Duration time of survival variable and basic group variable.
#' @param stage_list stage_list is numerical form of each scheme by using orther of stages information.
#' @param covar_list Covariate variables taking into consideration.
#' @param data Data set. 
#' @return Ranking of hazard consistency measurement and  standardized score.
#' @references Xu, W., et al. 'Refining evaluation methodology on TNM stage system: assessment on HPV-related oropharyngeal cancer.' Austin Biom Biostat 2 (2015): 1014.
#' @import  survival
#' @import  stats
#' @export
hz_cons_measure <- function(main_list, stage_list, covar_list, data) {
    
    n <- length(stage_list)
    Scheme <- substr(stage_list, 1, nchar(stage_list) - 6)
    Score <- rep(0, n)
    Std_Score <- rep(0, n)
    Rank <- c(1:1:n)
    
    
    col_num <- which(colnames(data) == main_list[3])
    data[, col_num] <- as.factor(data[, col_num])
    full_g <- nlevels(data[, col_num])
    
    cv_list_m <- paste(covar_list, collapse = "+")
    formu_b <- paste("Surv", "(", main_list[2], ",", main_list[1], ")", "~", main_list[3], "+")
    Full_L <- -2 * survival::coxph(as.formula(paste(formu_b, cv_list_m, collapse = "+")), data)$loglik[2]
    
    for (i in 1:n) {
        col_num <- which(colnames(data) == stage_list[i])
        data[, col_num] <- as.factor(data[, col_num])
        reduce_g <- nlevels(data[, col_num])
        
        formu_b <- paste("Surv", "(", main_list[2], ",", main_list[1], ")", "~", stage_list[i], "+")
        Reduce_L <- -2 * survival::coxph(as.formula(paste(formu_b, cv_list_m, collapse = "+")), data)$loglik[2]
        Score[i] <- (Reduce_L - Full_L)/(full_g - reduce_g)
    }
    
    table <- data.frame(Scheme, Score)[order(Score), ]
    diff <- table[n, 2] - table[1, 2]
    for (j in 1:n) {
        Std_Score[j] <- (table[j, 2] - table[1, 2])/diff
    }
    
    table <- cbind(table, Std_Score, Rank)
    colnames(table) <- c("Scheme", "Score", "Standardized Score", "Rank")
    mylist <- list(`Hazard Consistency Measurement` = table)
    return(mylist)
    
}

