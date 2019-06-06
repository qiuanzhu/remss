#' Hazard discrimination measurement.
#' 
#' \code{hz_dis_measure} returns hazard discrimination for the grouping scheme.  
#' 
#' @param main_list main_list includes survival indicator variable, Duration time of survival variable and basic group variable.
#' @param stage_list stage_list original of each scheme.
#' @param stage_list2 stage_list2 is numerical form of each scheme by using orther of stages information.
#' @param covar_list Covariate variables taking into consideration.
#' @param data Data set. 
#' @return Ranking of hazard discrimination measurement and its standardized score.
#' @references Xu, W., et al. 'Refining evaluation methodology on TNM stage system: assessment on HPV-related oropharyngeal cancer.' Austin Biom Biostat 2 (2015): 1014.
#' @import  survival
#' @import  stats
#' @export
hz_dis_measure <- function(main_list, stage_list, stage_list2, covar_list, data) {
    
    n <- length(stage_list)
    Scheme <- stage_list
    Score <- rep(0, n)
    Std_Score <- rep(0, n)
    Rank <- c(1:1:n)
    cv_list_m <- paste(covar_list, collapse = "+")
    
    for (i in 1:n) {
        col_num <- which(colnames(data) == stage_list[i])
        data[, col_num] <- as.factor(data[, col_num])
        reduce_g <- nlevels(data[, col_num])
        
        formu_1 <- paste("Surv", "(", main_list[2], ",", main_list[1], ")", "~", stage_list[i], "+")
        Reduce_L <- -2 * survival::coxph(as.formula(paste(formu_1, cv_list_m, collapse = "+")), data)$loglik[2]
        
        formu_2 <- paste("Surv", "(", main_list[2], ",", main_list[1], ")", "~", stage_list2[i], "+")
        Reduce_L2 <- -2 * survival::coxph(as.formula(paste(formu_2, cv_list_m, collapse = "+")), data)$loglik[2]
        
        Score[i] <- (Reduce_L2 - Reduce_L)/(reduce_g - 2)
    }
    
    table <- data.frame(Scheme, Score)[order(Score), ]
    diff <- table[n, 2] - table[1, 2]
    for (j in 1:n) {
        Std_Score[j] <- (table[j, 2] - table[1, 2])/diff
    }
    
    table <- cbind(table, Std_Score, Rank)
    colnames(table) <- c("Scheme", "Score", "Standardized Score", "Rank")
    mylist <- list(`Hazard Discrimination Measurement` = table)
    return(mylist)
    
}
