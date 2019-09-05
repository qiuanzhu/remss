#' Explained variation measurement.
#' 
#' \code{explain_var_measure} returns explained variation measurement for the grouping scheme.  
#' 
#' @param main_list main_list includes survival indicator variable, Duration time of survival variable and basic group variable.
#' @param stage_list stage_list original of each scheme.
#' @param stage_list_2 stage_list_2 is numerical form of each scheme by using other of stages information.
#' @param covar_list Covariate variables taking into consideration.
#' @param data Data set. 
#' @return Ranking of explained variation measurement and its standardized score.
#' @references Xu, W., et al. 'Refining evaluation methodology on TNM stage system: assessment on HPV-related oropharyngeal cancer.' Austin Biometrics and Biostatistics 2 (2015): 1014.
#' @import  survival
#' @import  stats
explain_var_measure <- function(main_list, stage_list, stage_list_2, covar_list, data) {
    
    n <- length(stage_list_2)
    Scheme <- stage_list
    Score <- rep(0, n)
    Std_Score <- rep(0, n)
    Rank <- c(1:1:n)
    cv_list_m <- paste(covar_list, collapse = "+")
    
    ioparms <- rep(0, 6)
    n1 <- length(covar_list) + 1
    n2 <- dim(data)[1]
    
    tol_list <- c(c(main_list, stage_list_2), covar_list)
    n0 <- length(tol_list)
    col_num1 <- rep(0, n0)
    
    for (i in 1:n0) {
        col_num1[i] <- which(colnames(data) == tol_list[i])
    }
    
    data <- data[order(data[, col_num1[2]]), ]
    
    S_km <- survival::Surv(data[, col_num1[2]], data[, col_num1[1]] != 0)
    km_es <- summary(survfit(S_km ~ 1, type = "kaplan-meier"))
    u_3 <- cbind(unlist(km_es[2], use.names = FALSE), unlist(km_es[6], use.names = FALSE))
    n3 <- dim(u_3)[1]
    
    S_bkm <- survival::Surv(data[, col_num1[2]], data[, col_num1[1]] == 0)
    bkm_es <- summary(survival::survfit(S_bkm ~ 1, type = "kaplan-meier"))
    u_4 <- cbind(unlist(bkm_es[2], use.names = FALSE), unlist(bkm_es[6], use.names = FALSE))
    n4 <- dim(u_4)[1]
    
    value <- rep(0, n2)
    for (i in 1:n2) {
        value[i] <- round(data[i, col_num1[2]], 5) %in% round(u_4[, 1], 5)
    }
    n5 <- min(which(value == 1))
    value2 <- rep(1, n2)
    
    for (i in n5:n2) {
        if (value[i] == 1) {
            k1 <- which(round(data[i, col_num1[2]], 5) == round(u_4[, 1], 5))
            value2[i] <- u_4[k1, 2]
        } else {
            value2[i] <- value2[i - 1]
        }
    }
    value3 <- rep(1, n3)
    
    for (i in 1:n3) {
        k2 <- which(round(u_3[i, 1], 5) == round(data[, col_num1[2]], 5))[1]
        value3[i] <- value2[k2]
    }
    
    ioparms[1] <- n2
    ioparms[2] <- n1
    ioparms[3] <- n3
    
    u <- cbind(u_3, value3)
    
    for (i in 1:n) {
        formu_1 <- paste("Surv", "(", main_list[2], ",", main_list[1], ")", "~", stage_list_2[i], "+")
        cox <- survival::coxph(as.formula(paste(formu_1, cv_list_m, collapse = "+")), data)
        par <- cox$coef
        cov_data <- colMeans(data[, c(col_num1[i + 2], col_num1[(3 + n):n0])])
        
        av_es <- summary(survival::survfit(cox, cov_data, type = "kalbfleisch-prentice"))
        av_es_1 <- cbind(unlist(av_es[2], use.names = FALSE), unlist(av_es[6], use.names = FALSE))
        u0 <- cbind(cbind(u[, 1], av_es_1[, 2]), u[, -1])
        if (n3 < n2) {
            zero_matrix <- matrix(0, n2 - n3, 4)
            u0 <- rbind(u0, zero_matrix)
        }
        
        m1 <- data[, c(col_num1[i + 2], col_num1[(3 + n):n0])]
        m2 <- m1 - matrix(1, dim(m1)[1], dim(m1)[2]) %*% diag(cov_data)
        
        datalines <- cbind(cbind(data[, col_num1[2:1]], m2), u0)
        Score[i] <- explvar(par, datalines, ioparms) * 100
        
    }
    
    table <- data.frame(Scheme, Score)[order(-Score), ]
    diff <- table[n, 2] - table[1, 2]
    
    for (j in 1:n) {
        Std_Score[j] <- (table[j, 2] - table[1, 2])/diff
    }
    
    table <- cbind(table, Std_Score, Rank)
    mylist <- list(`Explained Variance Measurement` = table)
    
    return(mylist)
}



