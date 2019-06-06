#' Refining evaluation methodology on stage system.
#' 
#' \code{rank} returns five measurements for the grouping scheme and its overall rank.  
#' 
#' @param os Survival indicator, 1 for death, 0 for censoring.
#' @param ostime Duration time of survival.
#' @param groupvar Basic group variable having the most number of stages.
#' @param scheme Different grouping scheme, which has less stages than the basic group variable.
#' @param order The orther of stages in each grouping, from  
#' @param covariate Covariate variables taking into consideration. 
#' @param weight Weight on five measurements of grouping scheme. 
#' @param data Data set. 
#' @return Ranking of five measurements, which are Hazard consistency, Hazard discrimination, Explained variation, Likelihood difference and Balance. By standardized each measurement score, we provides overall ranking of schemes.
#' @examples  
#' data(Rdata)
#' Scheme=c('Scheme.1','Scheme.2','Scheme.3')
#' Covar=c('Age','Treatment')
#' weight=c(1,1,0.5,0.5,1)
#' Order=list(c('I','II','III'),c('I','II','III','IV'),c('I','II','III','IV'))
#' rank(os='OS',ostime='survmonth',groupvar='Basic_group', scheme=Scheme, order=Order,
#'      covariate=Covar,weight=weight,data=Rdata)
#' @references Xu, W., et al. 'Refining evaluation methodology on TNM stage system: assessment on HPV-related oropharyngeal cancer.' Austin Biom Biostat 2 (2015): 1014.
#' @export
#' @import  survival
#' @import  stats
rank <- function(os, ostime, groupvar, scheme, order, covariate, weight, data) {
    
    total_var <- c(os, ostime, groupvar, scheme, covariate)
    total_n <- length(total_var)
    data_red <- data[, total_var]
    
    var_col <- rep(0, total_n)
    Num_ind <- rep(0, total_n)
    
    for (i in 1:total_n) {
        var_col[i] <- which(colnames(data_red) == total_var[i])
        Num_ind[i] <- is.numeric(data_red[, var_col[i]])
    }
    
    
    ob_n_raw <- dim(data_red)[1]
    data_red <- data_red[!(rowSums(is.na(data_red) | data_red == "")), ]
    ob_n_com <- dim(data_red)[1]
    
    
    if (length(order) != length(scheme)) {
        stop("All stage scheme order needed.")
    }
    
    if (length(weight) != 5 || sum(weight < 0) >= 1) {
        stop("Positive Weight for 5 measurements needed.")
    }
    
    
    var_name <- total_var
    num_var <- Num_ind
    Mean <- rep(list(c()), total_n)
    Sd <- rep(list(c()), total_n)
    Level <- list()
    n_levels <- rep(0, total_n)
    
    
    for (i in 1:total_n) {
        if (Num_ind[i] == 1) {
            Mean[i] <- mean(data_red[, var_col[i]])
            Sd[i] <- sd(data_red[, var_col[i]])
            Level[[length(Level) + 1]] <- "Null"
        } else {
            Mean[i] <- NA
            Sd[i] <- NA
            n_levels[i] <- length(levels(data_red[, var_col[i]]))
            if (n_levels[i] >= 3) {
                Level[[length(Level) + 1]] <- c(levels(data_red[, var_col[i]])[1:2], "...")
            } else {
                Level[[length(Level) + 1]] <- levels(data_red[, var_col[i]])
            }
        }
    }
    
    
    infor_t <- data.frame(cbind(var_name, num_var, Mean, Sd, Level))
    colnames(infor_t) <- c("Variable Name", "Numeric", "Mean", "Standard Deviation", "Level Set")
    
    mis_num <- ob_n_raw - ob_n_com
    infor_t2 <- data.frame(cbind(ob_n_raw, ob_n_com, mis_num))
    colnames(infor_t2) <- c("Total obs", "Complete obs", "Imcomplete obs")
    
    T0 <- list(`Variable Information` = infor_t, Observation = infor_t2)
    
    
    new_var_list <- c()
    for (i in 1:length(scheme)) {
        var_name <- paste(scheme[i], "numer", sep = "_")
        k <- 3 + i
        new_col <- rep(0, ob_n_com)
        if (Num_ind[k] == 1) {
            
            new_col <- data_red[, var_col[k]]
            
        } else {
            
            for (j in 1:ob_n_com) {
                new_col[j] <- match(data_red[, var_col[k]][j], unlist(order[i]))
            }
        }
        data_red <- cbind(data_red, new_col)
        colnames(data_red)[total_n + i] <- var_name
        new_var_list <- c(new_var_list, var_name)
    }
    
    
    new_var_list2 <- c()
    num_k <- 0
    for (i in 1:length(covariate)) {
        k <- 3 + length(scheme) + i
        if (n_levels[k] == 0) {
            new_col <- data_red[, var_col[k]]
            var_name <- paste(covariate[i], "numer", sep = "_")
        } else if (n_levels[k] == 1) {
            new_col <- rep(1, ob_n_com)
            var_name <- paste(covariate[i], "numer", sep = "_")
        } else {
            new_col <- matrix(0, ob_n_com, n_levels[k] - 1)
            var_name <- paste(covariate[i], "numer", 1:(n_levels[k] - 1), sep = "_")
            for (j in 1:ob_n_com) {
                mk1 = match(data_red[, var_col[k]][j], levels(data_red[, var_col[k]])[-1])
                if (is.na(mk1)) {
                  new_col[j, ] = rep(0, n_levels[k] - 1)
                } else {
                  new_col[j, ] = c(rep(0, mk1 - 1), 1, rep(0, n_levels[k] - 1 - mk1))
                }
            }
        }
        data_red <- cbind(data_red, new_col)
        num_i <- length(var_name)
        colnames(data_red)[(total_n + length(scheme) + num_k + 1):(total_n + length(scheme) + num_k + num_i)] <- var_name
        num_k <- num_k + num_i
        new_var_list2 <- c(new_var_list2, var_name)
    }
    
    
    main_list = c(os, ostime, groupvar)
    T1 = hz_cons_measure(main_list, new_var_list, covariate, data_red)
    T2 = hz_dis_measure(main_list, scheme, new_var_list, covariate, data_red)
    T3 = lik_diff_measure(main_list, scheme, covariate, data_red)
    T4 = explain_var_measure(main_list[1:2], scheme, new_var_list, new_var_list2, data_red)
    T5 = balance_measure(main_list[1], scheme, data_red)
    T6 = overall_rank(T1, T2, T3, T4, T5, weight)
    
    mylist <- c(T0, T1, T2, T3, T4, T5, T6)
    return(mylist)
}
