#' Balance measurement.
#' 
#' \code{balance_measure} returns balance measurement for the grouping scheme.  
#' 
#' @param OS_ind OS_ind is the survival indicator variable.
#' @param stage_list stage_list original of each scheme.
#' @param data Data set. 
#' @return Ranking of balance measurement and its standardized score.
#' @references Xu, W., et al. 'Refining evaluation methodology on TNM stage system: assessment on HPV-related oropharyngeal cancer.' Austin Biom Biostat 2 (2015): 1014.
#' @import  survival
#' @import  stats
#' @export
balance_measure <- function(OS_ind, stage_list, data) {
    
    n <- length(stage_list)
    Scheme <- stage_list
    Score <- rep(0, n)
    Std_Score <- rep(0, n)
    Rank <- c(1:1:n)
    col_OS <- which(colnames(data) == OS_ind)
    # total_case <- sum(data[, col_OS])
    total_case <- length(data[, col_OS])
    
    for (i in 1:n) {
        col_num <- which(colnames(data) == stage_list[i])
        p <- nlevels(data[, col_num])
        num_case <- rep(0, p)
        expect_n <- total_case/p
        for (j in 1:p) {
            # num_case[j] = sum(data[which(data[, col_num] == levels(data[, col_num])[j]), col_OS])
            num_case[j] = length(data[which(data[, col_num] == levels(data[, col_num])[j]), col_OS])
        }
        Score[i] <- mean(abs(num_case - expect_n)/expect_n)
    }
    
    table <- data.frame(Scheme, Score)[order(Score), ]
    diff <- table[n, 2] - table[1, 2]
    for (j in 1:n) {
        Std_Score[j] <- (table[j, 2] - table[1, 2])/diff
    }
    table <- cbind(table, Std_Score, Rank)
    colnames(table) <- c("Scheme", "Score", "Standardized Score", "Rank")
    mylist <- list(`Balance Measurement` = table)
    return(mylist)
    
}


