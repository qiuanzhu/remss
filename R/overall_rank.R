#' Overall Ranking.
#' 
#' \code{overall_rank} returns overall ranking for the grouping scheme.  
#' 
#' @param HCM HCM is the hazard consistency measurement results.
#' @param HDM HDM is the hazard discrimination measurement results.
#' @param LDM LDM is the likelihood difference measurement results.
#' @param EVM EVM is the explained variance measurement results.
#' @param BM  BM is the balance measurement results.
#' @param weight weight vector of five measurements. 
#' @return Overall score and overall ranking.
#' @references Xu, W., et al. 'Refining evaluation methodology on TNM stage system: assessment on HPV-related oropharyngeal cancer.' Austin Biom Biostat 2 (2015): 1014.
#' @import  stats
overall_rank <- function(HCM, HDM, LDM, EVM, BM, weight) {
    
    st_sc_HCM = HCM$`Hazard Consistency Measurement`[, c(1, 3)]
    st_sc_HDM = HDM$`Hazard Discrimination Measurement`[, c(1, 3)]
    st_sc_LDM = LDM$`Likelihood Difference Measurement`[, c(1, 3)]
    st_sc_EVM = EVM$`Explained Variance Measurement`[, c(1, 3)]
    st_sc_BM = BM$`Balance Measurement`[, c(1, 3)]
    
    st_sc_HCM = st_sc_HCM[order(st_sc_HCM$Scheme), ]
    st_sc_HDM = st_sc_HDM[order(st_sc_HDM$Scheme), ]
    st_sc_LDM = st_sc_LDM[order(st_sc_LDM$Scheme), ]
    st_sc_EVM = st_sc_EVM[order(st_sc_EVM$Scheme), ]
    st_sc_BM = st_sc_BM[order(st_sc_BM$Scheme), ]
    
    Scheme = st_sc_HCM$Scheme
    Overall_score = st_sc_HCM[2] * weight[1] + st_sc_HDM[2] * weight[2] + st_sc_LDM[2] * weight[3] + st_sc_EVM[2] * weight[4] + 
        st_sc_BM[2] * weight[5]
    
    table <- data.frame(Scheme, Overall_score)[order(Overall_score), ]
    Rank <- c(1:1:length(Scheme))
    colnames(table)[2] <- "Overall Score"
    table <- cbind(table, Rank)
    mylist <- list(`Overall Rank` = table)
    return(mylist)
}
