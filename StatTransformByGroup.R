StatTransformByGroup <- function(DATA,groups_corresponding,GroupColorMatrix,replicate_scheme)
#This function finds the median of the specified group in replicate scheme and creates a new data frame with the median values
{
    grouping_schemes <- rownames(groups_corresponding)
    grouping_schemes <- grouping_schemes[grouping_schemes!=replicate_scheme]
    groups_corresponding_of_scheme <- groups_corresponding[replicate_scheme,]
    unique_groups <- unique(groups_corresponding[replicate_scheme,])
    n_unique_groups <- length(unique_groups)
    samples <- colnames(DATA)
    GroupColorMatrix <- t(as.matrix(GroupColorMatrix[!rownames(GroupColorMatrix)==replicate_scheme,]))
    groups_corresponding <- t(as.matrix(groups_corresponding[!rownames(groups_corresponding)==replicate_scheme,]))
    n_grouping_schemes <- nrow(groups_corresponding)

    DATA2 <- data.frame(matrix(ncol=0,nrow=nrow(DATA)))
    rownames(DATA2) <- rownames(DATA)

    Color_DF <- data.frame(matrix(nrow=n_grouping_schemes,ncol=n_unique_groups))
    colnames(Color_DF) <- unique_groups
    rownames(Color_DF) <- grouping_schemes
    Groupings_DF <- data.frame(matrix(nrow=n_grouping_schemes,ncol=n_unique_groups))
    colnames(Groupings_DF) <- unique_groups
    rownames(Color_DF) <- grouping_schemes

    for (i in 1:n_unique_groups)
    {
        group <- unique_groups[i]
        member_indices <- groups_corresponding_of_scheme==group
        members <- samples[member_indices]
        median_vector <- apply(DATA[,members],1,median)
        DATA2[,group] <- median_vector
        group_color_column <- which(member_indices)[1]
        Color_DF[,group] <- GroupColorMatrix[,group_color_column]
        Groupings_DF[,group] <- groups_corresponding[,group_color_column]
    }

    GroupColorMatrix <- as.matrix(Color_DF)
    groups_corresponding <- as.matrix(Groupings_DF)
    DATA <- DATA2

    StatTransformByGroup_return <- list(DATA,groups_corresponding,GroupColorMatrix)
}
