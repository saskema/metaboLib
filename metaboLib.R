##
## metabolib.R
##
## This library is a reference of functions for untargeted metabolomics.
## It contains the following functions:
##
## - normalise (normalisation of feature matrix by different approaches)
## - batchCorrection (perform batch correction using QC samples)
## - volcano (significance of features based on foldchange
##            and Welch two-sample t-test)
## - evalANOVA (significance of features based von one-way-ANOVA)
## - evalKruskalWallis (significance of features based von Kruskal-Wallis-test)
## - cv (calculate coefficient of variation)
## - iqr (calculate inter quantile range)
## - norm.sum (normalise to sum)
## - norm.mean (normalise to mean)
## - norm.median (normalise to median)
## - norm.iqr (normalise to interquantile range)
## - norm.tic (normalise to total ion count)
## - norm.auto (normalise to zero mean and unit variance)
## - norm.paretro (normalise to zero mean and standard variation)
## - norm.range (normalise to distribution range)
## - ppm (calculate parts per million deviation of two values)
## - indentify.compound (identify a compound based on retention time and exact mass)
##
## The following packages are necessary for this library: ggplot2, gridExtra, ggrepel, MSnbase


##
## Load necessary libraries
##


library("ggplot2")
library("gridExtra")
library("ggrepel")
library("MSnbase")


##
## Functions for untargeted metabolomics
##


## Normalisation of features by different measures
normalise  <- function(matrix,
                       margin,
                       method = c("sum", "mean", "median", "iqr", "tic",
                                  "auto", "paretro", "range", "is"),
                       ref = NA,
                       plot = NA,
                       ...) {


    ## Initial checks
    if(method %in% c("tic", "is")) {

        if(is.na(ref[1])) {

            switch(method,
                   "tic" = stop("Raw files of the analyses are necessary for TIC normalisation!"),
                   "is" = stop("Please specify a row or column number that represents the internal standard!")
                   )
            
        }

    }

    ## Normalise data by method chosen
    matrix.normalised <- switch(method,
                                "sum" = apply(matrix, margin, norm.sum),
                                "mean" = apply(matrix, margin, norm.mean),
                                "median" = apply(matrix, margin, norm.median),
                                "iqr" = apply(matrix, margin, norm.iqr),
                                "tic" = mapply(norm.tic, matrix, ref),
                                "auto" = apply(matrix, margin, norm.auto),
                                "paretro" = apply(matrix, margin, norm.paretro),
                                "range" = apply(matrix, margin, norm.range),
                                "is" = mapply(norm.is, matrix, rep(ref, length(matrix[1,])))
                                )

    ## Reconstruct column names after mapply
    if(method %in% c("tic","is")) {

        rownames(matrix.normalised) <- rownames(matrix)

    }

    ## Plot native and corrected data points
    if (!is.na(plot)) {

        original.long <- as.data.frame(matrix)[1:plot,]
        original.long$name <- rownames(original.long)
        original.long <- reshape(original.long, direction = "long",
                                   varying = names(original.long)[names(original.long) != "name"],
                                   v.names = "area", idvar = "name", timevar = "analysis",
                                   times = names(original.long)[names(original.long) != "name"])

        normalised.long <- as.data.frame(matrix.normalised)[1:plot,]
        normalised.long$name <- rownames(normalised.long)
        normalised.long <- reshape(normalised.long, direction = "long",
                                   varying = names(normalised.long)[names(normalised.long) != "name"],
                                   v.names = "area", idvar = "name", timevar = "analysis",
                                   times = names(normalised.long)[names(normalised.long) != "name"])


        plot.headline <- switch(method,
                              "sum" = "Normalisation by sum of areas",
                              "mean" = "Mean-centering",
                              "median" = "Normalization by median of areas",
                              "iqr" = "Normalization by interquantile range of areas",
                              "tic" = "Normalization by total ion count of analysis",
                              "auto" = "Normalization by zero mean and unit variance ",
                              "paretro" = "Paretro scaling ",
                              "range" = "Range SCaling"
                              )

        
        p1 <- ggplot(data = original.long, aes(x = name, y = area)) +
            geom_boxplot(outlier.colour = "red") +
            ggtitle("Original Data") +
            xlab("Feature") +
            ylab("Area") +
            coord_flip() +
            theme_classic() +
            theme(axis.text.y = element_blank(),
                  axis.ticks = element_blank()
                  )

        p2 <- ggplot(data = normalised.long, aes(x = name, y = area)) +
            geom_boxplot(outlier.colour = "red") +
            ggtitle("Normalised Data") +
            xlab("Feature") +
            ylab("Area") +
            coord_flip() +
            theme_classic() +
            theme(axis.text.y = element_blank(),
                  axis.ticks = element_blank())

        grid.arrange(arrangeGrob(p1, p2, ncol = 2),
                     top = plot.headline)

    }

    return(matrix.normalised)
    
}

## Definition for batch correction
batchCorrection <- function(matrix, # data matrix to perform correction on
                            order, # order of the analysis
                            day, # day of the analysis to allocate subsets of the matrix
                            class, # specifies classes to allocate QC samples
                            sub = NA, # specifies the subset that needs to be corrected
                            poly = 3, # specifies the polynomial of the linear model
                            plot = NA) # specifiesl-pac
    the number of features that will be plottet
{ 

    ## Filter sub, if specified
    if(!is.na(sub)[1]) {

        if(length(sub[sub == TRUE]) == 0) {

            stop("Subset is Empty!")

        }

        
        matrix.uncorrected <- matrix[sub == FALSE,]
        matrix.to.correct <- matrix[sub == TRUE,]

    }

    ## Filter QC measurements
    matrix.pre.corr <- as.data.frame(cbind(order, t(matrix.to.correct)))
    matrix.qc <- matrix.pre.corr[class == "QC",]
    matrix.qc.median <- apply(matrix.qc[,-1], 2, FUN = median)

    ## Calculate difference of each data point to the median of the measurements
    matrix.qc.diff <- as.data.frame(matrix(ncol = length(matrix.qc) - 1,
                                           nrow = length(matrix.qc[,1])))
    for (i in 1:length(matrix.qc.diff)) {

        for (j in 1:length(matrix.qc.diff[,1])) {

            matrix.qc.diff[j,i] <-  matrix.qc.median[i] - (matrix.qc[,-1])[j,i]

        }

    }

    ## Calculate the correction factor for the deviation
    matrix.qc.factor <- as.data.frame(matrix(ncol = length(matrix.qc) - 1,
                                             nrow = length(matrix.qc[,1])))
    for (i in 1:length(matrix.qc.factor)) {

        for (j in 1:length(matrix.qc.factor[,1])) {

            matrix.qc.factor[j,i] <-  matrix.qc.diff[j,i] / matrix.qc.median[i]

        }

    }
    matrix.qc.factor <- cbind(matrix.qc$order, matrix.qc.factor)
    rownames(matrix.qc.factor) <- rownames(matrix.qc)
    colnames(matrix.qc.factor) <- colnames(matrix.qc)

    ## Determine correction factor of remaining data points by linear regression
    matrix.factor <- as.data.frame(matrix(ncol = length(matrix.pre.corr),
                                          nrow = length(matrix.pre.corr[,1])))
    rownames(matrix.factor) <- rownames(matrix.pre.corr)
    colnames(matrix.factor) <- colnames(matrix.pre.corr)
    matrix.factor$order <- order
    for (i in 1:(length(matrix.qc.factor) - 1)) {

        lm <- lm(matrix.qc.factor[,i+1] ~ poly(order, poly), data = matrix.qc.factor)
        matrix.factor[,i + 1] <- predict(lm, newdata = data.frame(order = order))

    }

    ## Correct all data points by their correction factor
    matrix.corrected <- as.data.frame(matrix(ncol = length(matrix.pre.corr),
                                             nrow = length(matrix.pre.corr[,1])))
    rownames(matrix.corrected) <- rownames(matrix.pre.corr)
    colnames(matrix.corrected) <- colnames(matrix.pre.corr)
    matrix.corrected$order <- order
    for (i in 1:(length(matrix.corrected) - 1)) {
        for (j in 1:length(matrix.corrected[,1])) {

            matrix.corrected[j, i + 1] <-  matrix.pre.corr[j, i + 1] +
                matrix.pre.corr[j, i + 1] *
                matrix.factor[j, i + 1]

        }

    }

    ## Plot native and corrected data points
    if (!is.na(plot)) {

        if(plot > length(sub[sub == TRUE])){

            plot  <- length(sub[sub == TRUE])
            
        }

        original.long <- matrix.pre.corr[class == "QC", 2:(1+plot)]
        original.long$Analysis <- rownames(original.long)
        original.long <- reshape(original.long, direction = "long",
                                 varying = names(original.long)[names(original.long) != "Analysis"],
                                 v.names = "area", idvar = "Analysis", timevar = "Feature",
                                 times = names(original.long)[names(original.long) != "Analysis"])

        corrected.long <- matrix.corrected[class == "QC", 2:(1+plot)]
        corrected.long$Analysis <- rownames(corrected.long)
        corrected.long <- reshape(corrected.long, direction = "long",
                                   varying = names(corrected.long)[names(corrected.long) != "Analysis"],
                                   v.names = "area", idvar = "Analysis", timevar = "Feature",
                                   times = names(corrected.long)[names(corrected.long) != "Analysis"])

        p1 <- ggplot(data = original.long, aes(x = Feature, y = area)) +
            geom_boxplot(outlier.colour = "red") +
            ggtitle("Original Data") +
            xlab("Feature") +
            ylab("Area") +
            coord_flip() +
            theme_classic() +
            theme(axis.text.y = element_blank(),
                  axis.ticks = element_blank()
                  )
        p2 <- ggplot(data = corrected.long, aes(x = Feature, y = area)) +
            geom_boxplot(outlier.colour = "red") +
            ggtitle("Batch-Corrected Data") +
            xlab("Feature") +
            ylab("Area") +
            coord_flip() +
            theme_classic() +
            theme(axis.text.y = element_blank(),
                  axis.ticks = element_blank())

        grid.arrange(arrangeGrob(p1, p2, ncol = 2),
                     top = paste0("Batch correction using poly = ", poly))

    }

    ## Return the complete matrix
    return(rbind(matrix.uncorrected, t(matrix.corrected[,-1])))

}

## Definition for univariate statistics
evalVolcano <- function(features,
                        names,
                        classes,
                        class.1,
                        class.2,
                        p.value,
                        foldchange,
                        plot = TRUE,
                        labels = FALSE) {

    features.class.1 <- as.matrix(features[, classes == class.1])
    features.class.2 <- as.matrix(features[, classes == class.2])
    
    summary <- as.data.frame(matrix(data = NA,
                                    nrow = length(features[, 1]),
                                    ncol = 6,
                                    dimnames = list(NULL, c("feature",
                                                            "foldchange",
                                                            "tstat",
                                                            "pvalue",
                                                            "Significant",
                                                            "label"))))
    summary$feature <- names
    
    for (i in 1:length(features[, 1])) {
        
        summary$foldchange[i] <- mean(features.class.1[i, ])/mean(features.class.2[i,])
        
    }
    for (i in 1:length(features[, 1])) {
        
        t <- t.test(features.class.1[i, ], features.class.2[i, ])
        summary$tstat[i] <- t$statistic
        summary$pvalue[i] <- t$p.value
        
    }

    summary$Significant <- ifelse(summary$pvalue < p.value &
                                       (summary$foldchange < -(foldchange) |
                                        summary$foldchange > foldchange), "TRUE", "FALSE")
    summary$label <- names
    summary$label[summary$Significant == "FALSE"] <- ""

    if(plot == TRUE) {

        coloursBoolean <- c("TRUE" = "red", "FALSE" = "grey")

        print(
            ggplot(summary, aes(x = log2(foldchange), y = -log10(pvalue))) +
            geom_point(size = 3, aes(color = Significant)) +
            scale_color_manual(values = coloursBoolean) +
            geom_text_repel(point.padding = 0.2,
                            data = subset(summary, summary$Significant == TRUE),
                            aes(label = summary$feature[summary$Significant == TRUE])) +
            ggtitle(label = paste0("Volcano Plot of ", class.1, " and ", class.2)) +
            xlab("log2(Fold Change)") +
            ylab("-log10(p-Value)") +
            theme_bw() +
            theme(legend.position = "bottom")
        )
        
    }
    
    return(summary)

}

evalANOVA <- function (features,
                       names,
                       classes,
                       p.value = 0.05,
                       plot = TRUE,
                       labels = FALSE) {

    summary <- as.data.frame(matrix(ncol = 4,
                                    nrow = length(names),
                                    dimnames = list(NULL, c("feature",
                                                            "pvalue",
                                                            "Significant",
                                                            "label"))
                                    )
                             )
    
    summary$feature <- names

    for(i in 1:length(features[1,])) {

        sum <- unlist(summary(aov(features[,i]~classes)))
        summary$pvalue[i] <- sum["Pr(>F)1"]
    }

    summary$Significant <- ifelse(summary$pvalue < p.value, "TRUE", "FALSE")
    summary$label <- names
    summary$label[summary$Significant == "FALSE"] <- ""

    if(plot == TRUE) {

        coloursBoolean <- c("TRUE" = "red", "FALSE" = "grey")
        
        if(labels == TRUE) {

            print(
                ggplot(data = summary, aes(x = feature, y = -log10(pvalue),
                                           label = label, colour = Significant)) +
                geom_point(size = 2) +
                scale_color_manual(values = coloursBoolean) +
                geom_hline(yintercept = -log10(p.value), colour = "red") +
                geom_text_repel(point.padding = 0.2, colour = "black") +
                ggtitle("ANOVA Results") +
                xlab("Features") +
                ylab("-log10(p-Values)") +
                theme_classic() +
                theme(axis.text.x = element_blank(),
                      axis.ticks = element_blank()
                      )
            )

        }else {

            print(
                ggplot(data = summary, aes(x = feature, y = -log10(pvalue),
                                           label = label, colour = Significant)) +
                geom_point(size = 2) +
                scale_color_manual(values = coloursBoolean) +
                geom_hline(yintercept = -log10(p.value), colour = "red") +
                ggtitle("ANOVA Results") +
                xlab("Features") +
                ylab("-log10(p-Values)") +
                theme_classic() +
                theme(axis.text.x = element_blank(),
                      axis.ticks = element_blank()
                      )
            )

        }


        }

    return(summary)

}

evalKruskalWallis <- function (features,
                               names,
                               classes,
                               p.value = 0.05,
                               plot = TRUE,
                               labels = FALSE) {

    summary <- as.data.frame(matrix(ncol = 4,
                                    nrow = length(names),
                                    dimnames = list(NULL, c("feature",
                                                            "pvalue",
                                                            "Significant",
                                                            "label"))
                                    )
                             )
    
    summary$feature <- names

    for(i in 1:length(features[1,])) {

        summary$pvalue[i] <- kruskal.test(features[,i], classes)$p.value

    }

    summary$Significant <- ifelse(summary$pvalue < p.value, "TRUE", "FALSE")
    summary$label <- names
    summary$label[summary$Significant == "FALSE"] <- ""

    if(plot == TRUE) {

        coloursBoolean <- c("TRUE" = "red", "FALSE" = "grey")

        if(labels == TRUE) {

            print(
                ggplot(data = summary, aes(x = feature, y = -log10(pvalue),
                                           label = label, colour = Significant)) +
                geom_point(size = 2) +
                scale_color_manual(values = coloursBoolean) +
                geom_hline(yintercept = -log10(p.value), colour = "red") +
                geom_text_repel(point.padding = 0.2, colour = "black") +
                ggtitle("Kruskal-Wallis-Test Results") +
                xlab("Features") +
                ylab("-log10(p-Values)") +
                theme_classic() +
                theme(axis.text.x = element_blank(),
                      axis.ticks = element_blank()
                      )
            )

        } else {

            print(
                ggplot(data = summary, aes(x = feature, y = -log10(pvalue),
                                           label = label, colour = Significant)) +
                geom_point(size = 2) +
                scale_color_manual(values = coloursBoolean) +
                geom_hline(yintercept = -log10(p.value), colour = "red") +
                ggtitle("Kruskal-Wallis-Test Results") +
                xlab("Features") +
                ylab("-log10(p-Values)") +
                theme_classic() +
                theme(axis.text.x = element_blank(),
                      axis.ticks = element_blank()
                      )
            )

        }

        return(summary)
        
    }
    
}

## Perform PCA and plot its characteristics
evalPCA <- function(matrix,
                    plot = TRUE,
                    annotate.scores = TRUE,
                    annotate.loadings = 10,
                    classes,
                    scale = FALSE,
                    center = TRUE) {

    pca <- prcomp(matrix, scale = scale, center = center)
    pca.summary <- summary(pca)
    var <- var(matrix)
    eig <- pca$sdev^2
    kaiser <- length(which(eig > mean(eig)))
    jolliffe <- length(which(eig > (0.7*mean(eig))))
    scree <- data.frame(PC = seq(from = 1, to = length(eig), by = 1),
                        Eigenvalue = eig/sum(eig))
    scree$Importance <- pca.summary$importance[3,]
    scree <- reshape(data = scree, direction = "long",
                     varying = c("Eigenvalue", "Importance"),
                     v.names = "Value", idvar = "PC",
                     timevar = "Characteristic",
                     times = c("Eigenvalue", "Importance"))

    ## Extract loadings to be highlighted in the plot
    loadings.pca <- as.data.frame(pca$rotation[,1:2])
    loadings.pca$abs <- sqrt((loadings.pca$PC1)^2 + (loadings.pca$PC2)^2)
    loadings.pca$label <- rownames(loadings.pca)
    loadings.pca <- arrange(loadings.pca, desc(abs))

    if(annotate.loadings < length(loadings.pca$label)) {

        loadings.pca$label[(annotate.loadings+1):length(loadings.pca$label)] <- "" 

    }

    if(plot == TRUE) {        
        colours_classes <- c("blue", "red", "green", "black")
        names(colours_classes) <- levels(classes)

        if(annotate.scores == TRUE) {

            ## Plot PCA scores with labels
            print(
                ggplot(data = as.data.frame(pca$x[,1:2]), aes(x = PC1, y = PC2)) +
                ##stat_ellipse(type = "norm") +
                geom_point(size = 3, aes(col = Classes)) +
                geom_text_repel(point.padding = 0.2, data = as.data.frame(pca$x[,1:2]),
                                aes(label = rownames(matrix))) +
                scale_color_manual(values = colours_classes) +
                geom_hline(yintercept = 0) +
                geom_vline(xintercept = 0) +
                ggtitle("PCA Scores") +
                xlab(paste0("Principal Component 1", " ", "(", round(pca.summary$importance[2,1]*100, 0), "%)")) +
                ylab(paste0("Principal Component 2", " ", "(", round(pca.summary$importance[2,2]*100, 0), "%)")) +
                theme_bw()
            )

        } else{

            ## Plot PCA scores without labels
            print(
                ggplot(data = as.data.frame(pca$x[,1:2]), aes(x = PC1, y = PC2)) +
                ##stat_ellipse(type = "norm") +
                geom_point(size = 3, aes(col = Classes)) +
                scale_color_manual(values = colours_classes) +
                geom_hline(yintercept = 0) +
                geom_vline(xintercept = 0) +
                ggtitle("PCA Scores") +
                xlab(paste0("Principal Component 1", " ", "(", round(pca.summary$importance[2,1]*100, 0), "%)")) +
                ylab(paste0("Principal Component 2", " ", "(", round(pca.summary$importance[2,2]*100, 0), "%)")) +
                theme_bw()
            )

        }

        ## Plot PCA Loadings
        print(
            ggplot(data = loadings.pca, aes(x = PC1, y = PC2)) +
            geom_point(size = 3, show.legend = FALSE) +
            geom_text_repel(point.padding = 0.2,
                            data = loadings.pca,
                            aes(label = loadings.pca$label)) +
            geom_hline(yintercept = 0) +
            geom_vline(xintercept = 0) +
            ggtitle("PCA Loadings") +
            xlab(paste0("Principal Component 1", " ", "(", round(pca.summary$importance[2,1]*100, 0), "%)")) +
            ylab(paste0("Principal Component 2", " ", "(", round(pca.summary$importance[2,2]*100, 0), "%)")) +
            theme_bw()
        )

        ## Screeplot of eigenvectors per component
        print(
            ggplot(data = scree, aes(x = PC, y = Value, group = Characteristic,
                                     colour = Characteristic)) +
            geom_line(size = 1.5) +
            geom_point(size = 3) +
            scale_x_continuous(breaks = unique(scree$PC), labels = colnames(pca$rotation)) +
            ggtitle(label = "Screeplot",
                    subtitle = paste0("Kaiser = ", kaiser, ", Jolliffe = ", jolliffe)) +
            xlab("Principal Components") +
            ylab("Eigenvalues") +
            theme_bw()
        )

    }

    return(pca)
    
}


##
## Univariate statistics functions
##


## Calculate coefficient of variation
cv <- function(data) {

    sd(data) * 100 / mean(data)
    
}

## Calculate interquantile range
iqr <- function(x) {

    return(quantile(x, probs = 0.75) - quantile(x, 0.25))

}


##
## Normalisation functions
##


## Normalise vector to its sum, assuming constant sum of 1000
## Adopted from https://github.com/xia-lab/MetaboAnalystR/blob/master/R/general_norm_utils.R
norm.sum <- function(x) {

    return(1000*x / sum(x))

}

## Normalise vector to its mean
norm.mean <- function(x) {

    return(x - mean(x))
    
}

## Normalise vector to its median
norm.median <- function(x) {

    return(x - median(x))

}


## Normalise vector to its interquantile range
norm.iqr <- function(x){

    return(x- iqr(x))

}

## Normalise vector to total ion count of the analysis
norm.tic <- function(x, raw.file) {

    raw <- readMSData(raw.file, msLevel = 1, mode  = "onDisk")

        return(x / sum(ionCount(raw)))

}

## Normalise vector to mean of zero and unit variance
norm.auto <- function(x) {

    return((x - mean(x)) / sd(x))

}

## Normalise vector to mean of zero and standard deviation
norm.paretro <- function(x) {

    return((x - mean(x) / sqrt(sd(x))))

}

## Normalise vector to zero mean and range
## Adopted from https://github.com/xia-lab/MetaboAnalystR/blob/master/R/general_norm_utils.R
norm.range <- function(x) {

    if(max(x) == min(x)){
        
        return(x)
        
    } else {
        
        return((x - mean(x))/(max(x)-min(x)))
        
    }
    
}

norm.is <- function(x, is) {

    return(x / x[is])

}

##
## Matrix-specific functions
##

## Perform log10 transformations after replacing zero values by surrogate LOD
log10.matrix <- function(matrix) {

    lod <- min(matrix[matrix != 0])
    matrix[matrix == 0] <- lod
    matrix <- log10(matrix)

}

##
## Identification of compounds in XCMS results
##

## Calculate parts per million
ppm <- function(x, # Comparable quantity
                y) # Measure
{

    ppm <- 10^6 * (x - y) / y

    return(ppm)

}

## Identify compounds in XCMS peaklists (e.g. internal standards)
identify.compound  <-  function(data, # XCMS peaklist
                                compound.name, # name(s) of compound(s)
                                compound.mz, # exact mass(es) of compound(s)
                                compound.rt, # retention time(s) of compound(s)
                                ppmlim = 5, # tolerated mass deviation in parts per million
                                rtlim = 10) # tolerated retention time in seconds
{

    data$Compound <- NA
    
    for(i in 1:length(data[,1])) {

        delta.ppm <- NULL
        delta.rt <- NULL

        for(j in 1:length(compound.name)){

            delta.ppm[j] <- abs(ppm(data$mz[i], compound.mz[j]))
            delta.rt[j] <- abs(data$rt[i] - compound.rt[j])
        }

        index <- which(delta.ppm < ppmlim & delta.rt < rtlim)


        if(length(index) != 0) {

            data$Compound[i] <- as.character(compound.name[index])

        }

    }

    return(data)
    
}
