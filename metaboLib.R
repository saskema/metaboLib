## This program is free software: you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation, either version 3 of the License, or (at your option) any later
## version.

## This program is distributed in the hope that it will be useful, but WITHOUT ANY
## WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
## PARTICULAR PURPOSE.  See the GNU General Public License for more details.

## You should have received a copy of the GNU General Public License along with
## this program.  If not, see https://www.gnu.org/licenses/.

##
##
## metabolib.R
##
##
## This library is a reference of functions for untargeted metabolomics.
## It contains the following functions:
##
## - normalise (normalisation of feature matrix by different approaches)
## - evalVolcano (significance of features based on foldchange and t-tests)
## - evalANOVA (significance of features based von one-way-ANOVA)
## - evalKruskalWallis (significance of features based von Kruskal-Wallis-test)
## - cv (calculate coefficient of variation)
## - iqr (calculate inter quantile range)
## - norm.single.batch (normalize batch drifts)
## - norm.sum (normalise to sum)
## - norm.mean (normalise to mean)
## - norm.median (normalise to median)
## - norm.iqr (normalise to interquantile range)
## - norm.tic (normalise to total ion count)
## - norm.auto (normalise to zero mean and unit variance)
## - norm.paretro (normalise to zero mean and standard variation)
## - norm.range (normalise to distribution range)
## - log10.matrix (log10 transformation after replacing 0s by surrogate LOD)
## - ppm (calculate parts per million deviation of two values)
## - annotate.compound (annotate a compound based on retention time and exact mass)
##
## This library depends on the following packages: ggplot2, gridExtra, ggrepel, MSnbase


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
normalise <- function(matrix, margin, method = c("sum", "mean",
                       "median", "iqr", "tic", "auto", "paretro",
                       "range", "is"), ref = NA, plot = NA, ...) {


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
                                "tic" = norm.tic(x = matrix, margin = margin, raw.files = ref),
                                "auto" = apply(matrix, margin, norm.auto),
                                "paretro" = apply(matrix, margin, norm.paretro),
                                "range" = apply(matrix, margin, norm.range),
                                "is" = norm.is(x = matrix, margin = margin, is = ref)
    )

    ## Reconstruct column names after mapply
    if(method %in%
       c("tic","is")) {

        rownames(matrix.normalised) <- rownames(matrix)

    }

    return(matrix.normalised)

}

## Definition for univariate statistics
evalVolcano <- function(features,
                        names,
                        classes,
                        class.1,
                        class.2,
                        param = TRUE,
                        p.value,
                        fc.thresh,
                        plot = TRUE,
                        labels = FALSE) {

    features.class.1 <- as.matrix(features[, classes == class.1])
    features.class.2 <- as.matrix(features[, classes == class.2])
    
    summary <- as.data.frame(matrix(data = NA,
                                    nrow = length(features[, 1]),
                                    ncol = 7,
                                    dimnames = list(NULL, c("feature",
                                                            "foldchange",
                                                            "log2.foldchange",
                                                            "tstat",
                                                            "pvalue",
                                                            "Significant",
                                                            "label"))))
    summary$feature <- names

    for (i in 1:length(features[, 1])) {
        
        summary$foldchange[i] <- mean(features.class.2[i, ])/mean(features.class.1[i,])
        
    }

    if(param == TRUE) {

        for (i in 1:length(features[, 1])) {


            t <- t.test(features.class.1[i, ], features.class.2[i, ])
            summary$tstat[i] <- t$statistic
            summary$pvalue[i] <- t$p.value

        }
        
    }
    if(param == FALSE) {

        for (i in 1:length(features[, 1])) {

            t <- t.test(features.class.1[i, ], features.class.2[i, ])
            summary$tstat[i] <- t$statistic
            summary$pvalue[i] <- t$p.value

        }
    }

    ## Adopted from MetaboAnalyst
    ## https://github.com/xia-lab/MetaboAnalystR/blob/e7a981fc1c123009ecaacfbbcec7749f0fb530d7/R/stats_univariates.R
    ##
    ## Check if foldchange threshold is greater 1
    
    fc.thresh = ifelse(fc.thresh > 1, fc.thresh, 1 / fc.thresh)
    max.fc.thresh <- log2(fc.thresh)
    min.fc.thresh <- log2(1/fc.thresh)

    summary$log2.foldchange <- log2(summary$foldchange)
    summary$Significant <- ifelse(summary$pvalue < p.value &
                                       (summary$log2.foldchange < min.fc.thresh |
                                        summary$log2.foldchange > max.fc.thresh), "TRUE", "FALSE")
    summary$label <- names
    summary$label[summary$Significant == "FALSE"] <- ""

    if(plot == TRUE) {

        coloursBoolean <- c("TRUE" = "red", "FALSE" = "grey")

        if(labels == "TRUE") {
            print(

                ggplot(summary, aes(x = log2.foldchange, y = -log10(pvalue))) +
                geom_point(size = 3, aes(color = Significant)) +
                scale_color_manual(values = coloursBoolean) +
                geom_hline(yintercept = -log10(p.value), linetype = "dashed") +
                geom_vline(xintercept = min.fc.thresh, linetype = "dashed") +
                geom_vline(xintercept = max.fc.thresh, linetype = "dashed") +
                geom_text_repel(point.padding = 0.2,
                                data = subset(summary, summary$Significant == TRUE),
                                aes(label = summary$feature[summary$Significant == TRUE])) +
                ggtitle(label = paste0("Volcano Plot of ", class.1, " and ", class.2)) +
                xlab("log2(Fold Change)") +
                ylab("-log10(p-Value)") +
                theme_bw() +
                theme(legend.position = "bottom")
            )
        } else {
            print(

                ggplot(summary, aes(x = log2(foldchange), y = -log10(pvalue))) +
                geom_point(size = 3, aes(color = Significant)) +
                scale_color_manual(values = coloursBoolean) +
                geom_hline(yintercept = -log10(p.value), linetype = "dashed") +
                geom_vline(xintercept = log2(fc.thresh), linetype = "dashed") +
                geom_vline(xintercept = log2(1/fc.thresh), linetype = "dashed") +
                ggtitle(label = paste0("Volcano Plot of ", class.1, " and ", class.2)) +
                xlab("log2(Fold Change)") +
                ylab("-log10(p-Value)") +
                theme_bw() +
                theme(legend.position = "bottom")
            )
            
        }

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

        plot

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
                    colours,
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

## Normalise batch drifts in single batch by qc samples
norm.single.batch <- function(matrix,
                                        # matrix; data matrix to perform correction on
                                        # samples in columns, features in rows
                              class,
                                        # character, vector; specifies classes to allocate
                                        # QC samples
                              order,
                                        # numeric, vector; specifies order to correct by
                              type = c("pooled", "surrogate"),
                                        # character; specifies wether to directly correct
                                        # features by QC drift (pooled) or use
                                        # median drift to correct all (surrogate)
                              standard,
                                        # boolean, vector; specifies which features are
                                        # qualified to act as a standard
                              poly,
                                        # numeric; specifies polynomial of the linear model
                              plot = NA)
                                        # numeric; specifies the number of features to plot

{
    
    ## Filter QC samples
    if(type == "pooled") {
        

        matrix.to.correct <- as.data.frame(t(matrix[standard == "TRUE",]))
        matrix.rest <- as.data.frame(t(matrix[standard == "FALSE",]))
        qc <- matrix.to.correct[class == "QC",]
        qc.median <- apply(qc, 2, FUN = median)

        ## Calculate the correction factor based on median
        qc.factor <- as.data.frame(apply(qc, 2, norm.factor))

        ## Determine correction factor of remaining data points by linear regression
        qc.factor$order <- order[class == "QC"]
        correction.factor <- as.data.frame(matrix(ncol = length(matrix.to.correct[1,]),
                                                  nrow = length(matrix.to.correct[,1])))
        rownames(correction.factor) <- rownames(matrix.to.correct)
        colnames(correction.factor) <- colnames(matrix.to.correct)
        correction.factor$order <- order
    


        for (i in 1:(length(qc.factor) - 1)) {

            lm <- lm(qc.factor[,i] ~ poly(order, poly), data = qc.factor)
            correction.factor[,i] <- predict(lm, newdata = data.frame(order = order))

        }

        matrix.corrected <- matrix.to.correct / (1 + correction.factor[names(correction.factor) != "order"])
        matrix.post.corr <- cbind(matrix.corrected, matrix.rest)

        ## Check number of features to plot
        if(plot > length(standard[standard == TRUE])){

            plot  <- length(standard[standard == TRUE])
            
        }
        
    }

    if(type == "surrogate") {
        

        matrix <- as.data.frame(t(matrix))
        qc <- matrix[class == "QC", standard == "TRUE"]
        qc.median <- apply(qc, 2, FUN = median)

        ## Calculate the correction factor based on median
        qc.factor <- as.data.frame(apply(qc, 2, norm.factor))

        ## Determine correction factor of remaining data points by linear regression
        qc.factor$order <- order[class == "QC"]
        correction.factor <- as.data.frame(matrix(ncol = length(matrix[1,standard == "TRUE"]),
                                                  nrow = length(matrix[,1])))
        rownames(correction.factor) <- rownames(matrix)
        colnames(correction.factor) <- colnames(matrix[, standard == "TRUE"])
        correction.factor$order <- order
    


        for (i in 1:(length(qc.factor) - 1)) {

            lm <- lm(qc.factor[,i] ~ poly(order, poly), data = qc.factor)
            correction.factor[,i] <- predict(lm, newdata = data.frame(order = order))

        }

        correction.factor.median <- apply(correction.factor[names(correction.factor) != "order"], 1,
                                          FUN = median)
        matrix.post.corr <- matrix / correction.factor.median

    }

    ## Plot native and corrected data points
    if (!is.na(plot)) {

        original.long <- as.data.frame(t(matrix)[,1:(1+plot)])
        original.long$Analysis <- rownames(original.long)
        original.long <- reshape(original.long, direction = "long",
                                 varying = names(original.long)[names(original.long) != "Analysis"],
                                 v.names = "area", idvar = "Analysis", timevar = "Feature",
                                 times = names(original.long)[names(original.long) != "Analysis"])

        corrected.long <- matrix.post.corr[,1:(1+plot)]
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


    ## Return corrected matrix
    return(t(matrix.post.corr))

}

## Normalise vector to its sum, assuming constant sum of 1000
## Adopted from https://github.com/xia-lab/MetaboAnalystR/blob/master/R/general_norm_utils.R
norm.sum <- function(x) {

    return(1000*x / sum(x))

}

## Normalise vector to its mean
norm.mean <- function(x) {

    return(abs(x - mean(x)))
    
}

## Normalise vector to its median
norm.median <- function(x) {

    return(abs(x - median(x)))

}


## Normalise vector to its interquantile range
norm.iqr <- function(x){

    return(x- iqr(x))

}

## Normalise vector to total ion count of the analysis
norm.tic <- function(x, margin, raw.files) {

    if(margin == 1){
        for(i in 1:nrow(x)) {

            raw <- readMSData(raw.files[i], msLevel = 1, mode  = "onDisk")
            tic <- sum(ionCount(raw))
            x[i, ] <- x[i, ] / tic

        }
    }

    if(margin == 2) {

        for(i in 1:ncol(x)) {

            raw <- readMSData(raw.files[i], msLevel = 1, mode  = "onDisk")
            tic <- sum(ionCount(raw))
            x[, i] <- x[, i] / tic

        }

    }

        return(x)

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

norm.is <- function(x, margin, is) {

    if(margin == 1){
        for(i in 1:nrow(x)) {

            x[i, ] <- x[i, ] / x[i,is]
            
        }
    }

    if(margin == 2) {

        for(i in 1:ncol(x)) {

            x[, i] <- x[, i] / x[is, i]

        }

    }

    return(x)

}

norm.factor <- function(x) {

    return((median(x) - x ) / median(x))

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
## Annotation of compounds in XCMS results
##

## Calculate parts per million
ppm <- function(x, # Comparable quantity
                y) # Measure
{

    ppm <- 10^6 * (x - y) / y

    return(ppm)

}

## Annotate compounds in XCMS peaklists (e.g. internal standards)
annotate.compound  <-  function(data, # XCMS peaklist
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
