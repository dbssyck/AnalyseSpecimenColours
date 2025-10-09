#' @import colorscience
#' @import MASS
#' @import ape
#' @import ggplot2
#' @import openxlsx
NULL

#' Calculates the CIEDE2000 values between every pair of samples in the dataset.
#'
#' @param colour.data Data frame containing sample information including Accession Number (or any unique identifier for each sample), sample grouping(s), and a column each for R, G and B values
#' @param accession.number Name of column with Accession Number
#' @param R.col Name of column containing R values
#' @param G.col Name of column containing G values
#' @param B.col Name of column containing B values
#' @param as.distance.matrix TRUE/FALSE. Defaults to FALSE. TRUE = outputs pairwise distance as R's distance matrix format. FALSE = outputs pairwise distances as mirrored matrix.
#' @export

pairwise.delta <- function (colour.data, accession.number, R.col, G.col, B.col, as.distance.matrix = FALSE) {

  sample.names <- colour.data [ , accession.number]
  n.samples <- length(sample.names)
  rgb.data <- colour.data [ , c (R.col, G.col, B.col) ]

  ## Create empty matrix
  distance.matrix <- matrix(NA, nrow = n.samples, ncol = n.samples,
                            dimnames = list(sample.names, sample.names))

  for (sample1 in 1:(n.samples-1)) {
    for (sample2 in (sample1+1):n.samples) {
      rgb1 <- as.vector(unlist(rgb.data[sample1, ]))
      rgb2 <- as.vector(unlist(rgb.data[sample2, ]))
      distance.matrix[sample1, sample2] <- deltaE2000(rgb1, rgb2)
      distance.matrix[sample2, sample1] <- deltaE2000(rgb1, rgb2) # fill symmetric element
    }
  }

  if ( as.distance.matrix ) { return(as.dist(distance.matrix)) } else { return(distance.matrix) } # convert to distance matrix format if needed
}



#' Check overlap between pairs based on linear value i.e. 75% subspecies rule by Amadon (1949). Default values of 0.75 or threshold.1 and 0.99 for threshold 2, change as required.
#'
#' @param data.1 Vector values from first dataset
#' @param data.1 Vector values from second dataset
#' @param name.data.1 Name of first dataset
#' @param name.data.2 Name of second dataset
#' @param threshold.1 default value 0.75. Threshold value of the percentage of population that must be separable from threshold.2 of the other population. To test the Patten and Unitt (2002) rule, change this value to 0.95 and keep threshold.2 as 0.99.
#' @param threshold.2 default value 0.99. Threshold value of the percentage of population that must be separable from threshold.1 of the other population.
#' @export

check.linear.overlap.Amadon <- function(data.1, data.2, name.data.1, name.data.2, threshold.1 = 0.75, threshold.2 = 0.99) {

  # Determine which group has the higher mean
  is.data.1.higher <- mean(data.1) > mean(data.2)
  higher.data <- if (is.data.1.higher) { data.1 } else { data.2 }
  lower.data  <- if (is.data.1.higher) { data.2 } else { data.1 }
  name.lower.data  <- if (is.data.1.higher) { name.data.2 } else { name.data.1 }

  threshold.lower.1 <- quantile (lower.data, threshold.1)
  threshold.higher.1 <- quantile (higher.data, 1 - threshold.2)
  lower.from.higher <- ifelse(threshold.lower.1 > threshold.higher.1, "No", "Yes")

  threshold.lower.2 <- quantile (lower.data, threshold.2)
  threshold.higher.2 <- quantile (higher.data, 1 - threshold.1)
  higher.from.lower <- ifelse(threshold.lower.2 > threshold.higher.2, "No", "Yes")

  diagnosable.both.ways <- ifelse( lower.from.higher == "No" | higher.from.lower == "No", "No", "Yes")

  diagnosable <- data.frame ("Group 1" = name.data.1,
                             "Group 2" = name.data.2,
                             "Lower Group" = name.lower.data,
                             "Lower from higher" = lower.from.higher,
                             "Higher from lower" = higher.from.lower,
                             "Diagnosable from each other" = diagnosable.both.ways)

  return(diagnosable)
}


#' Runs check.linear.overlap.Amadon() on all pairwise combinations of groups for one parameter. Default values of 0.75 or threshold.1 and 0.99 for threshold 2, change as required.
#'
#' @param colour.data Data frame containing sample information including Accession Number (or any unique identifier for each sample), sample grouping(s), and a column each for R, G and B values
#' @param accession.number Name of column with Accession Number
#' @param population Name of column with population assignment (e.g. subspecies, general locality etc.)
#' @param parameter.to.check.overlap name of column with values to test subspecies rule by Amadon (1949) on
#' @param threshold.1 default value 0.75. Threshold value of the percentage of population that must be separable from threshold.2 of the other population. To test the Patten and Unitt (2002) rule, change this value to 0.95 and keep threshold.2 as 0.99.
#' @param threshold.2 default value 0.99. Threshold value of the percentage of population that must be separable from threshold.1 of the other population.
#' @export

pairwise.group.Amadon <- function (colour.data, accession.number, population, parameter.to.check.overlap, threshold.1 = 0.75, threshold.2 = 0.99) {

  subspecies.list <- unique (colour.data[, population])
  subspecies.combination <- expand.grid(subspecies.list, subspecies.list) # get all pairs of subspecies
  subspecies.combination <- subspecies.combination [ subspecies.combination[,1] != subspecies.combination[,2], ] # remove repeated pairs

  diagnosable <- data.frame (matrix(data = NA, nrow = nrow(subspecies.combination), ncol = 6))
  colnames(diagnosable) <- c("Group 1", "Group 2", "Lower Group", "Lower from higher", "Higher from lower", "Diagnosable from each other")

  for (pair in 1:nrow(subspecies.combination)) {
    subspecies.name.1 <- as.vector(subspecies.combination [pair, 1])
    subspecies.name.2 <- as.vector(subspecies.combination [pair, 2])
    values.subspecies.1 <- colour.data [ colour.data[, population] %in% subspecies.combination[pair,1], parameter.to.check.overlap ]
    values.subspecies.2 <- colour.data [ colour.data[, population] %in% subspecies.combination[pair,2], parameter.to.check.overlap ]

    diagnosable[pair, ] <- check.linear.overlap.Amadon (values.subspecies.1, values.subspecies.2, subspecies.name.1, subspecies.name.2, threshold.1, threshold.2)
  }

  return(diagnosable)
}


#' Create neighbour-joining tree based on RGB input and produce phylogram
#'
#' @param colour.data Data frame containing sample information including Accession Number (or any unique identifier for each sample), sample grouping(s), and a column each for R, G and B values
#' @param accession.number Name of column with Accession Number
#' @param population Name of column with population assignment (e.g. subspecies, general locality etc.)
#' @param R.col Name of column containing R values
#' @param G.col Name of column containing G values
#' @param B.col Name of column containing B values
#' @param export.figures TRUE/FALSE.  Defaults to FALSE. TRUE exports the phylogram as an image file.
#' @param output.folder Path to folder where image file will be exported into. If left blank, the file will be created in present working directory. If a folder with the name designated here doesn't exist in the present workign directory, it will be created.
#' @param output.file.prefix Specify file name of phylogram image file if desired. Otherwise the file will be nj_phylogram.jpg
#' @export

neighbour.joining <- function(colour.data, accession.number, population, R.col, G.col, B.col,
                              export.figures = FALSE, output.folder = NULL, output.file.prefix = NULL) {

  # Create pairwise distance matrix (needs to be in the distance matrix format and not mirrored by a diagonal 0)
  pairwise.distance.m <- pairwise.delta (colour.data, accession.number, R.col, G.col, B.col, as.distance.matrix = 1)

  # Create neighbour-joining trees
  nj.tree <- njs(pairwise.distance.m) # input here is pairwise distance as a distance matrix

  # Add subspecies (or group identity) to specimen label
  nj.tree$tip.label <- paste (colour.data[[population]],
                              colour.data[[accession.number]],
                              sep = "_") # add subspecies to specimen label

  # Root at longest internal node
  root.longest.internal <- function (tree.to.root) {
    internal_edges <- which(tree.to.root$edge[, 2] > Ntip(tree.to.root))
    longest_internal_edge <- internal_edges[which.max(tree.to.root$edge.length[internal_edges])]
    nodes <- tree.to.root$edge[longest_internal_edge, ]
    tree.to.root <- root(tree.to.root, node = nodes[2], resolve.root = TRUE)

    return(tree.to.root)
  }

  nj.tree.internal.root <- root.longest.internal(nj.tree)

  if (export.figures) {
    # Determine file size to plot
    n.samples <- nrow(colour.data)
    height.per.sample <- 20  # in pixels or units
    min.height <- 800        # minimum height in pixels
    plot.height <- max(min.height, height.per.sample * n.samples) # Calculate height to plot

    # Create output file path
    if (is.null(output.folder)) { output.folder <- getwd() } # if no output.folder designated, create files in working directory
    if (!dir.exists(output.folder)) {dir.create(output.folder, recursive = TRUE)}

    # Create output filename
    if (is.null(output.file.prefix)) { output.file.prefix <- "nj_phylogram" }

    # Plot trees
    jpeg(filename = paste(output.folder,"/",output.file.prefix,".jpg",sep=""), width = 2000, height = plot.height)
    plot(nj.tree.internal.root, type = "phylogram")
    dev.off()
  } else { plot(nj.tree.internal.root, type = "phylogram")}
}


#' Convert to other colour models starting from RGB
#'
#' @param colour.data Data frame containing sample information including Accession Number (or any unique identifier for each sample), sample grouping(s), and a column each for R, G and B values
#' @param R.col Name of column containing R values
#' @param G.col Name of column containing G values
#' @param B.col Name of column containing B values

add.colour.models.from.RGB <- function(colour.data, R.col, G.col, B.col) {
  rgb <- colour.data[,c(R.col, G.col, B.col)]
  colnames(rgb) <- c("R", "G", "B")
  greyscale <- data.frame ( 0.299 * rgb[,R.col] + 0.587 * rgb[,G.col]  + 0.114 * rgb[,B.col])
  colnames(greyscale) <- "greyscale"
  cmy <- data.frame (RGB2CMY(rgb))
  hsl <- data.frame (RGB2HSL(rgb))
  hsv <- data.frame (RGB2HSV(rgb))

  return(list(rgb=rgb,
              cmy=cmy,
              hsl=hsl,
              hsv=hsv,
              greyscale = greyscale))
}


#' Test 75% subspecies rule after running PCA and LDA on colours. Can be done for all colour models if the input for colour.data.from.all.colour.models is the results of add.colour.models.from.RGB()
#'
#' @param colour.data Data frame containing sample information including Accession Number (or any unique identifier for each sample), sample grouping(s), and a column each for R, G and B values
#' @param accession.number Name of column with Accession Number
#' @param population Name of column with population assignment (e.g. subspecies, general locality etc.)
#' @param R.col Name of column containing R values
#' @param G.col Name of column containing G values
#' @param B.col Name of column containing B values
#' @param threshold.1 default value 0.75. Threshold value of the percentage of population that must be separable from threshold.2 of the other population. To test the Patten and Unitt (2002) rule, change this value to 0.95 and keep threshold.2 as 0.99.
#' @param threshold.2 default value 0.99. Threshold value of the percentage of population that must be separable from threshold.1 of the other population.
#' @param colour.data.from.all.colour.models If running the subspecies rule by Amadon (1949) for other colour models, the input for this will be the results from add.colour.models.from.RGB()
#' @param export.figures TRUE/FALSE.  Defaults to FALSE. TRUE exports the PCAs and LDAs as image files.
#' @param output.folder Path to folder where image file will be exported into. If left blank, the file will be created in present working directory. If a folder with the name designated here doesn't exist in the present workign directory, it will be created.
#' @export

pca.lda.overlap <- function(colour.data, accession.number, population, R.col, G.col, B.col,
                            threshold.1 = 0.75, threshold.2 = 0.99,
                            colour.data.from.all.colour.models = NULL,
                            export.figures = FALSE,
                            output.folder = NULL) {

  subspecies.list <- unique (colour.data[, population])
  subspecies.combination <- expand.grid(subspecies.list, subspecies.list) # get all pairs of subspecies (or other defined group)
  subspecies.combination <- subspecies.combination [ subspecies.combination[,1] != subspecies.combination[,2], ] # remove repeated pairs

  n.colour.models <- if (is.null(colour.data.from.all.colour.models)) {1} else {length(colour.data.from.all.colour.models)}

  pca.overlap.check <- vector (length = n.colour.models, mode = "list") # create list to store PCA diagnosability
  names(pca.overlap.check) <-  if (is.null(colour.data.from.all.colour.models)) {"rgb"} else {names(colour.data.from.all.colour.models)}

  lda.overlap.check <- pca.overlap.check # create identical list to store LDA diagnosability

  # Create output file path if needed
  if (export.figures) {
  if (is.null(output.folder)) { output.folder <- getwd() } # if no output.folder designated, create files in working directory
  if (!dir.exists(output.folder)) {dir.create(output.folder, recursive = TRUE)}
  }

  for (colour.model in 1:n.colour.models) {

    # Load colour data
    if (n.colour.models == 1) { colour.model.used <- colour.data[, c(R.col, G.col, B.col)] } else
      { colour.model.used <- colour.data.from.all.colour.models[[colour.model]] }
    name.colour.model.used <- ifelse (n.colour.models == 1, "rgb", names(colour.data.from.all.colour.models)[colour.model])

    if (ncol (colour.model.used) == 1) {break} # ignore greyscale because there's only 1 variable and PCA/LDA cannot be conducted

    # colour.model.diagnose.pca <- vector (length = ncol (colour.model.used), mode = "list") # DELETE IF NOT NEEDED
    # names(colour.model.diagnose.pca) <- colnames(colour.model.used) # DELETE IF NOT NEEDED

    # Run PCA and LDA on entire dataset for visualisation
    pca.colour.model.all.samples <- prcomp(colour.model.used, scale = T) # Run PCA on the colour components

    if (export.figures) {

      pca.data.for.plot <- data.frame(pca.colour.model.all.samples$x[,1:3], population = factor(colour.data[,population]))
      pca.summary <- summary( pca.colour.model.all.samples )

      pca.plot.of.all <- ggplot(pca.data.for.plot) +
        geom_point(aes(PC1, PC2, colour = population), size = 2.5) +
        labs (x = paste ("PC1 (", sprintf ("%.1f", signif (pca.summary$importance[2,1] * 100, 3)), "%)", sep = ""),
              y = paste ("PC2 (", sprintf ("%.1f", signif (pca.summary$importance[2,2] * 100, 3)), "%)", sep = "")) +
        stat_ellipse(data = pca.data.for.plot,
                     aes (x = PC1, y = PC2, fill = population), geom="polygon", level=0.95, alpha=0.2)

      pca.plot.of.all.pc2n3 <- ggplot(pca.data.for.plot) +
        geom_point(aes(PC1, PC2, colour = population), size = 2.5) +
        labs (x = paste ("PC1 (", sprintf ("%.1f", signif (pca.summary$importance[2,2] * 100, 3)), "%)", sep = ""),
              y = paste ("PC2 (", sprintf ("%.1f", signif (pca.summary$importance[2,3] * 100, 3)), "%)", sep = "")) +
        stat_ellipse(data = pca.data.for.plot,
                     aes (x = PC2, y = PC3, fill = population), geom="polygon", level=0.95, alpha=0.2)

      png(filename = paste (output.folder,"/PCA_allsamples_PC1n2_",
                            name.colour.model.used,"_",
                            population,
                            ".png", sep=""), width = 1000, height = 1000)

      print (pca.plot.of.all)

      dev.off()

      png(filename = paste (output.folder,"/PCA_allsamples_PC2n3_",
                            name.colour.model.used,"_",
                            population,
                            ".png", sep=""), width = 1000, height = 1000)

      print (pca.plot.of.all.pc2n3)

      dev.off()
    } # plot PCA for all samples

    colour.model.used.lda <- cbind (colour.data[ ,population], colour.model.used)
    colnames(colour.model.used.lda)[1] <- population
    formula.lda <- as.formula(paste(population, "~ .", sep = ""))
    lda.colour.model.all.samples <- lda(formula.lda, data = colour.model.used.lda)
    ld.prop <- vector (length = length(lda.colour.model.all.samples$svd))
    ld.prop <- 100 * (lda.colour.model.all.samples$svd^2 / sum(lda.colour.model.all.samples$svd^2))

    colour.lda.values <- predict(lda.colour.model.all.samples)

    if (export.figures) {
      lda.for.plot <- data.frame(population = colour.model.used.lda[,population],
                                 lda = colour.lda.values$x)

      lda.plot.of.all <- ggplot(lda.for.plot) +
        geom_point(aes(lda.LD1, lda.LD2, colour = population), size = 2.5) +
        labs (x = paste ("LD1 (", signif (ld.prop[1], 3), "%)", sep = ""),
              y = paste ("LD2 (", signif (ld.prop[2], 3), "%)", sep = "") ) +
        stat_ellipse(data = lda.for.plot,
                     aes (x = lda.LD1, y = lda.LD2, fill = population), geom="polygon", level=0.95, alpha=0.2)


      png(filename = paste (output.folder,"/LDA_allsamples_LD1n2_",
                            name.colour.model.used,"_",
                            population,
                            ".png", sep=""), width = 1000, height = 1000)

      print (lda.plot.of.all)

      dev.off()

      if("lda.LD3" %in% colnames(lda.for.plot)) {
        lda.plot.of.all.ld2n3 <- ggplot(lda.for.plot) +
          geom_point(aes(lda.LD2, lda.LD3, colour = population), size = 2.5) +
          labs (x = paste ("LD2 (", signif (ld.prop[2], 3), "%)", sep = ""),
                y = paste ("LD3 (", signif (ld.prop[3], 3), "%)", sep = "") ) +
          stat_ellipse(data = lda.for.plot,
                       aes (x = lda.LD1, y = lda.LD2, fill = population), geom="polygon", level=0.95, alpha=0.2)

        png(filename = paste (output.folder,"/LDA_allsamples_LD2n3_",
                              name.colour.model.used,"_",
                              population,
                              ".png", sep=""), width = 1000, height = 1000)

        print (lda.plot.of.all.ld2n3)

        dev.off()
      }
    } # plot LDA for all samples

    # Run PCA and LDA on every pair of groups and check if they are diagnosable according to Amadon (1949) criteria
    diagnosable <- data.frame (matrix(data = NA, nrow = nrow(subspecies.combination), ncol = 6))
    colnames(diagnosable) <- c("Group 1", "Group 2", "Lower Group", "Lower from higher", "Higher from lower", "Diagnosable from each other")

    pca.store <- rep (list(diagnosable), ncol (pca.colour.model.all.samples$x))
    names(pca.store) <- paste("PC", seq(1:ncol (pca.colour.model.all.samples$x)), sep = "")

    lda.store <- rep (list(diagnosable), 1)
    names(lda.store) <- "LD1"

    for (pair in 1:nrow(subspecies.combination)) {
      subspecies.name.1 <- as.vector(subspecies.combination [pair, 1])
      subspecies.name.2 <- as.vector(subspecies.combination [pair, 2])
      names.subspecies.for.index <- colour.data [colour.data[, population] %in% c(subspecies.name.1, subspecies.name.2),
                                                 population]
      colour.model.used.subspecies.pair <- colour.model.used [ colour.data[, population] %in% c(subspecies.name.1, subspecies.name.2), ] # get colour values of those relevant to species pair being tested

      pca.colour.model.subspecies.pair <- prcomp(colour.model.used.subspecies.pair, scale = T) # run PCA

      if (export.figures) {
        pca.data.for.plot <- data.frame(pca.colour.model.subspecies.pair$x[,1:3], population = factor(names.subspecies.for.index))
        pca.summary <- summary( pca.colour.model.subspecies.pair )

        pca.plot.of.pair <- ggplot(pca.data.for.plot) +
          geom_point(aes(PC1, PC2, colour = population), size = 2.5) +
          labs (x = paste ("PC1 (", sprintf ("%.1f", signif (pca.summary$importance[2,1] * 100, 3)), "%)", sep = ""),
                y = paste ("PC2 (", sprintf ("%.1f", signif (pca.summary$importance[2,2] * 100, 3)), "%)", sep = "")) +
          stat_ellipse(data = pca.data.for.plot,
                       aes (x = PC1, y = PC2, fill = population), geom="polygon", level=0.95, alpha=0.2)

        pca.plot.of.pair.pc2n3 <- ggplot(pca.data.for.plot) +
          geom_point(aes(PC1, PC2, colour = population), size = 2.5) +
          labs (x = paste ("PC1 (", sprintf ("%.1f", signif (pca.summary$importance[2,2] * 100, 3)), "%)", sep = ""),
                y = paste ("PC2 (", sprintf ("%.1f", signif (pca.summary$importance[2,3] * 100, 3)), "%)", sep = "")) +
          stat_ellipse(data = pca.data.for.plot,
                       aes (x = PC2, y = PC3, fill = population), geom="polygon", level=0.95, alpha=0.2)

        png(filename = paste (output.folder,"/PC1n2_",
                              name.colour.model.used,"_",
                              population,"_",
                              subspecies.name.1, "_", subspecies.name.2,
                              ".png", sep=""), width = 1000, height = 1000)

        print (pca.plot.of.pair)

        dev.off()

        png(filename = paste (output.folder,"/PC2n3_",
                              name.colour.model.used,"_",
                              population,"_",
                              subspecies.name.1, "_", subspecies.name.2,
                              ".png", sep=""), width = 1000, height = 1000)

        print (pca.plot.of.pair.pc2n3)

        dev.off()
      } # plot PCA for subspecies pair

      for (PC in 1:ncol(pca.colour.model.subspecies.pair$x)) {
        values.subspecies.1 <- pca.colour.model.subspecies.pair$x [ which (names.subspecies.for.index %in% subspecies.combination[pair,1]), PC ]
        values.subspecies.2 <- pca.colour.model.subspecies.pair$x [ which (names.subspecies.for.index %in% subspecies.combination[pair,2]), PC ]

        pca.store[[PC]][pair, ] <- check.linear.overlap.Amadon (values.subspecies.1, values.subspecies.2, subspecies.name.1, subspecies.name.2, threshold.1, threshold.2)
      } # Run diagnosability check on PC loadings

      colour.model.used.lda.subspecies.pair <- cbind (names.subspecies.for.index, colour.model.used.subspecies.pair)
      colnames(colour.model.used.lda.subspecies.pair)[1] <- population
      formula.lda <- as.formula(paste(population, "~ .", sep = ""))
      lda.colour.model.subspecies.pair <- lda(formula.lda, data = colour.model.used.lda.subspecies.pair)
      ld.prop.subspecies.pair <- vector (length = length(lda.colour.model.all.samples$svd))
      ld.prop.subspecies.pair <- 100 * (lda.colour.model.all.samples$svd^2 / sum(lda.colour.model.all.samples$svd^2))

      colour.lda.values.subspecies.pair <- predict(lda.colour.model.subspecies.pair)

      if (export.figures) {
        lda.for.plot.subspecies.pair <- data.frame(lda = colour.lda.values.subspecies.pair$x, population = colour.model.used.lda.subspecies.pair[,population])

        lda.plot.of.pair <- ggplot(lda.for.plot.subspecies.pair, aes(x = population, y = LD1, fill = population)) +
          geom_boxplot() +  # Create the box plot
          labs(title = "LDA LD1 Values by Group", x = "Group", y = paste ("LD1 (", signif (ld.prop.subspecies.pair[1], 3), "%)", sep = "")) +
          theme_minimal() +  # Clean minimal theme
          scale_fill_manual(values = c("blue", "red"))

        png(filename = paste (output.folder,"/LD1_",
                              name.colour.model.used,"_",
                              population,"_",
                              subspecies.name.1, "_", subspecies.name.2,
                              ".png", sep=""), width = 1000, height = 1000)

        print (lda.plot.of.pair)

        dev.off()
      } # plot LD1 loadings for subspecies pair

      for (LD in 1:ncol(colour.lda.values.subspecies.pair$x)) {
        values.subspecies.1 <- colour.lda.values.subspecies.pair$x [ which (names.subspecies.for.index %in% subspecies.combination[pair,1]), LD ]
        values.subspecies.2 <- colour.lda.values.subspecies.pair$x [ which (names.subspecies.for.index %in% subspecies.combination[pair,2]), LD ]

        lda.store[[LD]][pair, ] <- check.linear.overlap.Amadon (values.subspecies.1, values.subspecies.2, subspecies.name.1, subspecies.name.2, threshold.1, threshold.2)
      } # Run diagnosability check on LD loading (only one since a pair is being compared)

    }
    pca.overlap.check[[colour.model]] <- pca.store
    lda.overlap.check[[colour.model]] <- lda.store
  }

  pca.overlap.check <- Filter(Negate(is.null), pca.overlap.check)
  lda.overlap.check <- Filter(Negate(is.null), lda.overlap.check)

  return(list (pca.check = pca.overlap.check,
               lda.check = lda.overlap.check))
}

#' Helper function: do not export. Flatten nested lists
#'
#' @param list.to.flatten List with sublists that need to be flattened
#' @param parent.name Name to rename flattened lists to

flatten.list <- function(list.to.flatten, parent.name = "") {
  flat.list <- list()

  for (i in seq_along(list.to.flatten)) {
    name <- names(list.to.flatten)[i]
    if (is.null(name) || name == "") { name <- paste0("Unnamed", i) }   # For missing names
    new.name <- ifelse(parent.name == "", name, paste(parent.name, name, sep = "."))

    item <- list.to.flatten[[i]]

    if (is.data.frame(item)) { flat.list[[new.name]] <- item # Store dataframe with unique name
    } else if (is.list(item)) {
      flat.list <- c(flat.list, flatten.list(item, new.name)) # Continue working into sublists
    }
  }

  return(flat.list)
}


#' Helper function: do not export. Calculate overlap for individual colours based on all colour models
#' @param colour.data Data frame containing sample information including Accession Number (or any unique identifier for each sample), sample grouping(s), and a column each for R, G and B valuesusethis::use_git_config(user.name = "Your Name", user.email = "youremail@example.com")
#' @param colour.data.from.all.colour.models Results from add.colour.models.from.RGB()
#' @param accession.number Name of column with Accession Number
#' @param population Name of column with population assignment (e.g. subspecies, general locality etc.)
#' @param threshold.1 default value 0.75. Threshold value of the percentage of population that must be separable from threshold.2 of the other population. To test the Patten and Unitt (2002) rule, change this value to 0.95 and keep threshold.2 as 0.99.
#' @param threshold.2 default value 0.99. Threshold value of the percentage of population that must be separable from threshold.1 of the other population.

calculate.overlap.for.all.colours <- function(colour.data, colour.data.from.all.colour.models,
                                              accession.number, population, # need colour.data for metadata
                                              threshold.1 = 0.75, threshold.2 = 0.99) {

  colour.model.overlap.check <- vector (length = length (colour.data.from.all.colour.models), mode = "list") # create list to store colour diagnosability
  names(colour.model.overlap.check) <- names(colour.data.from.all.colour.models)

  for (colour.model in 1:length(colour.data.from.all.colour.models)) {

    colour.model.used <- colour.data.from.all.colour.models[[colour.model]]
    colour.model.used.name <- names (colour.data.from.all.colour.models) [colour.model]
    colour.model.diagnose <- vector (length = ncol (colour.model.used), mode = "list")
    names(colour.model.diagnose) <- colnames(colour.model.used)

    for (colour.parameter in 1:ncol(colour.model.used)) {
      column.name <- colnames(colour.model.used)[colour.parameter]
      dataframe.to.run.pairwise <- cbind (colour.data [, c(accession.number, population)], colour.model.used)
      colour.model.diagnose[[colour.parameter]] <- pairwise.group.Amadon (dataframe.to.run.pairwise,
                                                                          accession.number,
                                                                          population,
                                                                          column.name, threshold.1, threshold.2)
    }
    colour.model.overlap.check[[colour.model]] <- colour.model.diagnose
  }

  return(colour.model.overlap.check)
}


#' Helper function: do not export. Plot RGB and CYM colours
#' @param colour.data Data frame containing sample information including Accession Number (or any unique identifier for each sample), sample grouping(s), and a column each for R, G and B values
#' @param colour.data.from.all.colour.models Results from add.colour.models.from.RGB()
#' @param population Name of column with population assignment (e.g. subspecies, general locality etc.)
#' @param output.folder Path to folder where image file will be exported into. If left blank, the file will be created in present working directory. If a folder with the name designated here doesn't exist in the present workign directory, it will be created.

plot.colours <- function(colour.data, colour.data.from.all.colour.models, population, output.folder) {  # need colour.data for metadata also

  RGB.data <- colour.data.from.all.colour.models[["rgb"]]
  CYM.data <- colour.data.from.all.colour.models[["cmy"]]

  RGB.data.long <- data.frame(
    subspecies = rep(colour.data[,population], times = 3),
    value = c(RGB.data$R, RGB.data$G, RGB.data$B),
    colour = rep(c("R", "G", "B"), each = nrow(colour.data)))

  CYM.data.long <- data.frame(
    subspecies = rep(colour.data[,population], times = 3),
    value = c(CYM.data$C, CYM.data$Y, CYM.data$M),
    colour = rep(c("C", "Y", "M"), each = nrow(colour.data)))

  RGB.plot <- ggplot(RGB.data.long, aes(x = factor(colour, levels = c("R", "G", "B")), y = value, fill = colour)) +
    geom_boxplot() +
    facet_wrap(~subspecies) +  # Separate plots for each group
    labs(title = paste("Boxplot of CYM values by ",population, sep = ""), x = "Colour", y = "Value") +
    scale_fill_manual(values = c("R" = "red", "G" = "green", "B" = "blue")) +  # Set colors for R, G, B
    theme_minimal()

  CYM.plot <- ggplot(CYM.data.long, aes(x = factor(colour, levels = c("C", "Y", "M")), y = value, fill = colour)) +
    geom_boxplot() +
    facet_wrap(~subspecies) +  # Separate plots for each group
    labs(title = paste("Boxplot of CYM values by ",population, sep = ""), x = "Colour", y = "Value") +
    scale_fill_manual(values = c("C" = "cyan", "Y" = "yellow", "M" = "magenta")) +  # Set colors for R, G, B
    theme_minimal()

  png(filename = paste (output.folder,"/box_colours_RGB_",
                        population,
                        ".png", sep=""), width = 1000, height = 1000)

  print (RGB.plot)

  dev.off()

  png(filename = paste (output.folder,"/box_colours_CYM_",
                        population,
                        ".png", sep=""), width = 1000, height = 1000)

  print (CYM.plot)

  dev.off()
}



#' Conduct all analyses described in Teo et al. (2025). Choose between key analyses only or all analyses (i.e. including exploratory analyses)
#'
#' @param colour.data Data frame containing sample information including Accession Number (or any unique identifier for each sample), sample grouping(s), and a column each for R, G and B values
#' @param accession.number Name of column with Accession Number
#' @param population Name of column with population assignment (e.g. subspecies, general locality etc.)
#' @param R.col Name of column containing R values
#' @param G.col Name of column containing G values
#' @param B.col Name of column containing B values
#' @param run "key" or "all". "key" runs the analyses found to be important in Teo et al. (2025), "all" runs analyses including exploratory analyses that were conducted
#' @param threshold.1 default value 0.75. Threshold value of the percentage of population that must be separable from threshold.2 of the other population. To test the Patten and Unitt (2002) rule, change this value to 0.95 and keep threshold.2 as 0.99.
#' @param threshold.2 default value 0.99. Threshold value of the percentage of population that must be separable from threshold.1 of the other population.
#' @param export.colours.plot TRUE/FALSE. Defaults to FALSE. TRUE exports boxplots of individual colour values (e.g. R, G and B) for all colour models, for all pairwise group comparisons
#' @param export.figures TRUE/FALSE. Defaults to FALSE. TRUE exports the LDA and PCA figures into an output file.
#' @param export.excel TRUE/FALSE. TRUE exports results of the Amadon (1949) subspecies rule test as an .xlsx file.
#' @param output.folder Path to folder where image file will be exported into. If left blank, the file will be created in present working directory. If a folder with the name designated here doesn't exist in the present workign directory, it will be created.
#' @param excel.file.name Specify file name of excel file (if export.excel = TRUE). Otherwise the file will be named diagnosability.xlsx
#' @export

run.functions <- function(colour.data, accession.number, population, R.col, G.col, B.col,
                          run,
                          threshold.1 = 0.75, threshold.2 = 0.99,
                          export.colours.plot = FALSE, export.figures = FALSE, export.excel = FALSE,
                          output.folder = NULL, excel.file.name = NULL) {

  # Run key functions only or all functions?
  if (!(run %in% c("all","key"))) { stop ( "The run argument (the argument after B.col, which is the 7th argument) needs to be either \"key\" or \"all\" ") }

  # Create output file path if needed
  if (export.figures | export.colours.plot | export.excel ) {
    if (is.null(output.folder)) { output.folder <- getwd() } # if no output.folder designated, create files in working directory
    if (!dir.exists(output.folder)) {dir.create(output.folder, recursive = TRUE)}
  }

  # Generate other colour models if needed
  if (run == "all") {
    colour.data.from.all.colour.models <- add.colour.models.from.RGB(colour.data, R.col, G.col, B.col)
  }

  # Plot colours if needed
  if (run == "all") {
    if (export.colours.plot) { plot.colours(colour.data, colour.data.from.all.colour.models, population, output.folder) }
  }

  # Run individual colours overlap check
  if (run == "all") {
    overlap.colours <- calculate.overlap.for.all.colours (colour.data, colour.data.from.all.colour.models, accession.number, population,
                                                          threshold.1 = threshold.1, threshold.2 = threshold.2)
  }

  # Run PCA and LDA check
  if (run == "all") {
    if (export.figures) {
      overlap.pca.lda <- pca.lda.overlap(colour.data, accession.number, population,
                                         R.col, G.col, B.col,
                                         threshold.1 = threshold.1, threshold.2 = threshold.2,
                                         colour.data.from.all.colour.models = colour.data.from.all.colour.models,
                                         export.figures = export.figures, output.folder = output.folder) }
    else { overlap.pca.lda <- pca.lda.overlap(colour.data, accession.number, population,
                                              R.col, G.col, B.col,
                                              colour.data.from.all.colour.models = colour.data.from.all.colour.models,
                                              threshold.1 = threshold.1, threshold.2 = threshold.2) }
  }

  if (run == "key") {
    if (export.figures) {
      overlap.pca.lda <- pca.lda.overlap(colour.data, accession.number, population,
                                         R.col, G.col, B.col,
                                         threshold.1 = threshold.1, threshold.2 = threshold.2,
                                         export.figures = export.figures, output.folder = output.folder) }
    else { overlap.pca.lda <- pca.lda.overlap(colour.data, accession.number, population,
                                              R.col, G.col, B.col,
                                              threshold.1 = threshold.1, threshold.2 = threshold.2) }
  }

  # Calculate overlap for individual colours
  if (run == "all") {
    individual.colours.overlap <- calculate.overlap.for.all.colours (colour.data, colour.data.from.all.colour.models,
                                                                     accession.number, population,
                                                                     threshold.1 = threshold.1, threshold.2 = threshold.2)
  }

  # Flatten lists and clean up
  if (run == "all") {
    individual.colours.overlap.flattened <- flatten.list(individual.colours.overlap)
    individual.colours.overlap.summary <- lapply(individual.colours.overlap.flattened, function (list) {list$"Diagnosable from each other"})
    individual.colours.overlap.summary.df <- data.frame(do.call(cbind, individual.colours.overlap.summary))
    colnames(individual.colours.overlap.summary.df) <- names(individual.colours.overlap.flattened)
  }

  overlap.pca.lda.flattened <- flatten.list (overlap.pca.lda)
  names(overlap.pca.lda.flattened) <- gsub ("pca.lda.values.", "", names(overlap.pca.lda.flattened)) # sheet names cannot exceed 31
  pairwise.names <- data.frame (overlap.pca.lda.flattened[[1]][, c(1,2)])
  overlap.pca.lda.summary <- lapply(overlap.pca.lda.flattened, function (list) {list$"Diagnosable from each other"})
  overlap.pca.lda.summary.df <- data.frame(do.call(cbind, overlap.pca.lda.summary))
  colnames(overlap.pca.lda.summary.df) <- names(overlap.pca.lda.flattened)

  if (run == "all") {
    all.results.summary.df <- data.frame (cbind (pairwise.names, individual.colours.overlap.summary.df, overlap.pca.lda.summary.df))
    all.results <- c( all.results.summary = list(all.results.summary.df), individual.colours.overlap.flattened, overlap.pca.lda.flattened)
  }

  if (run == "key") {
    all.results.summary.df <- data.frame (cbind (pairwise.names, overlap.pca.lda.summary.df))
    all.results <- c( overlap.pca.lda.summary = list(all.results.summary.df), overlap.pca.lda.flattened)
  }

  if (export.excel) {
    if (is.null(excel.file.name)) {excel.file.name <- "diagnosability"}
    write.xlsx(all.results, paste (output.folder, "\\", excel.file.name, ".xlsx", sep = "") )
  }

  # Run NJ Phylogram

  neighbour.joining(colour.data, accession.number, population, R.col, G.col, B.col,
                    export.figures = export.figures, output.folder = output.folder, output.file.prefix = NULL)

  return( all.results )

}


# Conduct key analyses from Teo et al. (2025) all at once: 1) Amadon test after conducting LDA on RGB values, 2) Neighbour-joining trees after calculating pairwise CIEDE2000 on RGB values (NOT USED ANYMORE, COMBINED INTO ONE SCRIPT run.functions()) ####
# run.key.functions <- function(colour.data, accession.number, population, R.col, G.col, B.col,
#                               threshold.1 = 0.75, threshold.2 = 0.99,
#                               export.figures = FALSE, export.excel = FALSE,
#                               output.folder = NULL, excel.file.name = NULL) {
#
#   # Create output file path if needed
#   if (export.figures | export.excel) {
#     if (is.null(output.folder)) { output.folder <- getwd() } # if no output.folder designated, create files in working directory
#     if (!dir.exists(output.folder)) {dir.create(output.folder, recursive = TRUE)}
#   }
#
#   # Run PCA and LDA check
#   if (export.figures) {
#     overlap.pca.lda <- pca.lda.overlap(colour.data, accession.number, population,
#                                        R.col, G.col, B.col,
#                                        threshold.1 = threshold.1, threshold.2 = threshold.2,
#                                        export.figures = export.figures, output.folder = output.folder) }
#   else { overlap.pca.lda <- pca.lda.overlap(colour.data, accession.number, population,
#                                             R.col, G.col, B.col,
#                                             threshold.1 = threshold.1, threshold.2 = threshold.2) }
#
#   # Flatten lists and clean up
#   overlap.pca.lda.flattened <- flatten.list (overlap.pca.lda)
#   names(overlap.pca.lda.flattened) <- gsub ("pca.lda.values.", "", names(overlap.pca.lda.flattened)) # sheet names cannot exceed 31
#   pairwise.names <- data.frame (overlap.pca.lda.flattened[[1]][, c(1,2)])
#   overlap.pca.lda.summary <- lapply(overlap.pca.lda.flattened, function (list) {list$"Diagnosable from each other"})
#   overlap.pca.lda.summary.df <- data.frame(do.call(cbind, overlap.pca.lda.summary))
#   colnames(overlap.pca.lda.summary.df) <- names(overlap.pca.lda.flattened)
#
#   all.results.summary.df <- data.frame (cbind (pairwise.names, overlap.pca.lda.summary.df))
#
#   all.results <- c( overlap.pca.lda.summary = list(all.results.summary.df), overlap.pca.lda.flattened)
#
#   if (export.excel) {
#     if (is.null(excel.file.name)) {excel.file.name <- "diagnosability"}
#     write.xlsx(all.results, paste (output.folder, "\\", excel.file.name, ".xlsx", sep = "") )
#   }
#
#   # Run NJ Phylogram
#   neighbour.joining (colour.data, accession.number, population, R.col, G.col, B.col,
#                      export.figures = export.figures, output.folder = output.folder, output.file.prefix = NULL)
#
#
#   return( all.results )
# }
#
#
# # Conduct all analyses described in Teo et al. (2025) including exploratory analyses (NOT USED ANYMORE, COMBINED INTO ONE SCRIPT run.functions())
# run.all.functions <- function(colour.data, accession.number, population, R.col, G.col, B.col,
#                               threshold.1 = 0.75, threshold.2 = 0.99,
#                               export.colours.plot = FALSE, export.figures = FALSE, export.excel = FALSE,
#                               output.folder = NULL, excel.file.name = NULL) {
#
#   # Function to calculate overlap for individual colours based on all colour models
#   calculate.overlap.for.all.colours <- function(colour.data, colour.data.from.all.colour.models,
#                                                 accession.number, population, # need colour.data for metadata
#                                                 threshold.1 = 0.75, threshold.2 = 0.99) {
#
#     colour.model.overlap.check <- vector (length = length (colour.data.from.all.colour.models), mode = "list") # create list to store colour diagnosability
#     names(colour.model.overlap.check) <- names(colour.data.from.all.colour.models)
#
#     for (colour.model in 1:length(colour.data.from.all.colour.models)) {
#
#       colour.model.used <- colour.data.from.all.colour.models[[colour.model]]
#       colour.model.used.name <- names (colour.data.from.all.colour.models) [colour.model]
#       colour.model.diagnose <- vector (length = ncol (colour.model.used), mode = "list")
#       names(colour.model.diagnose) <- colnames(colour.model.used)
#
#       for (colour.parameter in 1:ncol(colour.model.used)) {
#         column.name <- colnames(colour.model.used)[colour.parameter]
#         dataframe.to.run.pairwise <- cbind (colour.data [, c(accession.number, population)], colour.model.used)
#         colour.model.diagnose[[colour.parameter]] <- pairwise.group.Amadon (dataframe.to.run.pairwise,
#                                                                               accession.number,
#                                                                               population,
#                                                                               column.name, threshold.1, threshold.2)
#       }
#       colour.model.overlap.check[[colour.model]] <- colour.model.diagnose
#     }
#
#     return(colour.model.overlap.check)
#   }
#
#   # Function to plot RGB and CYM colours
#   plot.colours <- function(colour.data, colour.data.from.all.colour.models, population, output.folder) {  # need colour.data for metadata also
#
#     RGB.data <- colour.data.from.all.colour.models[["rgb"]]
#     CYM.data <- colour.data.from.all.colour.models[["cmy"]]
#
#     RGB.data.long <- data.frame(
#       subspecies = rep(colour.data[,population], times = 3),
#       value = c(RGB.data$R, RGB.data$G, RGB.data$B),
#       colour = rep(c("R", "G", "B"), each = nrow(colour.data)))
#
#     CYM.data.long <- data.frame(
#       subspecies = rep(colour.data[,population], times = 3),
#       value = c(CYM.data$C, CYM.data$Y, CYM.data$M),
#       colour = rep(c("C", "Y", "M"), each = nrow(colour.data)))
#
#     RGB.plot <- ggplot(RGB.data.long, aes(x = factor(colour, levels = c("R", "G", "B")), y = value, fill = colour)) +
#       geom_boxplot() +
#       facet_wrap(~subspecies) +  # Separate plots for each group
#       labs(title = paste("Boxplot of CYM values by ",population, sep = ""), x = "Colour", y = "Value") +
#       scale_fill_manual(values = c("R" = "red", "G" = "green", "B" = "blue")) +  # Set colors for R, G, B
#       theme_minimal()
#
#     CYM.plot <- ggplot(CYM.data.long, aes(x = factor(colour, levels = c("C", "Y", "M")), y = value, fill = colour)) +
#       geom_boxplot() +
#       facet_wrap(~subspecies) +  # Separate plots for each group
#       labs(title = paste("Boxplot of CYM values by ",population, sep = ""), x = "Colour", y = "Value") +
#       scale_fill_manual(values = c("C" = "cyan", "Y" = "yellow", "M" = "magenta")) +  # Set colors for R, G, B
#       theme_minimal()
#
#     png(filename = paste (output.folder,"/box_colours_RGB_",
#                           population,
#                           ".png", sep=""), width = 1000, height = 1000)
#
#     print (RGB.plot)
#
#     dev.off()
#
#     png(filename = paste (output.folder,"/box_colours_CYM_",
#                           population,
#                           ".png", sep=""), width = 1000, height = 1000)
#
#     print (CYM.plot)
#
#     dev.off()
#   }
#
#   # Create output file path if needed
#   if (export.figures | export.colours.plot | export.excel ) {
#     if (is.null(output.folder)) { output.folder <- getwd() } # if no output.folder designated, create files in working directory
#     if (!dir.exists(output.folder)) {dir.create(output.folder, recursive = TRUE)}
#   }
#
#   colour.data.from.all.colour.models <- add.colour.models.from.RGB(colour.data, R.col, G.col, B.col)
#
#
#   # Plot colours if needed
#   if (export.colours.plot) { plot.colours(colour.data, colour.data.from.all.colour.models, population, output.folder) }
#
#   # Run individual colours overlap check
#   overlap.colours <- calculate.overlap.for.all.colours (colour.data, colour.data.from.all.colour.models, accession.number, population, threshold.1 = threshold.1, threshold.2 = threshold.2)
#
#   # Run PCA and LDA check
#   if (export.figures) {
#     overlap.pca.lda <- pca.lda.overlap(colour.data, accession.number, population,
#                                        R.col, G.col, B.col,
#                                        threshold.1 = threshold.1, threshold.2 = threshold.2,
#                                        colour.data.from.all.colour.models = colour.data.from.all.colour.models,
#                                        export.figures = export.figures, output.folder = output.folder) }
#   else { overlap.pca.lda <- pca.lda.overlap(colour.data, accession.number, population,
#                                             R.col, G.col, B.col,
#                                             colour.data.from.all.colour.models = colour.data.from.all.colour.models,
#                                             threshold.1 = threshold.1, threshold.2 = threshold.2) }
#
#   # Calculate overlap for individual colours
#   individual.colours.overlap <- calculate.overlap.for.all.colours (colour.data, colour.data.from.all.colour.models,
#                                                                    accession.number, population,
#                                                                    threshold.1 = threshold.1, threshold.2 = threshold.2)
#
#   # Flatten lists and clean up
#   individual.colours.overlap.flattened <- flatten.list(individual.colours.overlap)
#   individual.colours.overlap.summary <- lapply(individual.colours.overlap.flattened, function (list) {list$"Diagnosable from each other"})
#   individual.colours.overlap.summary.df <- data.frame(do.call(cbind, individual.colours.overlap.summary))
#   colnames(individual.colours.overlap.summary.df) <- names(individual.colours.overlap.flattened)
#
#   overlap.pca.lda.flattened <- flatten.list (overlap.pca.lda)
#   names(overlap.pca.lda.flattened) <- gsub ("pca.lda.values.", "", names(overlap.pca.lda.flattened)) # sheet names cannot exceed 31
#   pairwise.names <- data.frame (overlap.pca.lda.flattened[[1]][, c(1,2)])
#   overlap.pca.lda.summary <- lapply(overlap.pca.lda.flattened, function (list) {list$"Diagnosable from each other"})
#   overlap.pca.lda.summary.df <- data.frame(do.call(cbind, overlap.pca.lda.summary))
#   colnames(overlap.pca.lda.summary.df) <- names(overlap.pca.lda.flattened)
#
#   all.results.summary.df <- data.frame (cbind (pairwise.names, individual.colours.overlap.summary.df, overlap.pca.lda.summary.df))
#
#   all.results <- c( all.results.summary = list(all.results.summary.df), individual.colours.overlap.flattened, overlap.pca.lda.flattened)
#
#   if (export.excel) {
#     if (is.null(excel.file.name)) {excel.file.name <- "diagnosability"}
#     write.xlsx(all.results, paste (output.folder, "\\", excel.file.name, ".xlsx", sep = "") )
#   }
#
#   # Run NJ Phylogram
#
#   neighbour.joining(colour.data, accession.number, population, R.col, G.col, B.col,
#                     export.figures = export.figures, output.folder = output.folder, output.file.prefix = NULL)
#
#   return( all.results )
#
# }

