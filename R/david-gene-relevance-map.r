library(FNN)
library(BiocParallel)
library(ggplot2)

#' Compute gene relevances for entire data set
#' 
#' Parallelised.
#' 
#' @param coord: (matrix cells x number of dimensions in
#' cell similarity space)
#' Coordinates of each cell in given space.
#' Not required if dm is given.
#' Can be used to feed in non-diffusion map projections.
#' @param matExpression: (matrix cells x genes)
#' Expression values for each cellxgenes.
#' Not required if dm is given.
#' Can be used to feed in non-diffusion map projections.
#' @param dm: (diffusion map)
#' Object based on which gene relevance is computed.
#' Not required if coord and matExpression given.
#' @param coorddim: (str vector)
#' Columns (dimensions) of projection space to use.
#' @paramNProc: (scalar)
#' Number of threads for parallelisation.
#' 
#' @return Either
#' dm (Diffusion map) with relevanceMap added
#' or relevanceMap as list if dm was not supplied.
#' @export
david_compute_gene_relevance <- function(coord=NULL,
																	 matExpression=NULL,
																	 dm=NULL,
																	 k=20,
																	 coorddim=c("DC1","DC2"),
																	 NProc=1){
	
	# Set parallelisation
	if(NProc > 1){ register(MulticoreParam(workers=NProc)) 
	} else { register(SerialParam()) }
	
	# Select input: dm object if called post-wrapper
	# or individual data objects if called from wrapper.
	# Extract accordingly.
	if(!is.null(dm)){
		matExpression <- t(dm@data_env$data@assayData$exprs) # cells x genes
		coord <- dm@eigenvectors[,coorddim]
	} else {
		coord <- coord[,coorddim]
	}
	nn.index <- get.knn(matExpression, k, algorithm = 'cover_tree')$nn.index
	
	scaK <- dim(nn.index)[2]
	scaNGenes <- dim(matExpression)[2]
	scaNCells <- dim(coord)[1]
	scaDCDims <- dim(coord)[2] # Diffusion components to compute partial derivatives for
	arr3DPartials <- array(NA, dim=c(scaNGenes, scaNCells, scaDCDims))
	
	for(dc in seq(1,scaDCDims)){
		# Compute partial derivatives in direction of current 
		# dimension (e.g. current DC)
		# matPartialDC: genes x cells
		
		# Pre-compute change in value of function: DC coordinate
		# !Normalised by range of dc
		#scaRangeDC <- max(coord[,dc]) - min(coord[,dc])
		matDiffCoord <- apply(nn.index, 2, function(nn){
			(coord[nn,dc]-coord[,dc])#/scaRangeDC
		})
		
		# Parallise over genes:
		arr3DPartials[,,dc] <- do.call(rbind, bplapply(seq(1,scaNGenes), function(gene){
			# Compute change in expression
			# Do not compute if reference is zero, could be drop-out
			vecMaskedImputed <- matExpression[,gene]
			vecMaskedImputed[vecMaskedImputed==0] <- NA
			matDiffExpr <- apply(nn.index, 2, function(nn){
				matExpression[nn,gene]-vecMaskedImputed
			})
			matDiffExpr[matDiffExpr==0] <- NA # Cannot evaluate partial
			
			# Compute median of difference quotients to NN
			vecCellPartials <- apply( matDiffCoord/matDiffExpr, 1, function(cell){
				# Only compute gradient if at least two observations are present!
				# Otherwise very unstable
				if(sum(!is.na(cell)) > 1)	return(median(cell, na.rm=TRUE))
				else return(NA)
			})
			return(vecCellPartials)
		}))
	}
	
	# Compute norm over partial derivates: Frobenius
	matNormPartials <- apply(arr3DPartials, c(1,2), function(z) sqrt(sum( z^2, na.rm=TRUE)) )
	
	# Find outlier cells: Not in NN of more than 1 other cell
	# Remove these as they tend to receive very large norms
	#vecboolOutlier <- sapply(seq(1,scaNCells), function(cell){
	#  sum(nn.index==cell) > 1
	#})
	#matNormPartials[,vecboolOutlier] <- NA
	
	# Prepare output
	rownames(matNormPartials) <- colnames(matExpression)
	dimnames(arr3DPartials)[[1]] <- colnames(matExpression)
	dimnames(arr3DPartials)[[3]] <- coorddim
	relevanceMap <- list(matNormPartials=matNormPartials,
											 arr3DPartials=arr3DPartials,
											 nn.index=nn.index)
	if(is.null(dm)){
		return(relevanceMap)
	} else {
		dm@relevanceMap <- relevanceMap
		return(dm)
	}
}


#' Plot gradient maps for selected genes
#' 
#' Based on diffusion map object with fitted gene relevances.
#' 
#' @param dm: (diffusion map)
#' Object with fitted gene relevances.
#' @param geneIDs: (str vector)
#' Genes to create gradient maps forl
#' 
#' @return list: List of ggplot2 objects.
#' One object (plot) per gene, ie one gradient
#' map per gene given.
#' @export
david_plot_gradient_maps <- function(dm,
															 geneIDs){
	
	counts <- dm@data_env$data@assayData$exprs
	lsGplotsGradientMaps <- list()
	for(strGene in geneIDs){
		vecExprValue <- dm@data_env$data@assayData$exprs[strGene,]
		vecGenePartialNorm <- dm@relevanceMap$matNormPartials[strGene,]
		# Select highest vectors in neighbourhoods
		vecboolTopNorm <- apply(cbind(seq(1, dim(dm@relevanceMap$nn.index)[1]) , dm@relevanceMap$nn.index), 1, function(cell) which.max(vecGenePartialNorm[cell])==1 )
		vecboolTopNorm[sapply(vecboolTopNorm, function(cell) length(cell)==0 )] <- FALSE
		vecboolTopNorm <- unlist(vecboolTopNorm)
		
		vecDC1Coord <- dm$DC1[vecboolTopNorm]
		vecDC2Coord <- dm$DC2[vecboolTopNorm]
		vecPartialDC1Gene <- dm@relevanceMap$arr3DPartials[strGene,vecboolTopNorm,1]
		vecPartialDC2Gene <- dm@relevanceMap$arr3DPartials[strGene,vecboolTopNorm,2]
		# Scale magnitude of partial derivates
		scaDeltaDC1 <- max(vecDC1Coord, na.rm=TRUE)-min(vecDC1Coord, na.rm=TRUE)
		scaDeltaDC2 <- max(vecDC2Coord, na.rm=TRUE)-min(vecDC2Coord, na.rm=TRUE)
		scaVectorLength <- 0.05 # Fraction of overall DC variability
		vecPartialDC1Gene <- vecPartialDC1Gene/max(abs(vecPartialDC1Gene), na.rm=TRUE)*scaVectorLength*scaDeltaDC1
		vecPartialDC2Gene <- vecPartialDC2Gene/max(abs(vecPartialDC2Gene), na.rm=TRUE)*scaVectorLength*scaDeltaDC2
		
		# Plot gradient vectors
		dfScatter <- data.frame(
			DC1=dm$DC1, 
			DC2=dm$DC2, 
			Expression=vecExprValue
		)
		dfGradientVectorsPositive <- data.frame(x1 = vecDC1Coord, 
																						y1 = vecDC2Coord, 
																						x2 = vecDC1Coord + vecPartialDC1Gene, 
																						y2 = vecDC2Coord + vecPartialDC2Gene,
																						partialNorm=vecGenePartialNorm[vecboolTopNorm])
		dfGradientVectorsNegative <- data.frame(x1 = vecDC1Coord, 
																						y1 = vecDC2Coord, 
																						x2 = vecDC1Coord - vecPartialDC1Gene, 
																						y2 = vecDC2Coord - vecPartialDC2Gene,
																						partialNorm=vecGenePartialNorm[vecboolTopNorm])
		gplotGradientMap <- ggplot() +
			geom_point(data=dfScatter, aes(x=DC1, y=DC2, colour=Expression), alpha=1) + 
			scale_colour_gradient(low="white", high="red") + 
			geom_segment(data = dfGradientVectorsPositive, 
									 aes(x = x1, y = y1, xend = x2, yend = y2, alpha=partialNorm)) +
			geom_segment(data = dfGradientVectorsNegative, 
									 aes(x = x1, y = y1, xend = x2, yend = y2, alpha=partialNorm)) +
			labs(title=strGene)
		lsGplotsGradientMaps[[length(lsGplotsGradientMaps)+1]] <- gplotGradientMap
	}
	return(lsGplotsGradientMaps)
}

#' Plot gene relevance map
#' 
#' Either based on global gradient norm maxima or for givne list of gene IDs.
#' 
#' @param matNormPartials: (matrix cells x genes)
#' Norm of partial derivatives wrt to considered dimensions
#' in reduced space.
#' Not required if dm is given.
#' Can be used to feed in non-diffusion map projections.
#' @param coord: (matrix cells x number of dimensions in
#' cell similarity space)
#' Coordinates of each cell in given space.
#' Not required if dm is given.
#' Can be used to feed in non-diffusion map projections.
#' @param dm: (diffusion map)
#' Object with fitted gene relevances.
#' @param smoothingiter: (scalar)
#' Number of label smoothing iterations to perform
#' on relevance map - the higher the more homogenous and 
#' the less local structure.
#' @param ids: (str vector)
#' Genes to based relevance map on.
#' @param nids: (scalar)
#' Number of genes to use if ids was not supplied.
#' 
#' @return list:
#' relMap: ggplot2 plot, the relevance map.
#' relMapData: Data frame for relMap to explore other dimension combinations.
#' ids: IDs used.
#' @export
david_plot_gene_relevance_map <- function(matNormPartials=NULL,
																		coord=NULL,
																		dm=NULL,
																		smoothingiter=2,
																		ids=NULL,
																		nids=5){
	
	if(!is.null(dm)){
		matNormPartials <- dm@relevanceMap$matNormPartials
		coord <- dm@eigenvectors[,dimnames(dm@relevanceMap$arr3DPartials)[[3]]]
	}
	# Select most often occurring genes with maximal norm of gradients at a cell.
	if(is.null(ids)){
		vecGeneIDsMaxNormPerCell <- rownames(matNormPartials)[unlist(apply(matNormPartials, 2, function(cell) which.max(cell) ))]
		vecHistGeneMaxNormPerCell <- sapply(unique(vecGeneIDsMaxNormPerCell), function(id) sum(vecGeneIDsMaxNormPerCell==id,na.rm=TRUE) )
		vecGeneIDsMap <- names(vecHistGeneMaxNormPerCell)[sort(vecHistGeneMaxNormPerCell, decreasing=TRUE, index.return=TRUE)$ix[1:nids]]
	} else { 
		vecboolFound <- sapply(ids, function(id) length(grep(id, rownames(matNormPartials))) > 0 )
		vecGeneIDsMap <- ids[vecboolFound]
	}
	# Plot a single map with cells coloured by gene which has 
	# the highest gradient norm of all genes considered.
	
	# Add more than two DC and return data frame so that user
	# can easily rebuild relevance map on other DC combination than 1 and 2.
	dfRelevanceMap <- data.frame(
		DC1=dm$DC1,
		DC2=dm$DC2,
		DC3=dm$DC3,
		DC4=dm$DC4,
		DC5=dm$DC5,
		Gene=NA
	)
	lsidxMaxGene <- apply(matNormPartials[vecGeneIDsMap,], 2, function(cell) which.max(cell) )
	lsidxMaxGene[sapply(lsidxMaxGene, function(cell) length(cell)==0 )] <- NA
	vecMaxGene <- vecGeneIDsMap[unlist(lsidxMaxGene)]
	# Label smoothing through graph structure
	if(smoothingiter>0){
		vecMaxGeneSmoothed <- vecMaxGene
		for(i in seq(1,smoothingiter)){
			vecMaxGeneSmoothed <- apply(dm@relevanceMap$nn.index, 1, function(cell){
				vecMaxGenesInNN <- unique(vecMaxGeneSmoothed[cell])
				vecHistMaxGenesInNN <- sapply(vecMaxGenesInNN, function(gene) sum(gene==vecMaxGeneSmoothed[cell], na.rm=TRUE))
				strMaxGeneInNN <- names(vecHistMaxGenesInNN)[which.max(vecHistMaxGenesInNN)]
				return(strMaxGeneInNN)
			})
		}
	}
	dfRelevanceMap$Gene <- as.factor(vecMaxGeneSmoothed)
	
	gplotRelevanceMap <- ggplot() +
		geom_point(data=dfRelevanceMap, aes(x=DC1, y=DC2, colour=Gene), alpha=0.8) + 
		labs(title="Gene relevance map: DC1 vs DC2")
	
	return(list(relMap=gplotRelevanceMap,
							relMapData=dfRelevanceMap,
							ids=vecGeneIDsMap))
}
