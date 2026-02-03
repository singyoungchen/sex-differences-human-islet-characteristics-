
#' S4 object holding the results of grouping GO terms.
#'
#' @slot members Character vector of the names of the gene sets.
#' @slot gene_uni Character vector of all the genes in the gene sets. This can be either the full list of all the genes associated with the list of gene sets or a reduced list defined by the user (e.g. if only signifiant genes are used for grouping).
#' @slot gene_set List representing gene sets.
#' @slot clusters List of GOrge objects each represent a subgroup of GO terms.
#' @examples
#' # example code
#'
setClass(
      Class = 'GOrge',
      slots = c(
            # character vector of GO terms to be summarized.
            members = 'character',

            # character vector of gene names for Jaccard index calculation.
            gene_uni = 'character',

            # GO gene sets as list.
            gene_set = 'list',

            # Most frequent phrases counted from the names of GO terms (members).
            #phraseSum = 'table',

            # A list of GOrge objects that represents clusters.
            clusters = 'list',

            # list of clusters returned by .cluster_leiden.
            cluster_list = 'list',
            cluster_n = 'numeric',
            # edge_list returned by .jaccardIndex.
            edge_list = 'data.frame'
      )
)

#' Group a list of GO terms (or other sort of gene sets).
#'
#' @param gene_sets A list representing gene sets.
#' @param gene_uni An optional character vector of gene names. If supplied, only genes in gene_uni will be used for Jaccard index calculation.
#' @param edge_list An optional data.frame returned by .jaccardIndex()
#' @param recursive If true, rerun the function on the individual groups found to further divide them until either the iter_max or iter_count is reached.
#' @param iter_count Don't touch this.
#' @param iter_max The maximum number of iterations to run if recursive is true.
#' @returns a 'GOrge' object.
#' @examples
#' # example code
#' @export
setGeneric(name = "CreateGOrge",
           def = function(
            gene_sets,
            gene_uni = NULL,
            edge_list = NULL,
            recursive = T,
            iter_count = 0,
            iter_max = 2,
            verbose = T
           ){standardGeneric("CreateGOrge")},
           signature = c('gene_sets')
)
setMethod(
      f = "CreateGOrge",
      signature = signature('list'),
      definition = function(
            gene_sets,
            gene_uni = NULL,
            edge_list = NULL,
            recursive = T,
            iter_count = 0,
            iter_max = 2,
            verbose = T
      ){
            # remove genes not in gene_uni if gene_uni is supplied.
            if (!is.null(gene_uni)){
                  gene_sets |>
                        lapply(function(g){
                              g[g %in% gene_uni]
                        }) -> gene_sets
            }

            # Calculate jaccard index.
            if (is.null(edge_list)){
                  # in the case of the first iteration. calculate jaccard index.
                  edge_list <- .jaccardIndex(gene_sets)
            } else {
                  # in the case of later iteration, simply use the already calculated jaccard index.
                  edge_list <- edge_list[edge_list$V1 %in% names(gene_sets) & edge_list$V2 %in% names(gene_sets),]
            }

            # Find communities.
            .cluster_leiden(edge_list = edge_list,
                            trimEdgeWeight = 0.01,
                            resolution_parameter = 1, verbose = F) -> cluster_list

            new(Class = 'GOrge',
                members = names(gene_sets),
                gene_uni = unique(unlist(gene_sets)),
                gene_set = gene_sets,
                #phraseSum = phraseSum,
                clusters = list(),
                cluster_list = cluster_list,
                cluster_n = length(cluster_list),
                edge_list = edge_list
            ) -> OUT

            if (all(recursive, OUT@cluster_n > 1, iter_count <= iter_max)){
                  lapply(1:OUT@cluster_n, \(i){
                        if (verbose) message(paste0('Creating summaries for ', i, ' of ', OUT@cluster_n, ' clusters.'))
                        CreateGOrge(
                              gene_uni = unique(unlist(OUT@gene_set)),
                              gene_set = gene_sets[OUT@cluster_list[[i]]],
                              edge_list = OUT@edge_list,
                              iter_count = (iter_count+1),
                              iter_max = iter_max,
                              recursive = T)
                  }) -> OUT@clusters
            } else if (OUT@cluster_n > 1 & iter_count <= 1){
                  lapply(1:OUT@cluster_n, \(i){
                        if (verbose) message(paste0('Creating summaries for ', i, ' of ', OUT@cluster_n, ' clusters.'))
                        CreateGOrge(
                              gene_uni = unique(unlist(OUT@gene_set)),
                              gene_set = gene_sets[OUT@cluster_list[[i]]],
                              edge_list = OUT@edge_list,
                              iter_count = (iter_count+1),
                              recursive = F)
                  }) -> OUT@clusters
            }

            return(OUT)
      }
)


# ---- .jaccardIndex ----
# Calculate jaccard index between a list of gene sets.
.jaccardIndex <- function(
            gene_set,
            gene_uni = NULL,
            na.rm = T,
            nan.rm = T,
            inf.rm = T,
            zero.rm = T,
            verbose = T

){
      if (!is.null(gene_uni)){
            if (verbose){message("Trimming gene sets according to gene_uni.")}
            lapply(gene_set, function(g){
                  g[g %in% gene_uni]
            }) -> gene_set
      }


      combn(names(gene_set), 2) -> c

      if (verbose){
            message(
                  paste0("Calculating Jaccard index for ",
                         ncol(c),
                         " comparisons. This may take a few minutes.")
                  )
            }

      lapply(1:ncol(c), \(i){
            data.frame(
                  V1 = c[1,i],
                  V2 = c[2,i],
                  Jaccard_index =
                        length(
                              intersect(gene_set[[c[1,i]]], gene_set[[c[2,i]]])
                        )/
                        length(
                              union(gene_set[[c[1,i]]], gene_set[[c[2,i]]])
                        )
            )
      }) |> (\(x){do.call(rbind,x)})() -> OUT

      if (na.rm){
            OUT <- OUT[!is.na(OUT$Jaccard_index),]
      }
      if (nan.rm){
            OUT <- OUT[!is.nan(OUT$Jaccard_index),]
      }
      if (inf.rm){
            OUT <- OUT[!is.infinite(OUT$Jaccard_index),]
      }
      if (zero.rm){
            OUT <- OUT[OUT$Jaccard_index>0,]
      }

      return(OUT)
}

#---- .cluster_leiden ----
# clustering a list of gene_sets using igraph's cluster_leiden.
.cluster_leiden <- function(
            edge_list,
            trimEdgeWeight = NULL,
            resolution_parameter = 1,
            verbose = T
){
      # trim edges.
      if (!is.null(trimEdgeWeight)){
            if (verbose){message(paste0(
                  'Edges with weight lower than ',
                  trimEdgeWeight,
                  ' are trimmed.'
            ))}
            edge_list <- edge_list[edge_list$Jaccard_index >= trimEdgeWeight, ]
      }
      # igraph::cluster_leiden
      if(verbose){message('Detecting community using igraph::cluster_leiden().')}
      igraph::graph_from_data_frame(d = edge_list, directed = F) -> g
      igraph::cluster_leiden(graph = g,
                             objective_function = 'modularity',
                             weights = igraph::E(g)$Jaccard_index,
                             n_iterations = 5000,
                             resolution_parameter = resolution_parameter
      ) -> com
      if (verbose){
            message(paste0(com$nb_clusters, ' clusters detected.'))
      }
      # arrange communities into list.
      lapply(1:com$nb_clusters, function(x){
            com$names[com$membership == x]
      })
}





