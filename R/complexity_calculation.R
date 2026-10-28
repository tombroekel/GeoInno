#' Import needed functions
#'
#'
#' Artificial data set to illustrate the calculation of the technological complexity of patents
#' A data set illustrating the form and style how patent data should be structured to be used by the function \emph{structural_complexity}.
#' @format ## `pat_df`
#' A data frame with 1,000 rows and 4 columns:
#' \describe{
#'   \item{appln_id}{application number of patent.}
#'   \item{cpc}{IPC/CPC code assigned to a patent.}
#'   \item{tech}{Technology code, usually 2 or 4-digit CPC code.}
#'   \item{year}{Year of observation.}
#' }
#' @source None
"pat_df"


#' Artificial data set to illustrate the calculation of the complexity frontier.
#' A data set illustrating the form and style how the data should be structured to be used by the function \emph{complexity_frontier}.
#'
#' @format ## `complexity_dat`
#' A data frame with 995 rows and 7 columns:
#' \describe{
#'   \item{geo}{The ids of the spatial units, e.g. regions.}
#'   \item{year}{The year for which the calculation is done.}
#'   \item{appln_id}{Unique patent document identifier.}
#'   \item{tech}{Technology identifier with less digits than patent class.}
#'   \item{structural}{Random number between 1 and 12 representing the value of structural complexity.}
#' }
#' @source None
"complexity_dat"

#'
#' Giant component
#'
#' @description Helper function for the calculation of the network diversity score of \insertCite{Emmert-Streib2012;textual}{GeoInno} and the structural diversity complexity measure of \insertCite{Broekel2019;textual}{GeoInno}.
#' @param graph An igraph object from which the giant component needs to be extracted.
#' @return The function returns an igraph object representing the giant component of a network
#' @export
#'
#' @examples
#' my.graph <- igraph::random.graph.game(p.or.m = 1/10, n=10)
#' giant.component(my.graph)
giant.component <- function(graph)
{
  cl <- igraph::clusters(graph)
  sub.g <- igraph::induced_subgraph(graph, which(cl$membership == which.max(cl$csize)))
  return(sub.g)
}

#' individual network diversity score (iNDS)
#'
#' @md
#' @description Historical fixed-step individual NDS calculation retained for reproducibility. The walk uses `reps_i` steps, resets the seed internally, and converts undefined sample scores to zero. Used by [structural_diversity_legacy()].
#' @importFrom Rdpack reprompt
#' @import tidyverse
#' @import widyr
#' @import future.apply
#' @import progressr
#' @import future
#' @import data.table
#' @importFrom igraph random_walk induced.subgraph cluster_walktrap vcount ecount induced_subgraph clusters random.graph.game graph_from_edgelist simplify sizes
#' @importFrom stats median quantile var window
#' @param g The igraph object. Usually, a binarized version of the combinatorial network of CPC classes co-occurring on patents.
#' @param node_x The randomly sampled node (position index in igraph object) for which the partial network is to be extracted by means of a random walk.
#' @param reps_i Number of random-walk steps, default 200; not a bootstrap repetition count.
#'
#' @return Returns the value of the iNDS measure.
#' @export
#'
#' @examples
#' my.graph <- igraph::random.graph.game(p.or.m = 1/10, n=10)
#' NDS.intern(node_x=1, g = my.graph, reps_i=10)
NDS.intern<-function(node_x=x, g, reps_i)
  {
  x<-NULL
  set.seed(123)
  if(missing(reps_i)){reps_i=200}
  sample_vertex <- igraph::random_walk(graph=g, start=node_x, steps = reps_i, mode="all")
  sample_net <- igraph::induced.subgraph(graph = g, vids = sample_vertex)
  set.seed(123)
  modules_g <- igraph::cluster_walktrap(sample_net, steps = 4, merges = TRUE, modularity = TRUE, membership = TRUE)
  m <- data.frame(igraph::sizes(modules_g))
  lap <- eigen(igraph::laplacian_matrix(sample_net, normalized=F))$values
  graph.3 <-  igraph::count_motifs(sample_net, 3)
  graph.4 <-  igraph::count_motifs(sample_net, 4)
  a.module <- length(modules_g)/igraph::vcount(sample_net)
  v.module <- var(m$Freq)/mean(m$Freq)
  v.l <- var(lap)/mean(lap)
  r.motif <- graph.3/graph.4
  nds.s <- (a.module*r.motif)/(v.module*v.l)
  nds.s <- ifelse(is.infinite(nds.s)==T | is.na(nds.s) ==T, 0, nds.s)
  return(nds.s)
}

#' Function nds
#' @md
#' @description Historical fixed-step NDS aggregator retained for reproducibility. The network must be binary and connected (largest component). The function samples up to `s` starting vertices, calls [NDS.intern()], substitutes 1 for a zero aggregate, caps values above 1, and returns `-log()` of the result.
#' @param g An igraph network main component.
#' @param reps_i Number of random-walk steps in each sampled subnetwork, default 200; not a bootstrap repetition count.
#' @param s Number of starting vertices (sampled subnetworks), default 125.
#'
#' @return A numeric value of the NDS-score, after a log-transformation and multiplication with -1.
#' @export
#'
#' @examples
#' my.graph <- igraph::random.graph.game(p.or.m = 1/10, n=10)
#' nds(g = my.graph, s = 125, reps_i = 10)
nds<-function(g,s,reps_i)
  {
  if(missing(s)){s=125}
  if(missing(reps_i)){reps_i=200}
  vsample<-ifelse(igraph::vcount(g) < s, igraph::vcount(g), s)
  set.seed(123)
  start_vertex<-sample(igraph::vcount(g), vsample, replace=F)
  comp<-sort(unlist(lapply(start_vertex, NDS.intern, g, reps_i)))
  structural <- mean(comp, na.rm=T)
  structural <- ifelse(structural == 0, 1, structural)
  structural <- ifelse(structural > 1, 1, structural)
  structural <- log(structural) * -1
  return(structural)
  }

#' Historical structural diversity (legacy fixed-step implementation)
#'
#' @md
#' @description This is the pre-2.0 GeoInno implementation, retained to reproduce earlier analyses. It uses fixed-length random walks, internal seed resets, historical treatment of undefined iNDS, and historical capping/transformation. Its numerical algorithm is unchanged. Patent input has `appln_id`, `year`, `tech`, and `cpc` columns.
#' @param p.dat A data frame in the form of \emph{pat_df}.
#' @param mw Parameter setting the length of moving window. Default set to 3 implying that patent data will be pooled across three years (t, t+1, t+2).
#' @param node.sample Number of starting vertices/sample subnetworks, default 125.
#' @param reps Number of random-walk steps per sampled subnetwork, default 200; not a bootstrap repetition count.
#' @param core.workers The number of cores used in parallelization. Default set to 1.
#'
#' @return A data frame with `year`, `tech`, annual and moving-window patent/CPC
#'   counts, historical transformed `structural` score, and largest-component
#'   `edges` and `nodes` counts. The historical score replaces non-finite
#'   sample-level iNDS with zero, substitutes 1 for a zero aggregate, caps raw
#'   values above 1 at 1, and then applies `-log()`.
#' @export
#'
#'@references
#' \insertAllCited{}
#'
#'
#' @examples
#' structural_diversity_legacy(pat_df)
structural_diversity_legacy <-function (p.dat = pat_df, mw = 3, node.sample = 125, reps = 200,
                                  core.workers = 1)
{
  pat_df <- NULL
  year <- NULL
  tech <- NULL
  appln_id <- NULL
  geo <- NULL
  structural <- NULL
  item1 <- NULL
  item2 <- NULL
  id <- NULL
  cpc <- NULL
  p.dat <- p.dat %>% dplyr::rowwise() %>% dplyr::mutate(window = paste(rep(year, mw) + c(0:(mw - 1)), collapse = "_"))
  p.dat <- p.dat %>% tidyr::separate_rows(window, sep = "_")
  p.dat <- p.dat %>% dplyr::mutate(window = as.numeric(window))

  results_year <- p.dat %>% dplyr::group_by(year, tech) %>%
                            dplyr::summarise(patents_year = dplyr::n_distinct(appln_id),
                            cpcs_year = dplyr::n_distinct(cpc)) %>% dplyr::ungroup()

  results_window <- p.dat %>% dplyr::group_by(window, tech) %>%
                              dplyr::summarise(patents_window = dplyr::n_distinct(appln_id),
                              cpcs_window = dplyr::n_distinct(cpc)) %>% dplyr::ungroup()

  progressr::handlers("progress")
  p.dat <- p.dat %>% dplyr::mutate(tech_pats = paste(p.dat$tech,
                                                     p.dat$window, sep = "_")) %>% dplyr::distinct()

  split_data <- base::split(x = p.dat, f = p.dat$tech_pats)

  future::plan(multisession, workers = core.workers)
  with_progress({
    p <- progressr::progressor(steps = length(split_data))
    dfs <- future.apply::future_lapply(split_data, future.seed = NULL,
                                       FUN = function(x, future.label = T, ...) {
                                         p()
                                         Sys.sleep(0.2)
                                         pats_window <- x %>% dplyr::pull(appln_id) %>%
                                                              unique()

                                         if (length(pats_window) > 0) {
                                           pat_tech <- p.dat %>% dplyr::filter(appln_id %in%
                                                                                 pats_window) %>% widyr::pairwise_count(item = cpc,
                                                                                                                        feature = appln_id, diag = F)
                                           pat_net <- pat_tech %>% dplyr::select(item1,
                                                                                 item2) %>% as.matrix() %>% igraph::graph_from_edgelist(directed = FALSE) %>%
                                             igraph::simplify()
                                           if (igraph::ecount(pat_net) > 0) {
                                             g_inv <- giant.component(pat_net)
                                             node_count <- igraph::vcount(g_inv)
                                             edge_count <- igraph::ecount(g_inv)
                                             if (edge_count >= 1) {
                                               tibble::tibble(structural = nds(g = g_inv,
                                                                               s = node.sample,
                                                                               reps_i = reps),
                                                                               edges = edge_count,
                                                                               nodes = node_count)
                                             }
                                             else {
                                               tibble::tibble(structural = NA, edges = edge_count,
                                                              nodes = node_count)
                                             }
                                           }
                                           else {
                                             tibble::tibble(structural = NA, edges = NA,
                                                            nodes = NA)
                                           }
                                         }
                                         else {
                                           tibble::tibble(structural = NA, edges = NA,
                                                          nodes = NA)
                                         }
                                       })
  })
  dfs_df <- data.table::rbindlist(dfs, fill = T)
  dfs_df <- dplyr::bind_cols(id = names(dfs), dfs_df)
  dfs_df <- dfs_df %>% tidyr::separate(id, sep = "_", into = c("tech",
                                                               "window")) %>% dplyr::mutate(window = as.numeric(window))
  results_window <- dplyr::left_join(results_window, dfs_df,
                                     by = c("tech", "window"), na_matches = "never")
  results_window <- dplyr::left_join(results_year, results_window,
                                     by = c("tech", year = "window"), na_matches = "never")
  return(results_window)
}

#' Structural diversity with fixed-unique-node sampling
#'
#' @md
#' @description Calculates Structural Diversity/NDS from patent co-occurrence networks.
#' For each technology-year, the largest connected component is sampled from up to
#' `n_starts` starting vertices. Each walk stops once `n_nodes` distinct vertices
#' (including the start) have been visited. If the component has at most `n_nodes`
#' vertices, the whole component is used for every start. Sample-level iNDS values
#' are averaged after excluding non-finite values. This is the recommended GeoInno
#' 2.0 estimator. [structural_diversity_legacy()] retains the historical fixed-step
#' estimator for reproducibility.
#' @param p.dat Patent data with `appln_id`, `year`, `tech`, and `cpc` columns;
#'   defaults to the package example `pat_df`.
#' @param mw Moving-window length, default 3 (years `t` through `t + mw - 1`).
#' @param n_starts Maximum number of starting vertices/sample subnetworks, default 125.
#' @param n_nodes Number of distinct vertices required in each sampled subnetwork,
#'   default 150. This is not a fixed number of walk steps.
#' @param seed Top-level random seed. Each technology-year receives a reproducible
#'   independent RNG stream through `future.apply`; distinct per-start seeds are
#'   drawn from that stream. Repeated runs and worker counts produce identical
#'   results for the same seed.
#' @param core.workers Number of parallel workers, default 1.
#' @return A data frame with technology (`tech`), year, patent/CPC counts,
#'   network edge/node counts, `Structural_raw` (mean finite sample-level iNDS),
#'   `Structural_neglog` (`-log(Structural_raw)` for positive finite raw values),
#'   `structural` (an alias of `Structural_neglog`), `n_samples`,
#'   `n_valid_samples`, and `share_valid_samples`. Undefined sample iNDS values
#'   are missing, never zero-imputed. The raw statistic is not capped at 1;
#'   values above 1 yield legitimate negative `Structural_neglog` values.
#' @export
#' @examples
#' complexity <- structural_diversity(pat_df)
#' complexity$Structural_raw
#' complexity$Structural_neglog
structural_diversity <- function(p.dat = pat_df, mw = 3, n_starts = 125,
                                 n_nodes = 150, seed = 123, core.workers = 1)
{
  if (length(n_starts) != 1L || !is.finite(n_starts) || n_starts < 1 ||
      n_starts != as.integer(n_starts)) stop("n_starts must be a positive integer.")
  if (length(n_nodes) != 1L || !is.finite(n_nodes) || n_nodes < 1 ||
      n_nodes != as.integer(n_nodes)) stop("n_nodes must be a positive integer.")
  if (length(core.workers) != 1L || !is.finite(core.workers) || core.workers < 1 ||
      core.workers != as.integer(core.workers)) stop("core.workers must be a positive integer.")
  if (length(seed) != 1L || !is.finite(seed) || seed < 0 ||
      seed != as.integer(seed)) stop("seed must be a non-negative integer.")
  force(p.dat)
  pat_df <- NULL
  year <- NULL
  tech <- NULL
  appln_id <- NULL
  geo <- NULL
  structural <- NULL
  item1 <- NULL
  item2 <- NULL
  id <- NULL
  cpc <- NULL
  p.dat <- p.dat %>% dplyr::rowwise() %>% dplyr::mutate(window = paste(rep(year, mw) + c(0:(mw - 1)), collapse = "_"))
  p.dat <- p.dat %>% tidyr::separate_rows(window, sep = "_")
  p.dat <- p.dat %>% dplyr::mutate(window = as.numeric(window))

  results_year <- p.dat %>% dplyr::group_by(year, tech) %>%
    dplyr::summarise(patents_year = dplyr::n_distinct(appln_id),
                     cpcs_year = dplyr::n_distinct(cpc)) %>% dplyr::ungroup()

  results_window <- p.dat %>% dplyr::group_by(window, tech) %>%
    dplyr::summarise(patents_window = dplyr::n_distinct(appln_id),
                     cpcs_window = dplyr::n_distinct(cpc)) %>% dplyr::ungroup()

  progressr::handlers("progress")
  p.dat <- p.dat %>% dplyr::mutate(tech_pats = paste(p.dat$tech,
                                                     p.dat$window, sep = "_")) %>% dplyr::distinct()
  split_data <- base::split(x = p.dat, f = p.dat$tech_pats)

  old_plan <- future::plan()
  # Keep the small corrected helper family together when futures serialize a
  # development checkout (an installed package resolves them by namespace).
  helper_env <- new.env(parent = asNamespace("stats"))
  for (name in c("sample_fixed_nodes", "NDS.intern_fixed",
                 "aggregate_sample_nds", "nds_fixed")) {
    helper_env[[name]] <- get(name, envir = environment(structural_diversity))
    environment(helper_env[[name]]) <- helper_env
  }
  score_fixed <- helper_env$nds_fixed
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(multisession, workers = core.workers)
  with_progress({
    p <- progressr::progressor(steps = length(split_data))
    dfs <- future.apply::future_lapply(split_data, future.seed = seed,
                                       FUN = function(x, future.label = T, ...) {
                                         p()
                                         pats_window <- x %>% dplyr::pull(appln_id) %>% unique()

                                         if (length(pats_window) > 0) {
                                           pat_tech <- p.dat %>% dplyr::filter(appln_id %in%
                                                                                 pats_window) %>% widyr::pairwise_count(item = cpc,
                                                                                                                        feature = appln_id, diag = F)
                                           pat_net <- pat_tech %>% dplyr::select(item1,
                                                                                 item2) %>% as.matrix() %>% igraph::graph_from_edgelist(directed = FALSE) %>%
                                             igraph::simplify()
                                           if (igraph::ecount(pat_net) > 0) {
                                             components <- igraph::components(pat_net)
                                             g_inv <- igraph::induced_subgraph(
                                               pat_net, which(components$membership == which.max(components$csize)))
                                             node_count <- igraph::vcount(g_inv)
                                             edge_count <- igraph::ecount(g_inv)
                                             if (edge_count >= 1) {
                                               score <- score_fixed(g = g_inv, n_starts = n_starts,
                                                                    n_nodes = n_nodes)
                                               tibble::tibble(Structural_raw = score$Structural_raw,
                                                              Structural_neglog = score$Structural_neglog,
                                                              structural = score$Structural_neglog,
                                                              n_samples = score$n_samples,
                                                              n_valid_samples = score$n_valid_samples,
                                                              share_valid_samples = score$share_valid_samples,
                                                              edges = edge_count, nodes = node_count)
                                             }
                                             else {
                                               tibble::tibble(Structural_raw = NA_real_, Structural_neglog = NA_real_,
                                                              structural = NA_real_, n_samples = 0L,
                                                              n_valid_samples = 0L, share_valid_samples = NA_real_,
                                                              edges = edge_count, nodes = node_count)
                                             }
                                           }
                                           else {
                                             tibble::tibble(Structural_raw = NA_real_, Structural_neglog = NA_real_,
                                                            structural = NA_real_, n_samples = 0L,
                                                            n_valid_samples = 0L, share_valid_samples = NA_real_,
                                                            edges = NA_real_, nodes = NA_real_)
                                           }
                                         }
                                         else {
                                           tibble::tibble(Structural_raw = NA_real_, Structural_neglog = NA_real_,
                                                          structural = NA_real_, n_samples = 0L,
                                                          n_valid_samples = 0L, share_valid_samples = NA_real_,
                                                          edges = NA_real_, nodes = NA_real_)
                                         }
                                       })
  })
  dfs_df <- data.table::rbindlist(dfs, fill = T)
  dfs_df <- dplyr::bind_cols(id = names(dfs), dfs_df)
  dfs_df <- dfs_df %>% tidyr::separate(id, sep = "_", into = c("tech",
                                                               "window")) %>% dplyr::mutate(window = as.numeric(window))
  results_window <- dplyr::left_join(results_window, dfs_df,
                                     by = c("tech", "window"), na_matches = "never")
  results_window <- dplyr::left_join(results_year, results_window,
                                     by = c("tech", year = "window"), na_matches = "never")
  return(results_window)
}

# Sample an induced subgraph by walking until exactly n_nodes distinct vertices
# have been visited. A connected component no larger than n_nodes is used whole.
sample_fixed_nodes <- function(g, node_x, n_nodes)
{
  if (igraph::vcount(g) <= n_nodes) return(g)
  visited <- rep(FALSE, igraph::vcount(g))
  current <- as.integer(node_x)
  visited[current] <- TRUE
  count <- 1L
  while (count < n_nodes) {
    adjacent <- as.integer(igraph::neighbors(g, current, mode = "all"))
    current <- adjacent[sample.int(length(adjacent), 1L)]
    if (!visited[current]) {
      visited[current] <- TRUE
      count <- count + 1L
    }
  }
  igraph::induced_subgraph(g, which(visited))
}

# Corrected counterpart of NDS.intern(): same iNDS formula, fixed-node sample,
# no internal seed reset, and non-finite sample values are missing.
NDS.intern_fixed <- function(node_x, g, n_nodes)
{
  sample_net <- sample_fixed_nodes(g, node_x, n_nodes)
  modules_g <- igraph::cluster_walktrap(sample_net, steps = 4, merges = TRUE,
                                       modularity = TRUE, membership = TRUE)
  m <- data.frame(igraph::sizes(modules_g))
  lap <- eigen(igraph::laplacian_matrix(sample_net, normalization = "unnormalized"))$values
  graph.3 <- igraph::count_motifs(sample_net, 3)
  graph.4 <- igraph::count_motifs(sample_net, 4)
  a.module <- length(modules_g)/igraph::vcount(sample_net)
  v.module <- var(m$Freq)/mean(m$Freq)
  v.l <- var(lap)/mean(lap)
  r.motif <- graph.3/graph.4
  nds.s <- (a.module*r.motif)/(v.module*v.l)
  if (!is.finite(nds.s)) return(NA_real_)
  nds.s
}

# Corrected counterpart of nds(): mean finite sample-level iNDS with diagnostics.
aggregate_sample_nds <- function(sample_scores)
{
  n_samples <- length(sample_scores)
  valid <- is.finite(sample_scores)
  n_valid_samples <- sum(valid)
  if (!any(valid)) {
    raw_nds <- NA_real_
  } else {
    raw_nds <- mean(sample_scores[valid])
  }
  transformed <- if (is.finite(raw_nds) && raw_nds > 0) -log(raw_nds) else NA_real_
  list(Structural_raw = raw_nds, Structural_neglog = transformed,
       n_samples = n_samples, n_valid_samples = n_valid_samples,
       share_valid_samples = if (n_samples > 0) n_valid_samples/n_samples else NA_real_)
}

nds_fixed <- function(g, n_starts = 125, n_nodes = 150)
{
  n_samples <- min(igraph::vcount(g), n_starts)
  start_vertex <- sample.int(igraph::vcount(g), n_samples, replace = FALSE)
  walk_seeds <- sample.int(.Machine$integer.max, n_samples, replace = FALSE)
  sample_scores <- vapply(seq_along(start_vertex), function(i) {
    set.seed(walk_seeds[i])
    NDS.intern_fixed(start_vertex[i], g = g, n_nodes = n_nodes)
  }, numeric(1))
  c(aggregate_sample_nds(sample_scores),
    list(sample_scores = sample_scores, walk_seeds = walk_seeds))
}


#' Complexity frontier
#'
#' @description The function calculates the complexity frontier for a range of spatial units following the procedure described in \insertCite{Mewes2022;textual}{GeoInno}. That is, for each spatial unit and year, the median of the x-percentil (top) most complex knowledge components is calculated. Knowledge components are represented by distinct occurences of technologies on patents. The latter means that each \emph{appln_id} is associated to the same technology only once.
#' @param c.dat A dataframe with columns: \emph{geo}, \emph{appln_id}, \emph{tech}, \emph{year}, and \emph{structural}. Only distinct technology occurrences are to be included per \emph{appln_id}. The example data set \emph{complexity_dat} illustrates the dataframe structure.
#' @param top Percentile for which frontier is to be estimated. E.g. top=0.05 corresponds to the consideration of the 5 percent most complex knowledge components.
#' @param mw Number of years considered in the mowing window in the type of \emph{t:(t+mw-1)}.
#' @return The function returns a data.frame with the name complexity.frontier including the following information:
#' \itemize{
#'   \item{\emph{geo}} {The spatial unit for which the frontier is calculated.}
#'   \item{\emph{year}} {The year for which the frontier is calculated.}
#'   \item{\emph{vfrontier}} {The value of the complexity frontier.}
#'   }
#' @export
#'@references
#' \insertAllCited{}
#'
#' @examples
#' complexity_frontier(c.dat = complexity_dat, top = 0.05, mw = 1)
complexity_frontier <- function(c.dat = complexity_dat, top = 0.05, mw = 1)
  {
  complexity_dat <- NULL
  year <- NULL
  tech <- NULL
  appln_id <- NULL
  geo <- NULL
  structural <- NULL
  c.dat <- c.dat %>% dplyr::distinct(geo, appln_id, tech, year, structural)
  if(mw>0)
    {
    c.dat <- c.dat %>% dplyr::rowwise() %>% dplyr::mutate(window=paste(rep(year, mw)+c(0:(mw-1)),collapse="_"))
    c.dat <- c.dat %>%
      tidyr::separate_rows(window, sep="_") %>%
      dplyr::mutate(window = as.numeric(window)) %>%
      dplyr::arrange(geo, window, dplyr::desc(structural)) %>% dplyr::ungroup()
    cfrontier <- c.dat %>% dplyr::group_by(geo, window) %>% dplyr::slice(1:quantile(c(1:dplyr::n()),probs=top)) %>% dplyr::summarise(vfrontier=median(structural,na.rm=T)) %>% dplyr::ungroup()
    cfrontier <- cfrontier %>% dplyr::rename(year=window)
    }
  return(cfrontier)
}

#' Create sample data
#'
#' @param num.pat Number of patents.
#' @param geos Number of spatial units.
#'
#' @return The function returns a data.frame suited to be used with the function \emph{structural_complexity}.
#' \itemize{
#'   \item{geo}{The ids of the spatial unit, e.g. regions.}
#'   \item{year}{The year for which the calculation is done.}
#'   \item{appln_id}{Unique patent document identifier.}
#'   \item{cpc}{Patent class identifier with more digits than technology.}
#'   \item{tech}{Technology identifier with less digits than patent class.}
#'   \item{structural}{Random number between 1 and 12 representing the value of structural complexity.}
#'   }
#' @export
#'
#' @examples
#' create_sample_data(num.pat=200, geos=10)
create_sample_data <- function(num.pat=200, geos=10)
{
  set.seed(123)
  regs <- sample(paste0("GEO_",1:geos), size=num.pat*5, replace = T)
  appln_id <- sample(num.pat, size=num.pat*5, replace = T)
  cpc.1 <- sample(LETTERS[1:10], size=num.pat*5, replace = T)
  cpc.2 <- sample(LETTERS[1:25], size=num.pat*5, replace = T)
  cpc <- paste0(cpc.1,cpc.2)

  pat.df <- data.frame(geo=regs, appln_id=appln_id, cpc=cpc) %>% dplyr::arrange(appln_id)
  pat.df <- pat.df %>% dplyr::mutate(tech = substr(cpc,1,1))
  year <- sample(c(1999:2005), size=length(unique(appln_id)), replace = T)
  pat.df <- pat.df %>% dplyr::mutate(year = year[appln_id])
  pat.df <- pat.df %>% stats::na.omit()
  return(pat.df)
}

#complexity.dat <- complex.dat %>% select(geo, appln_id, tech, year) %>% distinct()
#helper <- complex.dat %>% select(tech,year) %>% distinct()  %>% mutate(structural = runif(n(),min=1,max=12))
#complexity.dat <- left_join(complexity.dat,helper,by=c("tech","year"))
#complexity.dat <- complexity.dat %>% distinct()
#save(complexity.dat, file="complexity_dat.rda")

#complex.dat %>% filter(year==1999, geo=="GEO_1") %>% arrange(desc(structural)) %>% slice(1:quantile(c(1:n()),probs=0.05)) %>% pull(structural) %>% median()



