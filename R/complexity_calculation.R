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
#' @description NDS.intern() calculates the individual network diversity score as defined by \insertCite{Emmert-Streib2012;textual}{GeoInno}. It is in the calculation of the structural diversity complexity measure of \insertCite{Broekel2019;textual}{GeoInno}.
#' @importFrom Rdpack reprompt
#' @import tidyverse
#' @import widyr
#' @import netdist
#' @import future.apply
#' @import progressr
#' @import future
#' @import data.table
#' @importFrom igraph random_walk induced.subgraph cluster_walktrap vcount ecount induced_subgraph clusters random.graph.game graph_from_edgelist simplify sizes
#' @importFrom stats median quantile var window
#' @param g The igraph object. Usually, a binarized version of the combinatorial network of CPC classes co-occurring on patents.
#' @param node_x The randomly sampled node (position index in igraph object) for which the partial network is to be extracted by means of a random walk.
#' @param reps_i The number of repetitions used in the bootstrap, default set to 200.
#'
#' @return Returns the value of the iNDS measure.
#' @export
#'
#' @examples
#' my.graph <- igraph::random.graph.game(p.or.m = 1/10, n=10)
#' NDS.intern(node_x=1, g = my.graph, reps_i=10)
NDS.intern<-function(node_x=1, g = g_sample, reps_i=200)
  {
  set.seed(123)
  sample_vertex <- igraph::random_walk(graph=g, start=node_x, steps = reps_i, mode="all")
  sample_net <- igraph::induced.subgraph(graph = g, vids = sample_vertex)
  set.seed(123)
  modules_g <- igraph::cluster_walktrap(sample_net, steps = 4, merges = TRUE, modularity = TRUE, membership = TRUE)
  m <- data.frame(igraph::sizes(modules_g))
  lap <- eigen(igraph::laplacian_matrix(sample_net, normalized=F))$values
  graphletis <- count_graphlets_for_graph(sample_net, 4) #to be replaced with igraph: count_motifs(g, 3)
  graph.3 <- sum(graphletis[paste("G",1:2,sep="")])
  graph.4 <- sum(graphletis[paste("G",3:8,sep="")])
  a.module <- length(modules_g)/vcount(sample_net)
  v.module <- var(m$Freq)/mean(m$Freq)
  v.l <- var(lap)/mean(lap)
  r.motif <- graph.3/graph.4
  nds.s <- (a.module*r.motif)/(v.module*v.l)
  nds.s <- ifelse(is.infinite(nds.s)==T | is.na(nds.s) ==T, 0, nds.s)
  return(nds.s)
}
if(getRversion() >= "2.15.1")  utils::globalVariables(c("g_sample"))

#' Function nds
#' @description Background function for the estimation of the structural diversity measure by \insertCite{Emmert-Streib2012;textual}{GeoInno}. It corresponds to the calculation of the network diversity score (NDS) of \insertCite{Broekel2019;textual}{GeoInno} for a binary network. The network must be binary and connected (giant component).
#' @param g An igraph network main component.
#' @param reps_i The number of repetition of the bootstrap.
#' @param s The number of nodes sampled in the Network Diversity Score calculation, set to 125, see \insertCite{Emmert-Streib2012;textual}{GeoInno}
#'
#' @return A numeric value of the NDS-score, after a log-transformation and multiplication with -1.
#' @export
#'
#' @examples
#' my.graph <- igraph::random.graph.game(p.or.m = 1/10, n=10)
#' nds(g = my.graph, s = 125, reps_i = 10)
nds<-function(g = g_sample, s = 125, reps_i = 200)
  {
  vsample<-ifelse(igraph::vcount(g) < s, igraph::vcount(g), s)
  set.seed(123)
  start_vertex<-sample(igraph::vcount(g), vsample, replace=F)
  comp<-sort(unlist(lapply(start_vertex, NDS.intern, g = g, reps_i = reps_i)))
  structural <- mean(comp, na.rm=T)
  structural <- ifelse(structural == 0, 1, structural)
  structural <- ifelse(structural > 1, 1, structural)
  structural <- log(structural) * -1
  return(structural)
  }

#' Structural diversity
#'
#' @description The function calculates the measure of structural diversity as defined by \insertCite{Broekel2019;textual}{GeoInno} from patent data. The primary input should is a data.frame in long-format with five columns. A column *appln_id lists* patents' ids numbers, with the same id appearing as many times as patents are associated to unique CPC (corporate patent classification) classes. Column *year* specifies the year of the patent application. Column *tech* contains the name of the aggregated technology, usually the 2 or 4 digit CPC code. Lastly, column *cpc* features the CPC code (10-digits) of the patent application. The function is set-up for parallel computing using the future.apply() approach. A progress bar from progressr is implemented. An example is provided below showcasing the way the data is to be structured.
#' @param p.dat A data frame in the form of \emph{pat_df}.
#' @param mw Parameter setting the length of moving window. Default set to 3 implying that patent data will be pooled across three years (t, t+1, t+2).
#' @param node.sample The number of nodes sampled in the Network Diversity Score calculation, set to 125, see \insertCite{Emmert-Streib2012;textual}{GeoInno}.
#' @param reps The number of repetitions used in the bootstrap. Default set to 200.
#' @param core.workers The number of cores used in parallelization. Default set to 1.
#'
#' @return The function returns a data.frame with the name complexity including the following information:
#' \itemize{
#'   \item{tech}{The focal technology.}
#'   \item{year}{The year for which the calculation is done.}
#'   \item{tech_nodes}{The number of nodes in the technology-specific binarized co-occurrence (combinatorial) network of the focal technology.}
#'   \item{tech_edges}{The number of edges in the technology-specific binarized co-occurrence (combinatorial) network of the focal technology.}
#'   \item{structural}{The values of the measure of structural diversity.}
#'   \item{patents}{The number of unique patent applications associated with this technology.}
#'   \item{cpcs}{The number of unique CPC codes associated with patents of this technology.}
#' }
#' @export
#'
#'@references
#' \insertAllCited{}
#'
#'
#' @examples
#' structural_diversity(pat_df)
structural_diversity <- function(p.dat = pat_df, mw = 3, node.sample = 125, reps = 200, core.workers = 1)
  {
  pat_df<-NULL
  year <- NULL
  tech <- NULL
  appln_id <- NULL
  geo <- NULL
  structural <- NULL
  item1 <- NULL
  item2 <- NULL
  id <- NULL
  cpc <- NULL
  #for year t, every patent with year in interval t:(t-mw+1) is duplicated into the year t, so that group_by for t includes also all older patents
  p.dat <- p.dat %>% dplyr::rowwise() %>% dplyr::mutate(window=paste(rep(year, mw)+c(0:(mw-1)),collapse="_"))
  p.dat <- p.dat %>% tidyr::separate_rows(window, sep="_")
  p.dat <- p.dat %>% dplyr::mutate(window=as.numeric(window))
  results_year <- p.dat %>% dplyr::group_by(year, tech) %>% dplyr::summarise(patents_year=dplyr::n_distinct(appln_id),
                                                             cpcs_year=dplyr::n_distinct(cpc))  %>%  dplyr::ungroup()

  results_window <- p.dat %>% dplyr::group_by(window, tech) %>% dplyr::summarise(patents_window=dplyr::n_distinct(appln_id),
                                                                    cpcs_window=dplyr::n_distinct(cpc))  %>% dplyr::ungroup()
  handlers(global = TRUE)
  handlers("progress")
  p.dat <- p.dat %>% dplyr::mutate(tech_pats=paste(p.dat$tech, p.dat$window, sep="_")) %>% dplyr::distinct()
  split_data <- split(x=p.dat, f=p.dat$tech_pats)
  plan(multisession, workers = core.workers)
  with_progress({
    p <- progressor(steps = length(split_data))
    dfs <- future_lapply(split_data, future.seed=NULL, FUN=function(x, future.label=T, ...)
    {
      p()
      Sys.sleep(.2)
      pats_window <- x %>% dplyr::pull(appln_id) %>% unique()
      if(length(pats_window)>0)
      {
        pat_tech <- p.dat %>% dplyr::filter(appln_id %in% pats_window) %>% widyr::pairwise_count(item = cpc, feature = appln_id, diag=F)
        pat_net <- pat_tech %>% dplyr::select(item1,item2) %>% as.matrix() %>% igraph::graph_from_edgelist(directed=FALSE) %>% igraph::simplify()
        if(igraph::ecount(pat_net)>0)
          {
          g_inv <- giant.component(pat_net)
          node_count <- igraph::vcount(g_inv)
          edge_count <- igraph::ecount(g_inv)
          if(edge_count>=1)
            {
            tibble::tibble(structural = nds(g = g_inv, s = node.sample, reps_i=reps), edges = edge_count, nodes = node_count)
            }else{
              tibble::tibble(structural=NA, edges = edge_count, nodes = node_count)
            }
          }else{
            tibble::tibble(structural=NA,edges=NA,nodes=NA)
          }
        }else{
          tibble::tibble(structural=NA,edges=NA,nodes=NA)
      }
      #complexity_estimation(x, patdat=p.dat, mw=mw, node.sample=node.sample,reps=reps)
    })
  })
  dfs_df <- data.table::rbindlist(dfs, fill = T)
  dfs_df <- dplyr::bind_cols(id=names(dfs), dfs_df)
  dfs_df <- dfs_df %>% tidyr::separate(id,sep = "_",into=c("tech","window")) %>%
                        dplyr::mutate(window = as.numeric(window))
  results_window <- dplyr::left_join(results_window, dfs_df, by=c("tech","window"), na_matches="never")
  results_window <- dplyr::left_join(results_year, results_window,by=c("tech","year"="window"), na_matches="never")
  return(results_window)
  }

#' Complexity frontier
#'
#' @description The function calculates the complexity frontier for a range of spatial units following the procedure described in \insertCite{Mewes2022;textual}{GeoInno}. That is, for each spatial unit and year, the median of the x-percentil (top) most complex knowledge components is calculated. Knowledge components are represented by distinct occurences of technologies on patents. The latter means that each \emph{appln_id} is associated to the same technology only once.
#' @param c.dat A dataframe with columns: \emph{geo}, \emph{appln_id}, \emph{tech}, \emph{year}, and \emph{structural}. Only distinct technology occurrences are to be included per \emph{appln_id}. The example data set \emph{complexity_dat} illustrates the dataframe structure.
#' @param top Percentile for which frontier is to be estimated. E.g. top=0.05 corresponds to the consideration of the 5 percent most complex knowledge components.
#' @param mw Number of years considered in the mowing window in the type of \emph{t:(t+mw-1)}.
#' @return The function returns a data.frame with the name complexity.frontier including the following information:
#' \itemize{
#'   \item{geo}{The spatial unit for which the frontier is calculated.}
#'   \item{year}{The year for which the frontier is calculated.}
#'   \item{vfrontier}{The value of the complexity frontier.}
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



