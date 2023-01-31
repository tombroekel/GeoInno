#' Artificial data set to illustrate GeoInno
#'
#' A data set illustrating the form and style how patent data should be structured to be used by functions in this package
#'
#' @format ## `pat.df`
#' A data frame with 1,000 rows and 4 columns:
#' \describe{
#'   \item{appln_id}{application number of patent}
#'   \item{cpc}{IPC/CPC code assigned to a patent}
#'   \item{tech}{Technology code, usually 2 or 4-digit CPC code}
#'   \item{year}{Year of observation}
#' }
#' @source None
"pat.df"

#' Giant component
#'
#' @description Helper function for the calculation of the network diversity score of \insertCite{Emmert-Streib2012;textual}{GeoInno} and the structural diversity complexity measure of \insertCite{Broekel2019;textual}{GeoInno}.
#'
#' @param graph An igraph object from which the giant component needs to be extracted.
#' @param ...
#'
#' @return The function returns an igraph object representing the giant component of a network
#' @export
#'
#' @examples
#' my.graph <- igraph::random.graph.game(p.or.m = 1/10, n=10)
#' giant.component(my.graph)
giant.component <- function(graph, ...)
{
  cl <- clusters(graph, ...)
  induced_subgraph(graph, which(cl$membership == which.max(cl$csize)))
}


#' individual network diversity score (iNDS)
#'
#'@description NDS.intern() calculates the individual network diversity score as defined by \insertCite{Emmert-Streib2012;textual}{GeoInno}. It is in the calculation of the structural diversity complexity measure of \insertCite{Broekel2019;textual}{GeoInno}.
#'
#' @param s The (random) sample of nodes (position indizes in igraph object) for which the partial networks are to be extracted by means of a random walk.
#' @param g The igraph object. Usually, the binarized version of the combinatorial network of CPC classes co-occurring on patents.
#' @param node.sample The number of nodes sampled in the Network Diversity Score calculation, set to 125, see \insertCite{Emmert-Streib2012;textual}{GeoInno}
#' @param reps The number of repetitions used in the bootstrap, default set to 200.
#'
#' @return Returns the value of the iNDS measure.
#' @export
#'
#' @examples
#' my.graph <- igraph::random.graph.game(p.or.m = 1/10, n=10)
#' NDS.intern(g = my.graph, s=c(1:10), node.sample=10, reps=10)
NDS.intern<-function(g, s, node.sample, reps)
{
  set.seed(2)
  sample_vertex <- random_walk(g, start=s, steps=reps, mode="all")
  sample_net <- induced.subgraph(graph = g, vids = sample_vertex)
  set.seed(2)
  modules_g <- cluster_walktrap(sample_net, steps = 4, merges = TRUE, modularity = TRUE, membership = TRUE)
  m <- data.frame(sizes(modules_g))
  lap <- eigen(laplacian_matrix(sample_net,norm=F))$values
  graphletis <- count_graphlets_for_graph(sample_net, 4)
  graph.3 <- sum(graphletis[paste("G",1:2,sep="")])
  graph.4 <- sum(graphletis[paste("G",3:8,sep="")])
  a.module <- length(modules_g)/vcount(sample_net)
  v.module <- var(m$Freq)/mean(m$Freq)
  v.l <- var(lap)/mean(lap)
  r.motif <- graph.3/graph.4
  nds.s <- (a.module*r.motif)/(v.module*v.l)
  ifelse(is.infinite(nds.s)==T | is.na(nds.s) ==T,0,nds.s)
}

nds<-function(g, node.sample, reps)
{
  vsample<-ifelse(vcount(g)<node.sample, vcount(g), node.sample)
  set.seed(2)
  start_vertex<-sample(vcount(g), vsample, replace=F)
  comp<-sort(unlist(lapply(start_vertex,NDS.intern,g=g,node.sample=node.sample, reps=reps)))
  structural<-mean(comp,na.rm=T)
  structural <- ifelse(structural == 0, 1, structural)
  structural <- ifelse(structural > 1, 1, structural)
  structural <- log(structural) * -1
  return(structural)
}

#selector<-function(x)
#{
#  which(substr(V(pat_net)$name,1,4)==x)
#}

complexity_estimation <- function(patdat, pats_window, ...)
  {
  pats_window<-pats_window %>% unique()
  if(length(pats_window)>0)
    {
    pat_tech <- patdat %>% filter(appln_id %in% pats_window) %>%
      widyr::pairwise_count(item=cpc,feature = appln_id, diag=F)

    pat_net <- pat_tech %>% select(item1,item2) %>% as.matrix() %>% igraph:::graph_from_edgelist(directed=FALSE) %>% igraph:::simplify()
    if(igraph:::ecount(pat_net)>0)
    {
      g_inv <- giant.component(pat_net)
      node_count <- igraph:::vcount(g_inv)
      edge_count <- igraph:::ecount(g_inv)
      if(edge_count>=1)
        {
        return(tibble(structural=nds(g_inv, node.sample=node.sample, reps=reps),edges=edge_count,nodes=node_count))
        }else{
              return(tibble(structural=NA, edges=edge_count, nodes=node_count))
              }
    }else{
      return(tibble(structural=NA,edges=NA,nodes=NA))
    }
  }else{
    return(tibble(structural=NA,edges=NA,nodes=NA))
  }
}

#' Structural diversity
#'
#' @description The main function in this script structural_diversity() calculates the measure of structural diversity as defined by \insertCite{Broekel2019;textual}{GeoInno} from patent data. The primary input should is a data.frame in long-format with five columns. A column *appln_id lists* patents' ids numbers, with the same id appearing as many times as patents are associated to unique CPC (corporate patent classification) classes. Column *year* specifies the year of the patent application. Column *tech* contains the name of the aggregated technology, usually the 2 or 4 digit CPC code. Lastly, column *cpc* features the CPC code (10-digits) of the patent application. An example is provided below showcasing the way the data is to be structured.
#' @param patdat A data frame in the form of pat.df .
#' @param mw Parameter setting the length of moving window , default set to 3 implying that patent data will be pooled across three years (t, t+1, t+2).
#' @param node.sample The number of nodes sampled in the Network Diversity Score calculation, set to 125, see \insertCite{Emmert-Streib2012;textual}{GeoInno}.
#' @param reps The number of repetitions used in the bootstrap, default set to 200.
#' @importFrom Rdpack reprompt
#' @import tidyverse
#' @import widyr
#' @import netdist
#'
#' @return The function returns a data.frame with the name complexity including the follwoing information
#' \itemize{
#'   \item{\emph{tech}}{focal technology}
#'   \item{\emph{year}}{year for which the calculation is done}
#'   \item{\emph{tech_nodes}}{number of nodes in the technology-specific binarized co-occurrence (combinatorial) network of the focal technology}
#'   \item{\emph{tech_edges}}{number of edges in the technology-specific binarized co-occurrence (combinatorial) network of the focal technology}
#'   \item{\emph{structural}}{the values of the measure of structural diversity}
#'   \item{\emph{patents}}{number of unique patent applications associated with this technology}
#'   \item{\emph{cpcs}}{number of unique CPC codes associated with patents of this technology.}
#' }
#'@references
#' \insertAllCited{}
#'
#' @export
#'
#' @examples
#' structural_diversity(pat.df)
structural_diversity <- function(patdat, mw=3, node.sample=125, reps=200)
  {
  #for year t, every patent with year in interval t:(t-mw+1) is duplicated into the year t, so that group_by for t includes also all older patents
  patdat <- patdat %>% rowwise() %>% mutate(window=paste(rep(year,mw)+c(0:(mw-1)),collapse="_"))
  patdat <- patdat %>% separate_rows(window, sep="_")
  patdat <- patdat %>% mutate(window=as.numeric(window))
  results_year <- patdat %>% group_by(year, tech) %>% summarise(patents_year=n_distinct(appln_id),
                                                             cpcs_year=n_distinct(cpc))  %>%    ungroup()

  results_window <- patdat %>% group_by(window, tech) %>% summarise(patents_window=n_distinct(appln_id),
                                                           cpcs_window=n_distinct(cpc),
                                                           GeoInno:::complexity_estimation(patdat,appln_id))  %>%
                                                            ungroup()
  results_window <- left_join(results_year,results_window,by=c("tech","year"="window"))
  return(results_window)
  }


#' Create sample data
#'
#' @param num.pat Number of patents
#'
#' @return The function returns a data.frame suited to be used with other functions in the package
#' @export
#'
#' @examples
#' create_sample_data(200)
create_sample_data <- function(num.pat=200)
{
  set.seed(123)
  appln_id <- sample(num.pat, size=num.pat*5, replace = T)
  cpc.1 <- sample(LETTERS[1:10], size=num.pat*5, replace = T)
  cpc.2 <- sample(LETTERS[1:25], size=num.pat*5, replace = T)
  cpc <- paste0(cpc.1,cpc.2)

  pat.df <- data.frame(appln_id=appln_id, cpc=cpc) %>% arrange(appln_id)
  pat.df <- pat.df %>% mutate(tech = substr(cpc,1,1))
  year <- sample(c(1999:2005), size=length(unique(appln_id)), replace = T)
  pat.df <- pat.df %>% mutate(year = year[appln_id])
  return(pat.df)
}
