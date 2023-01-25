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

giant.component <- function(graph, ...)
{
  cl <- clusters(graph, ...)
  induced_subgraph(graph, which(cl$membership == which.max(cl$csize)))
}

NDS.intern<-function(s, g, node.sample, reps)
{
  #select sample network of size 200 through random walk
  set.seed(2)
  sample_vertex <- random_walk(g,start=s,steps=reps,mode="all")
  sample_net <- induced.subgraph(graph = g, vids = sample_vertex)
  #Network diversity score
  set.seed(2)
  modules_g <- cluster_walktrap(sample_net, steps = 4, merges = TRUE, modularity = TRUE, membership = TRUE)
  m <- data.frame(sizes(modules_g))
  lap <- eigen(laplacian_matrix(sample_net,norm=F))$values
  graphletis <- count_graphlets_for_graph(sample_net,4)
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
  start_vertex<-sample(vcount(g),vsample,replace=F)
  comp<-sort(unlist(lapply(start_vertex,NDS.intern,g=g,node.sample=node.sample, reps=reps)))
  structural<-mean(comp,na.rm=T)
  structural <- ifelse(structural == 0, 1, structural)
  structural <- ifelse(structural > 1, 1, structural)
  structural <- log(structural) * -1
  return(structural)
}

selector<-function(x)
{
  which(substr(V(pat_net)$name,1,4)==x)
}

#' Structural diversity
#'
#' @description The main function in this script structural_diversity() calculates the measure of structural diversity as defined by \insertCite{Broekel2019;textual}{GeoInno} from patent data. The primary input should is a data.frame in long-format with five columns. A column *appln_id lists* patents' ids numbers, with the same id appearing as many times as patents are associated to unique CPC (corporate patent classification) classes. Column *year* specifies the year of the patent application. Column *tech* contains the name of the aggregated technology, usually the 2 or 4 digit CPC code. Lastly, column *cpc* features the CPC code (10-digits) of the patent application. An example is provided below showcasing the way the data is to be structured.
#' @param patdat a data frame in the form of pat.df
#' @param mw parameter setting the length of moving window , default set to 3 implying that patent data will be pooled across three years (t, t+1, t+2)
#' @param node.sample the number of nodes sampled in the Network Diversity Score calculation , set to 125, see \insertCite{Emmert-Streib2012;textual}{GeoInno}
#' @param reps the number of repetitions used in the bootstrap, default set to 200
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
  years <- patdat$year %>% unique() %>% sort()
  techs <- patdat$tech %>% unique() %>% sort()

  results<-data.frame("tech"="","year"=0, "tech_nodes"=0, "tech_edges"=0, "structural"=0, stringsAsFactors = F)
  results<-results %>% slice(-1)

  for(t in 1:length(years))
  {
    results.t<-data.frame("tech"=techs,"year"=years[t],"tech_nodes"=rep(0,length(techs)),"tech_edges"=rep(0,length(techs)),"structural"=rep(0,length(techs)),stringsAsFactors = F)

    #Extract moving window data
    upper <- years[t]
    lower <- ifelse((years[t]-mw)< min(years),min(years),years[t-(mw-1)])
    pat.t <- patdat %>% filter(year <= upper & year >= lower)

    #Number of patents and cpcs per tech
    pats <- pat.t %>% group_by(tech) %>% summarise(patents=n_distinct(appln_id),
                                                   cpcs=n_distinct(cpc)) %>%
      ungroup()

    results.t <- left_join(results.t, pats, by=c("tech"="tech"))

    for (i in 1:nrow(results.t))
    {
      tech_appln<-pat.t %>% filter(tech==techs[i]) %>%  pull(appln_id) %>% unique()
      if(length(tech_appln)>0)
      {
        #all patents with at least one CPC of tech and network without isolates
        pat_tech <- pat.t %>% filter(appln_id %in% tech_appln) %>%
          widyr::pairwise_count(item=cpc,feature = appln_id, diag=F)

        pat_net <- pat_tech %>% select(item1,item2) %>% as.matrix() %>% graph_from_edgelist(directed=FALSE) %>% igraph::simplify()
        if(ecount(pat_net)>0)
        {
          g_inv <- giant.component(pat_net)
          results.t$tech_nodes[i] <- vcount(g_inv)
          results.t$tech_edges[i] <- ecount(g_inv)
          if(ecount(g_inv)>=1)
          {
            results.t[i,"structural"]<-nds(g_inv, node.sample=node.sample, reps=reps)
          }
        }
      }
      print(paste(years[t]," ","Structural diversity for ","i=",i,sep=""))
    }
    results<-bind_rows(results,results.t)
  }
  return(results)
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
