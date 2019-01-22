

# Testing set:
# haplo.sum<- readRDS("inst/shiny/haPLOType/satrovirens_7_19_16.rds")  %>% mutate(id = as.character(id)) %>% tbl_df()
# if(ncol(haplo.sum)==11) colnames(haplo.sum) <- c("group","id", "locus", "haplo", "depth", "logP.call", "logP.miscall", "pos", "sumP.call","allele.balance","rank")
# if(ncol(haplo.sum)==10) colnames(haplo.sum) <- c("group","id", "locus", "haplo", "depth", "logP.call", "logP.miscall", "pos","allele.balance","rank")
#
# collect.data <- RunSrMicrohap(haplo.sum, "tag_id_914", 2000)
#
# haplo.tbl <- haplo.sum
# locus <- "tag_id_914"
# prior.model<-"uniform"
# n.sam <- 2000
# random.seed <- 43454
# min.read.depth <- 0
#@param n.burn integer. number of burn-in cycle. n.burn < n.sample. Default sets as 0. Optional.

#' Microhaplotype sampler
#'
#' This hierarchical Gibbs sampler infers the true haplotype of grouped individuals of a given locus.
#' @param haplo.tbl data frame. A haplotype data frame generated from running `runHaplo`.  The table contains "group" label, individual "id" label,
#'  "locus" label, observed "haplo"type, haplotype read "depth", "logP.call", "logP.miscall", "pos", "allele.balance","rank", etc.
#' @param locus string. the locus name. Required
#' @param n.sam integer. number of iterations. required
#' @param random.seed integer. Set the random seed number. Default sets as 43454. Optional
#' @param prior.model String. Choose two different prior models: "uniform" - prior set all prior haplotype weight to 1, or "empirical": the prior
#' alpha values are defined by the number of observed cases under refinement
#' @param min.read.depth integer. Set minimal read depth. Default as 0, as no minimal.
#' @export
#' @examples
#' # collect.data <- RunSrMicrohap(haplo.sum, "tag_id_1377", 2)
RunSrMicrohap <- function(haplo.tbl,
           locus,
           n.sam,
           random.seed = 43454,
           prior.model = "uniform",
           min.read.depth = 0)
  {
    set.seed(random.seed)

    # reformat haplot.tbl
    haplo.tbl <- tidyHaplo(haplo.tbl, locus, min.read.depth)

    #initialize parameters
    param <-
      setParam(haplo.tbl, locus, n.sam, random.seed, prior.model)

    if (is.null(param))
      return()

    #precompute all
    match.matrix <- PreComputeMatching(param, haplo.tbl)

    ##call updating for N.sam iterations
    for (i in 1:n.sam)
    {
      param$f <- UpdateF(param)
      param$pf <- UpdatePf(param, haplo.tbl)
      param$h <- UpdateH(param, match.matrix)

      param$save.freq[[i]] <- param$f
      param$save.pfreq[[i]] <- param$pf
      param$save.hap[[i]] <- param$h
    }

    return(param)
  }

tidyHaplo <- function(haplo.tbl, locus.select, min.read.depth) {
  haplo.tbl <- haplo.tbl %>% dplyr::filter(locus == locus.select)
  n.sites <- nchar(haplo.tbl$haplo[1])
  haplo.tbl %>%
    #tidyr::separate(haplo, paste0("haplo.",1:n.sites), sep="(?!^)", extra="drop", remove=F) %>%
    #tidyr::separate(logP.call, paste0("logC.",1:n.sites), sep=",") %>%
    #tidyr::separate(logP.miscall, paste0("logW.",1:n.sites), sep=",") %>%
    dplyr::mutate(uniq.id = as.numeric(factor(id, levels = unique(haplo.tbl$id)))) %>%
    dplyr::filter(depth >= min.read.depth)
}

setParam <-
  function(haplo.tbl,
           locus,
           n.sam,
           random.seed,
           prior.model) {
    param <- NULL

    param$locus.select <- locus
    param$n.sam <- n.sam
    param$random.seed <- random.seed
    param$prior.model <- prior.model
    param$indiv.id <- unique(haplo.tbl$id)

    param$group <- unique(haplo.tbl$group)
    param$n.group <- length(param$group)
    param$n.indiv <- length(unique(haplo.tbl$id))
    param$n.sites <- nchar(haplo.tbl$haplo[1])

    param$grp.assoc.indiv <- haplo.tbl %>%
      dplyr::group_by(uniq.id) %>%
      dplyr::summarise(grp.indx = which(group[1] == param$group)) %>%
      dplyr::arrange(uniq.id) %>% dplyr::select(grp.indx) %>% unlist

    # define the sample space of all possible true haplotypes:
    all.haplotype <- haplo.tbl %>%
      dplyr::filter(rank < 3 , allele.balance > 0.3, depth > 10) %>%
      #dplyr::filter(!grepl("[N]", haplo)) %>%
      dplyr::group_by(haplo) %>%
      dplyr::summarise(count = dplyr::n())

    if (is.null(all.haplotype$haplo) || dim(all.haplotype)[1] == 0) {
      cat("No good haplotype candidate")
      return()
    }

    # this section removes any true haplotype with unnecessary "N" variant sites
    haplo.char <-
      matrix(unlist(strsplit(unique(
        all.haplotype$haplo
      ), "")),
      ncol = param$n.sites,
      byrow = T)

    singlet.site <-
      (apply(haplo.char, 2, function(i)
        length(unique(i))) != 1) * 1
    dispose.indx <- ((haplo.char == "N") %*% singlet.site > 0)

    pass.haplo.list <- all.haplotype$haplo[!dispose.indx]
    pass.ct <- all.haplotype$count[!dispose.indx]
    fail.haplo.list <- all.haplotype$haplo[dispose.indx]
    fail.ct <- all.haplotype$count[dispose.indx]

    if (sum(dispose.indx) > 0) {
      for (i in 1:sum(dispose.indx)) {
        n.match <- sapply(1:length(pass.haplo.list), function(j) {
          sum(unlist(strsplit(fail.haplo.list[i], "")) ==
                unlist(strsplit(pass.haplo.list[j], "")))
        })
        best.match <- n.match == max(n.match)
        best.match <- best.match / length(best.match)
        pass.ct <- pass.ct + (fail.ct[i] * best.match)
      }
    }


    param$haplo <- pass.haplo.list
    param$haplo.ct <- pass.ct
    param$n.haplo <- length(all.haplotype$haplo)

    cat("number of haplotype pair:", param$n.haplo, "\n")


    if (param$n.haplo > 1) {
      param$haplo.pair <-
        rbind(t(combn(1:param$n.haplo, 2)), matrix(rep(1:param$n.haplo, 2), ncol =
                                                     2))
    } else {
      param$haplo.pair <- matrix(rep(1:param$n.haplo, 2), ncol = 2)
    }

    param$n.haplo.pair <- dim(param$haplo.pair)[1]

    ###cache the draw
    param$save.freq <- vector("list", n.sam)
    param$save.pfreq <- vector("list", n.sam)
    param$save.hap <- vector("list", n.sam)


    ## initialize all current parameters
    param$alpha <- param$haplo.ct
    if (prior.model == "uniform")
      param$alpha <- rep(1, param$n.haplo)
    param$f <- gtools::rdirichlet(1, param$alpha)

    param$pf <- t(sapply(1:param$n.group, function(i) {
      haplo.ct <- haplo.tbl %>%
        dplyr::filter(group == param$group[i]) %>%
        dplyr::group_by(haplo) %>%
        dplyr::summarise(count = dplyr::n()) %>%
        dplyr::right_join(.,
                          data.frame("haplo" = param$haplo, stringsAsFactors = F),
                          by = "haplo")

      haplo.ct[is.na(haplo.ct)] <- 0
      gtools::rdirichlet(1, haplo.ct$count + 1)
    }))

    if (param$n.haplo.pair == 1)
      param$pf <- t(param$pf)

    param$h <- array(0, dim = c(param$n.indiv, param$n.haplo))

    haplo.select <- haplo.tbl %>%
      dplyr::group_by(uniq.id) %>%
      dplyr::summarise(
        haplo.1.indx = sample(
          1:param$n.haplo,
          1,
          prob = tabulate(match(haplo, param$haplo),
                          nbins = param$n.haplo) +
            param$pf[which(group[1] == param$group), ]
        ),
        haplo.2.indx = sample(
          1:param$n.haplo,
          1,
          prob = tabulate(match(haplo, param$haplo),
                          nbins = param$n.haplo) +
            param$pf[which(group[1] == param$group), ]
        )
      )

    param$h[cbind(haplo.select$uniq.id,
                  haplo.select$haplo.1.indx)] <- 1
    param$h[cbind(haplo.select$uniq.id,
                  haplo.select$haplo.2.indx)] <-
      param$h[cbind(haplo.select$uniq.id,
                    haplo.select$haplo.2.indx)] +
      1

    param$indic.combn <-
      matrix(0, nrow = param$n.haplo , ncol = param$n.haplo.pair)

    param$indic.combn[cbind(param$haplo.pair[, 1],
                            1:param$n.haplo.pair)] <- 1
    param$indic.combn[cbind(param$haplo.pair[, 2],
                            1:param$n.haplo.pair)] <-
      param$indic.combn[cbind(param$haplo.pair[, 2],
                              1:param$n.haplo.pair)] + 1

    return(param)
  }

PreComputeMatching <- function(param, haplo.tbl) {
  sites.matrix <-
    strsplit(haplo.tbl$haplo, "") %>% unlist %>% matrix(ncol = param$n.sites, byrow =
                                                          T)
  logC.matrix <-
    strsplit(haplo.tbl$logP.call, ",") %>%
    unlist %>%
    as.numeric() %>%
    matrix(ncol = param$n.sites, byrow = T)

  # assume the phred call for reads with the same haplotype variant sites
  # more or less the same
  logC.matrix <- logC.matrix / haplo.tbl$depth

  # logI.matrix <-
  #   strsplit(haplo.tbl$logP.miscall, ",") %>% unlist %>% as.numeric() %>% matrix(ncol =
  #                                                                                  param$n.sites, byrow = T)
  #
  # logI.matrix <- logI.matrix /haplo.tbl$depth


  logI.matrix <- log(1-exp(logC.matrix))


  if("sumP.call" %in% colnames(haplo.tbl)) {
  sC.matrix <-
    strsplit(haplo.tbl$sumP.call, ",") %>%
    unlist %>%
    as.numeric() %>%
    matrix(ncol = param$n.sites, byrow = T)

  aveC.matrix <- sC.matrix / haplo.tbl$depth
  aveI.matrix <- 1-aveC.matrix

  # mock
  cat("create a mock\n")
  # aveC.matrix <- matrix(0.999,
  #                       ncol=param$n.sites,
  #                       nrow=dim(haplo.tbl)[1])
  # aveI.matrix <- 1-aveC.matrix
  #
  # logC.matrix <- log(aveC.matrix)
  # locI.matrix <- log(aveI.matrix)

  }
  # cat(logC.matrix[1,1],"\n")



  n.reads <- dim(haplo.tbl)[1]

  # compare read haplotype to reference haplotype
  # the way i have works as long as the reference haplotype does not contain any "N"
  logP.read.match.ref <- sapply(1:param$n.haplo, function(i) {
    ref.matrix <-
      strsplit(param$haplo[i], "") %>% unlist %>% rep(., n.reads) %>% matrix(ncol =
                                                                               param$n.sites, byrow = T)
    rowSums(
      logC.matrix * (sites.matrix == ref.matrix) +
        logI.matrix * (sites.matrix != ref.matrix)
    )
  })

  depth <- matrix(rep(haplo.tbl$depth, param$n.haplo.pair),
                  ncol=param$n.haplo.pair,
                  byrow = F)

  P.read.match.ref <- exp(logP.read.match.ref)

  # case where we marginalize over two possible states of derived haplotypes
  logP.coefficient <- param$indic.combn
  logP.coefficient[logP.coefficient == 1] <- 1/2
  logP.coefficient[logP.coefficient == 2] <- 1

  logP.match.by.indiv <- log((P.read.match.ref  %*% logP.coefficient)) * depth
  logP.match.by.indiv[is.infinite(logP.match.by.indiv)] <- min(logP.read.match.ref)

  # parsimonus assumption #1: Any observed haplotype has an exact match to the
  # underlying truth haplotype is derived from it.
  #instead of mixture prob from both haplotype

  # indx.to.first.hap <- matrix(0, nrow = param$n.haplo , ncol = param$n.haplo.pair)
  # indx.to.first.hap[cbind(param$haplo.pair[, 1],
  #                         1:param$n.haplo.pair)] <- 1
  # logP.match.first.hap <- logP.read.match.ref %*% indx.to.first.hap
  #
  # indx.to.second.hap <- matrix(0, nrow = param$n.haplo , ncol = param$n.haplo.pair)
  # indx.to.second.hap[cbind(param$haplo.pair[, 2],
  #                         1:param$n.haplo.pair)] <- 1
  # logP.match.second.hap <- logP.read.match.ref %*% indx.to.second.hap
  #
  # min.P.indic <- (logP.match.first.hap >= logP.match.second.hap)
  #
  # logP.sel.best.hap <- ( logP.match.first.hap * min.P.indic) +
  #   (logP.match.second.hap * (1-min.P.indic))
  #
  # logP.match.by.indiv <- logP.sel.best.hap
  # depth.from.first.hap <- min.P.indic * haplo.tbl$depth

  #   matrix(log(2) * (param$haplo.pair[,1]!=param$haplo.pair[,2]),
  #                                                   ncol=param$n.haplo.pair,
  #                                                   nrow=n.reads,
  #                                                   byrow = T)


  #logP.coefficient[logP.coefficient == 1] <- log(2)/2
  #logP.match.by.indiv <- (logP.read.match.ref %*% logP.coefficient)

  # impose that any observed haplotype that completely matches with one of the reference haplotype MUST
  # derived from that reference haplotype, thus ignore the possibility that observed haplotype
  # could also be a product of genotype error occurred in the other mismatch haplotype
  # for the case of heterzygous

  # NOTE: This assumption should be relaxed especially for unbalanced reads

  # exact.match <- sapply(1:param$n.haplo, function(i){
  #   haplo.tbl$haplo==param$haplo[i]
  # })
  #
  # exact.match.pair.haplo <- exact.match %*% param$indic.combn > 0
  # exact.match.pair.haplo.2 <- exact.match %*% param$indic.combn == 2
  # exact.match.pair.haplo.1 <- exact.match %*% param$indic.combn == 1
  #
  # logP.match.by.indiv <- (exact.match.pair.haplo.2 *
  #                           rowSums(logP.read.match.ref *exact.match) *
  #                           2 ) + (
  #                             exact.match.pair.haplo.1 *
  #                                rowSums(logP.read.match.ref *exact.match) *
  #                                2*log(2) )+ (
  #                             (!exact.match.pair.haplo) *
  #                               (logP.read.match.ref %*% logP.coefficient))


  # Summing all log P by individual
  index.read.to.indiv <-
    haplo.tbl %>% dplyr::mutate(indx = dplyr::row_number()) %>% dplyr::select(uniq.id, indx) %>% as.matrix(ncol =
                                                                                                      2, byrow = T)
  indic.matrix.read.by.indiv <-
    matrix(0, nrow = param$n.indiv, ncol = n.reads)
  indic.matrix.read.by.indiv[index.read.to.indiv] <- 1

  phred.prob <- indic.matrix.read.by.indiv %*% logP.match.by.indiv


  # n <- matrix(indic.matrix.read.by.indiv %*% haplo.tbl$depth,
  #             ncol=param$n.haplo.pair,
  #             nrow=param$n.indiv,
  #             byrow = F)
  #
  # k <- indic.matrix.read.by.indiv %*% depth.from.first.hap
  #
  # param$alpha <- 50
  # param$beta <- 50
  #
  # beta.balance <- lgamma(n+1) + lgamma(k+param$alpha) + lgamma(n-k+param$beta) +
  #   lgamma(param$alpha+param$beta) - lgamma(param$alpha) - lgamma(param$beta) -
  #   lgamma(k+1) - lgamma(n-k+1) - lgamma(n+param$alpha+param$beta)
  #
  # beta.balance[,colSums(param$indic.combn==2)==1] <- 0

  phred.prob #+ beta.balance

}



UpdateF <- function(param) {
  gtools::rdirichlet(1, param$alpha)
}

UpdatePf <- function(param, haplo.tbl) {
  if (param$n.haplo.pair > 1) {
    t(sapply(1:param$n.group, function(i) {
      gtools::rdirichlet(1, colSums(param$h[param$grp.assoc.indiv == i, ]) + param$f)
    }))
  } else {
    t(sapply(1:param$n.group, function(i) {
      gtools::rdirichlet(1, sum(param$h[param$grp.assoc.indiv == i, ]) + param$f)
    }))
  }
}

UpdateH <- function(param, match.matrix) {
  pf <- param$pf
  if (param$n.haplo.pair == 1)
    pf <- t(pf)

  logP.contrib.pf <- pf %*% param$indic.combn %>%
    .[param$grp.assoc.indiv, ] %>%
    log

  logP.h <- match.matrix + logP.contrib.pf

  # need to normalize to prevent numeric overflow
  prob.distrib <- t(apply(logP.h, 1, function(i) {
    exp(i - max(i))
  }))

  if (param$n.haplo.pair == 1)
    prob.distrib <- t(prob.distrib)

  pair.indx <- sapply(1:param$n.indiv, function(i) {
    sample(param$n.haplo.pair, 1, prob = prob.distrib[i, ])
  })

  h <- array(0, dim = c(param$n.indiv, param$n.haplo))
  h[cbind(1:param$n.indiv,
          param$haplo.pair[pair.indx, 1])] <- 1

  h[cbind(1:param$n.indiv,
          param$haplo.pair[pair.indx, 2])] <-
    h[cbind(1:param$n.indiv,
            param$haplo.pair[pair.indx, 2])] + 1

  h

}


# prior alpha parameter: set as all 1 (weak prior), or set # based on what's observed
