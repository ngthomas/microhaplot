

# quiets concerns of R CMD check re: the . and other column names
# that appear in dplyr chains
if (getRversion() >= "2.15.1")  {
  utils::globalVariables(
    c(
      ".",
      "V1",
      "V2",
      "allele.balance",
      "depth",
      "group",
      "grp.indx",
      "haplo",
      "id",
      "indx",
      "locus",
      "max.Phred.C",
      "max.uniq.hapl",
      "n.haplo.per.indiv",
      "n.indiv.per.locus",
      "sum.Phred.C",
      "summary.tbl",
      "uniq.id"
    )
  )
}

