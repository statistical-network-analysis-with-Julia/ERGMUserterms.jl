# Golden fixture: every BUNDLED USER TERM of ERGMUserterms.jl -- and the
# validation harness's brute-force reference -- against the statnet `ergm`
# terms they are mathematically identical to.
#
# Regenerate from the package root (well under a second; nothing here is
# Monte Carlo):
#
#   Rscript test/fixtures/r/userterms_examples.R > test/fixtures/userterms_examples.toml
#
# WHAT THIS FIXTURE PINS, AND WHY AGAINST R RATHER THAN A HAND LITERAL
#
# ERGMUserterms.jl ships five reference terms (`ExampleTerm`, `TemplateTerm`,
# `WeightedEdges`, `DyadCovTerm`, `InteractionTerm`) and a package template
# (`examples/MyTermPackage`'s `ReciprocatedHomophily`). Third-party authors
# copy them, so a wrong number in any of them propagates. Until the 2026-09
# sprint their tests compared them against hand-typed literals and against
# each other (the harness's `_brute_change_stat`), never against an
# independent implementation. Each of the six is, however, exactly a statnet
# term on a suitably prepared network, so their statistics AND every
# add-direction change statistic can be taken from `ergm` itself:
#
#   ERGMUserterms.jl term          statnet term (ergm 4.x)
#   -----------------------------  -----------------------------------------
#   ExampleTerm()                  nodecov("id") with id = vertex index
#                                  (sum over edges of i + j)
#   TemplateTerm(1.0)              edges
#   WeightedEdges()                edgecov(W), W = 2.5 on the 5 weighted arcs,
#                                  1.0 (WeightedEdges' `default`) elsewhere
#   DyadCovTerm(cov)               edgecov(cov), cov[i,j] = 2i + j (asymmetric)
#   InteractionTerm(:a, :b)        edgecov(M), M[i,j] = a_i*b_j + a_j*b_i
#   ReciprocatedHomophily(:group)  mutual(same = "group", diff = FALSE)
#
# What is recorded:
#
#   (a) summary() of the six terms on a deterministic 8-vertex DIRECTED
#       network (the template's own test network plus four arcs, chosen so
#       that every term's change statistic takes more than one distinct
#       value on it: the weight, covariate and interaction matrices vary
#       across dyads, and the mutual/reciprocity term is 1 on exactly the
#       dyads that would complete a same-group mutual pair).
#   (b) the FULL change-statistic array from `ergmMPLE(..., output="array")`
#       (`$predictor`, dims [tail, head, term]), one vector per term,
#       flattened ROW-MAJOR over the off-diagonal ordered dyads:
#         for i in 1..n, for j in 1..n, j != i  ->  56 values per term.
#       ergmMPLE's array is the add-direction change statistic
#       g(y+ij) - g(y-ij) for EVERY dyad including the ties present, i.e.
#       exactly the state-independent quantity ERGMUserterms' `change_stat`
#       contract demands (and `mutual` is 1 on both arcs of an existing
#       same-group mutual pair, not only on the arc that would complete it).
#   (c) the same summary and change-statistic rows for the five bundled
#       terms on the UNDIRECTED projection of that network (each unordered
#       pair with at least one arc becomes one edge; W and cov symmetrised
#       as the (min, max) entry, which is the canonical key the Julia terms
#       read on an undirected network). The undirected array is filled only
#       for tail < head, so 28 values per term, flattened row-major over
#       i < j. This is what pins the `(min, max)` canonical-key branches in
#       `WeightedEdges`/`DyadCovTerm`.
#
# The Julia testset compares (a)/(c) with `compute`, and (b)/(c) BOTH with
# each term's `change_stat` AND with the harness's `_brute_change_stat`
# (toggle-and-recompute). The second comparison is the point of pinning the
# arrays rather than the summaries alone: `validate_term`'s reference
# implementation is thereby itself validated against ergm's C change
# statistics, on every dyad, instead of being trusted.
#
# WHAT IS DELIBERATELY NOT PINNED: the coefficient NAMES R prints
# (`nodecov.id`, `edgecov.W`, `mutual.group`, ...). User terms keep their own
# labels (`example`, `weightedges.weight`, `recip_homophily.group`, ...);
# they are recorded under `r_names_*` for the record and the Julia testset
# asserts they are NOT claimed. The `ergm.userterms` R package itself is not
# a counterpart of anything shipped here (it is a C template, archived from
# CRAN) and is not used.

suppressMessages({
  .libPaths(c(path.expand("~/R/library"), .libPaths()))
  library(ergm)
})

seed <- 20260912   # nothing here is stochastic; recorded for provenance only
set.seed(seed)

n <- 8

# --- the directed network: the template's 7 arcs plus 4 more ------------------
# (2,5)/(5,2) is a mutual pair ACROSS groups (b, a) -> mutual.group must not
# count it; (3,6) and (8,1) are single arcs whose reverse is absent.
el <- rbind(c(1, 3), c(3, 1), c(2, 4), c(4, 2), c(1, 2), c(5, 7), c(6, 8),
            c(2, 5), c(5, 2), c(3, 6), c(8, 1))
net <- network.initialize(n, directed = TRUE)
add.edges(net, el[, 1], el[, 2])

a <- 1:n
b <- 9 - (1:n)
group <- ifelse((1:n) %% 2 == 1, "a", "b")
net %v% "id" <- 1:n
net %v% "a" <- a
net %v% "b" <- b
net %v% "group" <- group

# Weight matrix: 2.5 on five arcs with DISTINCT unordered pairs ((1,3),
# (2,4), (1,2), (5,7), (2,5) -- so the undirected projection also has five
# weighted edges), 1.0 (= WeightedEdges' default weight) on every other cell,
# so a dyad without a stored weight -- the reverse arcs (3,1), (4,2), (5,2)
# included -- contributes exactly what the Julia term's `default` gives it.
weighted <- c(1, 3, 5, 6, 8)
W <- matrix(1.0, n, n)
for (k in weighted) W[el[k, 1], el[k, 2]] <- 2.5
weighted_tails <- el[weighted, 1]
weighted_heads <- el[weighted, 2]

# Asymmetric dyadic covariate and the interaction matrix
cov <- outer(1:n, 1:n, function(i, j) 2 * i + j)
M <- outer(a, b) + t(outer(a, b))         # M[i,j] = a_i*b_j + a_j*b_i

f_dir <- net ~ nodecov("id") + edges + edgecov(W) + edgecov(cov) + edgecov(M) +
  mutual(same = "group", diff = FALSE)
sum_dir <- summary(f_dir)
arr_dir <- ergmMPLE(f_dir, output = "array")$predictor
stopifnot(dim(arr_dir) == c(n, n, 6))

# --- the undirected projection ------------------------------------------------
# One edge per unordered pair with at least one arc (deduplicated: add.edges
# would otherwise store (1,3) and (3,1) as two parallel edges).
pairs <- unique(cbind(pmin(el[, 1], el[, 2]), pmax(el[, 1], el[, 2])))
und <- network.initialize(n, directed = FALSE)
add.edges(und, pairs[, 1], pairs[, 2])
und %v% "id" <- 1:n
und %v% "a" <- a
und %v% "b" <- b
sym_minmax <- function(X) {
  Y <- X
  for (i in 1:n) for (j in 1:n) Y[i, j] <- X[min(i, j), max(i, j)]
  Y
}
Wu <- sym_minmax(W)
covu <- sym_minmax(cov)
f_und <- und ~ nodecov("id") + edges + edgecov(Wu) + edgecov(covu) + edgecov(M)
sum_und <- summary(f_und)
arr_und <- ergmMPLE(f_und, output = "array")$predictor
stopifnot(dim(arr_und) == c(n, n, 5))

# --- flatteners --------------------------------------------------------------
# `x + 0` turns the `-0` ergm's C code leaves in some change statistics
# into a plain 0 (TOML accepts `-0`, but it reads as a bug).
num <- function(x) paste(sprintf("%.17g", x + 0), collapse = ", ")
ints <- function(x) paste(sprintf("%d", as.integer(x)), collapse = ", ")
strs <- function(x) paste(sprintf('"%s"', x), collapse = ", ")
# Row-major over ordered off-diagonal dyads (directed) / i < j (undirected)
flat_dir <- function(A, k) {
  out <- numeric(0)
  for (i in 1:n) for (j in 1:n) if (i != j) out <- c(out, A[i, j, k])
  stopifnot(length(out) == n * (n - 1), !anyNA(out))
  out
}
flat_und <- function(A, k) {
  out <- numeric(0)
  for (i in 1:n) for (j in 1:n) if (i < j) out <- c(out, A[i, j, k])
  stopifnot(length(out) == n * (n - 1) / 2, !anyNA(out))
  out
}
dir_keys <- c("example", "template", "weightededges", "dyadcov", "interaction",
              "recip_homophily")
und_keys <- dir_keys[1:5]

cat('name = "userterms_examples"\n\n')

cat("[provenance]\n")
cat(sprintf('r_version = "%s"\n', as.character(getRversion())))
cat(sprintf('ergm_version = "%s"\n', as.character(packageVersion("ergm"))))
cat(sprintf('network_version = "%s"\n', as.character(packageVersion("network"))))
cat(sprintf("seed = %d\n", seed))
cat('script = "test/fixtures/r/userterms_examples.R"\n')
cat(sprintf('date = "%s"\n', format(Sys.Date())))
cat('datasets = "8-vertex directed network from a literal edge list (examples/MyTermPackage test network: (1,3),(3,1),(2,4),(4,2),(1,2),(5,7),(6,8) plus (2,5),(5,2),(3,6),(8,1)); vertex attributes id = 1:8, a = 1:8, b = 9 - (1:8), group = a,b alternating; W = 2.5 on the arcs (1,3),(2,4),(1,2),(5,7),(2,5) and 1.0 elsewhere; cov[i,j] = 2i + j; M[i,j] = a_i*b_j + a_j*b_i; and its undirected projection (one edge per unordered pair, W and cov symmetrised to their (min,max) entry)"\n')
cat('model_directed = "net ~ nodecov(\\"id\\") + edges + edgecov(W) + edgecov(cov) + edgecov(M) + mutual(same=\\"group\\", diff=FALSE) -- summary() and ergmMPLE(output=\\"array\\")$predictor; no estimation"\n')
cat('model_undirected = "und ~ nodecov(\\"id\\") + edges + edgecov(Wu) + edgecov(covu) + edgecov(M) -- summary() and ergmMPLE(output=\\"array\\")$predictor; no estimation"\n')
cat('term_map = "ExampleTerm() = nodecov(id); TemplateTerm(1.0) = edges; WeightedEdges() = edgecov(W); DyadCovTerm(cov) = edgecov(cov); InteractionTerm(:a, :b) = edgecov(M); ReciprocatedHomophily(:group) = mutual(same=group, diff=FALSE)"\n')
cat("\n")

cat("[tolerance]\n")
cat("# Every value is a DETERMINISTIC function of a fixed graph and fixed\n")
cat("# integer/rational covariates (integers, halves and their sums): no\n")
cat("# estimator, no Monte Carlo, no optimizer. R's C code and the Julia terms\n")
cat("# must agree to machine precision; 1e-9 is a formality. Any disagreement\n")
cat("# is a bug in a term formula, in a canonical (min,max) key, or in the\n")
cat("# harness's brute-force reference -- do not loosen.\n")
cat("summary_directed = 1e-9\n")
cat("summary_undirected = 1e-9\n")
for (k in dir_keys) cat(sprintf("change_%s_directed = 1e-9\n", k))
for (k in und_keys) cat(sprintf("change_%s_undirected = 1e-9\n", k))
cat("\n")

cat("[values]\n")
cat("# --- the network (tail -> head, 1-based), so the Julia test rebuilds it\n")
cat("# from the TOML alone; the (weighted_tails, weighted_heads) arcs carry\n")
cat("# `weight`, every other arc has no stored weight ----------------------------\n")
cat(sprintf("n = %d\n", n))
cat(sprintf("tails = [%s]\n", ints(el[, 1])))
cat(sprintf("heads = [%s]\n", ints(el[, 2])))
cat(sprintf("weighted_tails = [%s]\n", ints(weighted_tails)))
cat(sprintf("weighted_heads = [%s]\n", ints(weighted_heads)))
cat("weight = 2.5\n")
cat(sprintf("a = [%s]\n", ints(a)))
cat(sprintf("b = [%s]\n", ints(b)))
cat(sprintf("group = [%s]\n", strs(group)))
cat("\n# --- (a) directed summary, in the order ExampleTerm, TemplateTerm(1.0),\n")
cat("# WeightedEdges(), DyadCovTerm(cov), InteractionTerm(:a,:b),\n")
cat("# ReciprocatedHomophily(:group). `r_names_directed` are R's labels: they\n")
cat("# are NOT the Julia terms' names and the testset asserts they are not\n")
cat("# claimed ---------------------------------------------------------------\n")
cat(sprintf("r_names_directed = [%s]\n", strs(names(sum_dir))))
cat(sprintf("summary_directed = [%s]\n", num(sum_dir)))
cat("\n# --- (b) directed add-direction change statistics from ergmMPLE's array,\n")
cat("# one vector per term, ROW-MAJOR over ordered dyads: for i in 1:8, for j\n")
cat("# in 1:8, j != i (56 values) -----------------------------------------------\n")
for (k in seq_along(dir_keys))
  cat(sprintf("change_%s_directed = [%s]\n", dir_keys[k], num(flat_dir(arr_dir, k))))
cat("\n# --- (c) undirected projection: edge list (i < j), summary, and change\n")
cat("# statistics ROW-MAJOR over unordered dyads i < j (28 values) --------------\n")
cat(sprintf("und_tails = [%s]\n", ints(pairs[, 1])))
cat(sprintf("und_heads = [%s]\n", ints(pairs[, 2])))
cat(sprintf("r_names_undirected = [%s]\n", strs(names(sum_und))))
cat(sprintf("summary_undirected = [%s]\n", num(sum_und)))
for (k in seq_along(und_keys))
  cat(sprintf("change_%s_undirected = [%s]\n", und_keys[k], num(flat_und(arr_und, k))))
