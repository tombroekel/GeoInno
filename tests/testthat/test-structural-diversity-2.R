test_that("legacy functions reproduce the master snapshot", {
  baseline <- readRDS(test_path("fixtures", "master_snapshot.rds"))
  expect_identical(baseline$commit, "e4819aeded4bb2112e987dd443a5f01030252043")
  graphs <- list(
    ring8 = igraph::make_ring(8),
    star9 = igraph::make_star(9, mode = "undirected"),
    complete6 = igraph::make_full_graph(6),
    lattice = igraph::make_lattice(c(4, 4))
  )
  helpers <- suppressWarnings(lapply(graphs, function(g) list(
    intern = NDS.intern(node_x = 1, g = g, reps_i = 10),
    nds = nds(g = g, s = 5, reps_i = 10)
  )))
  expect_equal(helpers, baseline$helpers, tolerance = 1e-9)
  expect_equal(
    suppressWarnings(structural_diversity_legacy(pat_df, mw = 3, node.sample = 5,
                                                 reps = 10, core.workers = 1)),
    baseline$example, tolerance = 1e-9)
})

test_that("fixed-node sampling visits exactly the requested unique vertices", {
  g <- igraph::make_lattice(c(10, 10))
  igraph::V(g)$name <- as.character(seq_len(igraph::vcount(g)))
  set.seed(42)
  one <- GeoInno:::sample_fixed_nodes(g, 1, 8)
  two <- GeoInno:::sample_fixed_nodes(g, 1, 8)
  expect_equal(igraph::vcount(one), 8L)
  expect_equal(length(unique(igraph::V(one)$name)), 8L)
  expect_false(identical(sort(igraph::V(one)$name), sort(igraph::V(two)$name)))

  small <- igraph::make_ring(6)
  expect_identical(GeoInno:::sample_fixed_nodes(small, 1, 8), small)
})

test_that("undefined iNDS is missing and finite samples are averaged", {
  complete <- igraph::make_full_graph(6)
  expect_true(is.na(GeoInno:::NDS.intern_fixed(1, complete, n_nodes = 8)))
  all_missing <- GeoInno:::nds_fixed(complete, n_starts = 4, n_nodes = 8)
  expect_true(is.na(all_missing$Structural_raw))
  expect_true(is.na(all_missing$Structural_neglog))
  expect_equal(all_missing$n_samples, 4)
  expect_equal(all_missing$n_valid_samples, 0)
  expect_equal(all_missing$share_valid_samples, 0)

  mixed <- GeoInno:::aggregate_sample_nds(c(2, NA_real_, 4, Inf))
  expect_equal(mixed$Structural_raw, 3)
  expect_equal(mixed$Structural_neglog, -log(3))
  expect_equal(mixed$n_valid_samples, 2)
  expect_equal(mixed$share_valid_samples, .5)
  expect_true(is.na(GeoInno:::aggregate_sample_nds(c(NA_real_, NaN))$Structural_raw))
})

test_that("raw values above one remain and negative logs are valid", {
  set.seed(4)
  result <- GeoInno:::nds_fixed(igraph::make_ring(20), n_starts = 5, n_nodes = 7)
  expect_gt(result$Structural_raw, 1)
  expect_lt(result$Structural_neglog, 0)
  expect_equal(result$Structural_neglog, -log(result$Structural_raw))
  expect_equal(length(unique(result$walk_seeds)), result$n_samples)
  set.seed(4)
  expect_identical(GeoInno:::nds_fixed(igraph::make_ring(20),
                                      n_starts = 5, n_nodes = 7), result)
})

test_that("public estimator is repeatable and uses the published data structure", {
  small_data <- subset(pat_df, year == min(year) & tech %in% c("A", "B", "C"))
  first <- structural_diversity(small_data, mw = 1, n_starts = 3,
                                n_nodes = 8, seed = 123, core.workers = 1)
  again <- structural_diversity(small_data, mw = 1, n_starts = 3,
                                n_nodes = 8, seed = 123, core.workers = 1)
  expect_identical(first, again)
  expect_true(all(c("Structural_raw", "Structural_neglog", "n_samples",
                    "n_valid_samples", "share_valid_samples") %in% names(first)))
  valid <- is.finite(first$Structural_raw) & first$Structural_raw > 0
  expect_equal(first$Structural_neglog[valid], -log(first$Structural_raw[valid]))
  expect_equal(nrow(structural_diversity(mw = 1, n_starts = 2, n_nodes = 5,
                                        seed = 123, core.workers = 1)), 70L)
})

test_that("one and two workers produce identical estimates", {
  socket_ok <- tryCatch({
    cl <- parallel::makePSOCKcluster(2L)
    parallel::stopCluster(cl)
    TRUE
  }, error = function(e) FALSE)
  if (!socket_ok) skip("Local sandbox does not permit parallel worker sockets")
  small_data <- subset(pat_df, year == min(year) & tech %in% c("A", "B", "C"))
  sequential <- structural_diversity(small_data, mw = 1, n_starts = 3,
                                     n_nodes = 8, seed = 123, core.workers = 1)
  parallel <- structural_diversity(small_data, mw = 1, n_starts = 3,
                                   n_nodes = 8, seed = 123, core.workers = 2)
  expect_identical(sequential, parallel)
})
