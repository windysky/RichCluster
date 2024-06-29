# tests/testthat/test-merge_richsets.R

test_that("merge_richsets works correctly", {
  # Create sample data
  rr1 <- data.frame(Term = "Term1", GeneID = "Gene1", Pvalue = 0.01, Padj = 0.05)
  rr2 <- data.frame(Term = "Term1", GeneID = "Gene2", Pvalue = 0.02, Padj = 0.04)
  rr3 <- data.frame(Term = "Term1", GeneID = "Gene3", Pvalue = 0.03, Padj = 0.03)

  richsets <- list(rr1, rr2, rr3)
  result <- merge_richsets(richsets)

  expect_equal(result$Term, "Term1")
  expect_true(grepl("Gene1,Gene2,Gene3", result$GeneID))
  expect_equal(result$Pvalue, mean(c(0.01, 0.02, 0.03)))
  expect_equal(result$Padj, mean(c(0.05, 0.04, 0.03)))
})
