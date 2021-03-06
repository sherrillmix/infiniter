context("infiniter functions")
test_that("Test generateFakeClusters",{
  expect_equal(length(generateFakeClusters()), 5)
  expect_equal(names(generateFakeClusters()), c('ipds','centers','seqs','pwms','ids'))
  tmp<-generateFakeClusters(nIpd=3,nPoints=40,nBase=10,nCluster=4)
  expect_equal(dim(tmp[['ipds']]), c(4*40,3))
  expect_equal(length(tmp[['seqs']]), 4*40)
  expect_equal(nchar(tmp[['seqs']]), rep(10,4*40))
  expect_equal(length(tmp[['pwms']]), 4)
  expect_equal(unlist(unique(lapply(tmp[['pwms']],dim))), c(10,4))
  expect_equal(dim(tmp[['centers']]), c(4,3))
  expect_equal(length(tmp[['ids']]), 4*40)
  expect_equal(unique(table(tmp[['ids']])), 40)
  tmp<-generateFakeClusters(nIpd=4,nPoints=30,nCluster=5,withinSd=0,betweenSd=0)
  expect_equal(unique(as.vector(tmp[['ipds']])), 0)
  expect_equal(unique(as.vector(tmp[['centers']])), 0)
})
