test_that("estprop and estprop slow give the same answer", {
    result=estprop(c(1,1),c(2,2),diag(2),diag(2))
    result_slow=estprop(c(1,1),c(2,2),diag(2),diag(2))
    expect_equal(result$result, result_slow$result)
})


test_that("adjust_LD requires n0, n1; gives same answer if sample size constant", {
    library(coloc)
    data(coloc_test_data)
    with(coloc_test_data, {
        LD=D1$LD
        dimnames(LD)=list(D1$snp,D1$snp)
        expect_error(colocPropTest:::adjust_LD(D1,LD))
        D1$type="cc"
        D1$s=.5
        D1$n1 = D1$n0 = D1$N * .5
        expect_error(colocPropTest:::adjust_LD(D1,LD))
        D1$n1 = D1$n0 = rep(D1$N * .5, length(D1$snp))
        aLD=colocPropTest:::adjust_LD(D1,LD)
        expect_identical(LD,aLD)
    })
})
    
test_that("run_proptests produces expected output", {
    library(colocPropTest)
    library(coloc)
    data(coloc_test_data)
    results=with(coloc_test_data, {
        run_proptests(D1,D2,LD=D1$LD,topsnps=D1$snp,maxtests=100)
    })
    expect_s3_class(results, "data.table")
    expect_equal(nrow(results), 100)
    expect_equal(min(results$fdr), 1)
})

          
          
