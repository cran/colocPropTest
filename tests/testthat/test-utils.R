test_that("tag requires matrix", {
    expect_error(tag(numeric(5)))
})
test_that("keep_from_S warns if no match", {
    library(coloc)
    data(coloc_test_data)
    with(coloc_test_data, expect_warning(keep_from_S(D1, "nosnphasthissillyname")))
})
test_that("logp gives same answer as standard way", {
    expect_equal(lp(1,1), -log10(pnorm(-1)*2))
})
test_that("plot ellipses checks input data dimensions", {
    expect_error(plot_ellipses(b1=1,vb1=1,b2=1,vb2=1))
    expect_error(plot_ellipses(b1=c(1,1),vb1=1,b2=c(1,1),vb2=1))
    expect_error(plot_ellipses(b1=c(1,1),vb1=diag(3),b2=c(1,1),vb2=diag(3)))
    retval=plot_ellipses(b1=c(1,1),vb1=diag(2),b2=c(1,1),vb2=diag(2))
    expect_identical(class(retval),"list")
})

        
    
