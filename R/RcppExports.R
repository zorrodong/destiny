# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

censoring_impl <- function(data, sigmas, dists, thr_or_null, uncertain_or_null, missing_or_null, callback) {
    .Call('_destiny_censoring_impl', PACKAGE = 'destiny', data, sigmas, dists, thr_or_null, uncertain_or_null, missing_or_null, callback)
}

predict_censoring_impl <- function(data, data2, thr, uncertain, missing, sigma) {
    .Call('_destiny_predict_censoring_impl', PACKAGE = 'destiny', data, data2, thr, uncertain, missing, sigma)
}

knn_cross <- function(data, query, k, distance = "euclidean") {
    .Call('_destiny_knn_cross', PACKAGE = 'destiny', data, query, k, distance)
}

knn_asym <- function(data, k, distance = "euclidean") {
    .Call('_destiny_knn_asym', PACKAGE = 'destiny', data, k, distance)
}

rank_mat <- function(x) {
    .Call('_destiny_rank_mat', PACKAGE = 'destiny', x)
}

