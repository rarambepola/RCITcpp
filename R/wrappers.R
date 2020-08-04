RIT_wrapper <- function(x, y, n_rff=25, n_bs=100, get_ts=FALSE){
  if(!is.matrix(x)) x <- matrix(x)
  if(!is.matrix(y)) y <- matrix(y)
  #make w
  b <- runif(n_rff)
  w <- matrix(rnorm(dim(x)[2] * n_rff), nrow=n_rff)
  wy <- matrix(rnorm(dim(x)[2] * n_rff), nrow=n_rff)


  x <- normalise_cpp2(x)
  y <- normalise_cpp2(y)

  return(RIT_cpp(x = x,
                 y = y,
                 median_x = median(c(t(dist(x)))),
                 median_y = median(c(t(dist(y)))),
                 w, b, return_ts = get_ts))
}


RCIT_wrapper <- function(x, y, z, n_rff=5, n_rffz=100, n_bs=100, get_ts=FALSE){
  if(!is.matrix(x)) x <- matrix(x)
  if(!is.matrix(y)) y <- matrix(y)
  if(!is.matrix(z)) z <- matrix(z)
  #make w
  b <- runif(n_rff) * 2 * pi
  bz <- runif(n_rffz) * 2 * pi
  by <- runif(n_rff) * 2 * pi
  w <- matrix(rnorm(dim(x)[2] * n_rff), nrow=n_rff)


  x <- normalise_cpp2(x)
  y <- normalise_cpp2(y)
  z <- normalise_cpp2(z)


  y <- cbind(y, z)
  wz <- matrix(rnorm(dim(z)[2] * n_rffz), nrow=n_rffz)
  wy <- matrix(rnorm(dim(y)[2] * n_rff), nrow=n_rff)

  return(RCIT_cpp(x = x,
                  y = y,
                  z = z,
                  median_x = median(c(t(dist(x)))),
                  median_y = median(c(t(dist(y)))),
                  median_z = median(c(t(dist(z)))),
                  w=w,
                  wy=wy,
                  wz=wz,
                  b=b,
                  by=by,
                  bz=bz,
                  return_ts = get_ts,
                  n_bs=n_bs))
}

get_poly_start_end <- function(vals_list){
  vals_list <- lapply(vals_list, function(m) if(!is.matrix(m)) return(matrix(m)))
  x_size <- sapply(vals_list, function(x) dim(x)[1])
  x_size_sum <- cumsum(x_size)

  polygon_start_index <- c(0, x_size_sum)
  polygon_start_index <- polygon_start_index[-length(polygon_start_index)]
  polygon_end_index <- x_size_sum - 1
  return(list(polygon_start_index,
              polygon_end_index))
}

RCIT_disag_wrapper <- function(x_poly, y, z,
                               polygon_start_index,
                               polygon_end_index,
                               n_rff=5, n_rffz=100, n_bs=100, get_ts=FALSE,
                               population=NULL){
  x <- x_poly
  if(!is.matrix(x)) x <- matrix(x_poly)
  if(!is.matrix(y)) y <- matrix(y)
  if(!is.matrix(z)) z <- matrix(z)

  population_supplied <- !is.null(population)
  if(!population_supplied) population <- c(0, 0)

  #make w
  b <- runif(n_rff) * 2 * pi
  bz <- runif(n_rffz) * 2 * pi
  by <- runif(n_rff) * 2 * pi
  w <- matrix(rnorm(dim(x)[2] * n_rff), nrow=n_rff)


  x <- normalise_cpp2(x)
  y <- normalise_cpp2(y)
  z <- normalise_cpp2(z)


  y <- cbind(y, z)
  wz <- matrix(rnorm(dim(z)[2] * n_rffz), nrow=n_rffz)
  wy <- matrix(rnorm(dim(y)[2] * n_rff), nrow=n_rff)

  cpp_version <- RCIT_disag_cpp(x = x,
                                y = y,
                                z = z,
                                median_x = median(c(t(dist(x)))),
                                median_y = median(c(t(dist(y)))),
                                median_z = median(c(t(dist(z)))),
                                w=w,
                                wy=wy,
                                wz=wz,
                                b=b,
                                by=by,
                                bz=bz,
                                polygon_start_index = polygon_start_index,
                                polygon_end_index = polygon_end_index,
                                return_ts = get_ts,
                                n_bs=n_bs,
                                population=population,
                                population_supplied=population_supplied)
  return(cpp_version)
}

RIT_disag_wrapper <- function(x_poly, y,
                              polygon_start_index,
                              polygon_end_index,
                              n_rff=5, n_bs=100, get_ts=FALSE,
                              population=NULL){
  x <- x_poly
  if(!is.matrix(x)) x <- matrix(x_poly)
  if(!is.matrix(y)) y <- matrix(y)

  population_supplied <- !is.null(population)
  if(!population_supplied) population <- c(0, 0)

  #make w
  b <- runif(n_rff) * 2 * pi
  w <- matrix(rnorm(dim(x)[2] * n_rff), nrow=n_rff)


  x <- normalise_cpp2(x)
  y <- normalise_cpp2(y)



  cpp_version <- RIT_disag_cpp(x = x,
                               y = y,
                               median_x = median(c(t(dist(x)))),
                               median_y = median(c(t(dist(y)))),
                               w=w,
                               b=b,
                               polygon_start_index = polygon_start_index,
                               polygon_end_index = polygon_end_index,
                               return_ts = get_ts,
                               n_bs=n_bs,
                               population=population,
                               population_supplied=population_supplied)
  return(cpp_version)
}
