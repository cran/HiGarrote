############################################### Simple functions written in R
## transfer all the columns of D into numeric vectors
#' @noRd
check_D <- function(D) {
  D.colname <- colnames(D)
  D <- apply(D, 2, as.numeric)
  colnames(D) <- D.colname
  D <- data.frame(unique(D))
  return(D)
}

## check if columns are equally spaced
#' @noRd
evenly_spaced <- function(quanti_uni_level) {
  result <- lapply(quanti_uni_level, function(j) {
    length(unique(diff(j))) == 1
  })
  result <- unlist(result)
  return(result)
}


## calculate the distance between two points
#' @noRd
h_dist <- function(x, my_contrast, two_level, qualitative) {
  # x: vector
  # two_level: TRUE/FALSE
  # qualitative: TRUE/FALSE
  
  if(two_level) {
    h <- as.matrix(dist(x))
    h <- ifelse(h != 0, 1, 0)
    return(h)
    
  }else if(qualitative) {
    x <- as.factor(x)
    if(!is.null(my_contrast)) {
      contrasts(x) <- my_contrast
    } else {
      contrasts(x) <- contr.helmert(levels(x))
    }
    m <- model.matrix(~., data.frame(x))
    m <- m[,-1]
    h <- lapply(1:ncol(m), function(i) {
      a <- as.matrix(dist(m[,i]))
      a <- ifelse(a != 0, 1, 0)
      return(a)
    })
    return(h)
    
  }else {
    h <- as.matrix(dist(x))
    return(h)
  }
}



## weak heredity
#' @noRd
gweak <- function(U) {
  effects.name <- colnames(U)
  # effects id
  me.idx <- which(!stringr::str_detect(effects.name, ":"))
  hoe.idx <- which(stringr::str_detect(effects.name, ":"))
  # effects name
  me.names <- effects.name[me.idx]
  hoe.names <- stringr::str_split(colnames(U)[hoe.idx], ":")
  # effects num
  m.eff.num <- length(me.idx)
  h.eff.num <- length(hoe.idx)
  mat = mat.or.vec(m.eff.num, h.eff.num)
  if(h.eff.num != 0) {
    for(i in 1:h.eff.num){
      mat[,i] <- as.numeric(me.names %in% hoe.names[[i]])
    }
  }
  return(cbind(-1,diag(m.eff.num+h.eff.num),
               rbind(mat,-diag(h.eff.num))))
}

## strong heredity
#' @noRd
gstrong <- function(U) {
  effects.name <- colnames(U)
  # effects id
  me.idx <- which(!stringr::str_detect(effects.name, ":"))
  hoe.idx <- which(stringr::str_detect(effects.name, ":"))
  # effects name
  me.names <- effects.name[me.idx]
  hoe.names <- stringr::str_split(colnames(U)[hoe.idx], ":")
  hoe.names.unls <- unlist(hoe.names)
  # effects num
  m.eff.num <- length(me.idx)
  h.eff.num <- length(hoe.idx)
  mat <- mat.or.vec((m.eff.num+h.eff.num), length(hoe.names.unls))
  h.each.num <- lengths(hoe.names)
  h.cum.num <- c(0,cumsum(h.each.num))
  if(h.eff.num != 0) {
    for(i in seq_along(hoe.idx)){
      mat[hoe.idx[i], (h.cum.num[i]+1):(h.cum.num[i+1])] <- -1
    }
    for(i in seq_along(hoe.names.unls)){
      a <- which(me.names %in% hoe.names.unls[i])
      mat[a,i] <- 1
    }
  }
  return(cbind(-1,diag(m.eff.num+h.eff.num), mat))
}


