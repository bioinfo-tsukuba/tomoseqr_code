generateExpressionPattern <- function(
    xLen,
    yLen,
    zLen,
    func,
    tValue = 1,
    fValue = 0
) {
  mask <- array(0, dim = c(xLen, yLen, zLen))
  indexList <- expand.grid(1:xLen, 1:yLen, 1:zLen)
  for (i in seq_along(indexList[, 1])) {
    x <- indexList[i, 1]
    y <- indexList[i, 2]
    z <- indexList[i, 3]
    if (func(x, y, z) == TRUE) {
      mask[x, y, z] <- tValue
    } else {
      mask[x, y, z] <- fValue
    }
  }
  return(mask)
}

maskFunc <- function(x, y, z) {
  return(
    (x - 25)^2 + 5 * (y - 25)^2 + 10 * (z - 25)^2 <= 25^2
  )
}

func1 <- function(x, y, z) {
  return(
    10 * (x - 25)^2 + 50 * (y - 25)^2 + 100 * (z - 25)^2 <= 25^2
  )
}

func3 <- function(x, y, z) {
  return(
    (x - 12.5)^2 + (y - 25)^2 + (z - 25)^2 <= (25 / 16)^2
  )
}
convertToTomoSeqData <- function(expPattern, geneID) {
  tomoX <- apply(expPattern, 1, sum)
  tomoY <- apply(expPattern, 2, sum)
  tomoZ <- apply(expPattern, 3, sum)
  tomoMatrix <- rbind(tomoX, tomoY, tomoZ)
  colnames(tomoMatrix) <- str_c("section", seq_along(tomoX))
  retTibble <- tibble("geneID" = rep(geneID, 3)) %>%
    bind_cols(as_tibble(tomoMatrix))
  return(retTibble)
}

simSeqToTibble <- function(simSeqMat, mask, axis) {
  maskSum <- apply(mask, axis, sum)
  devided <- simSeqMat *
    t(matrix(
      rep(maskSum / mean(maskSum), length(simSeqMat[, 1])),
      nrow = length(maskSum)
    ))
  devided %>%
    as_tibble() %>%
    mutate("geneID" = str_c("gene", seq_along(simSeqMat[, 1]) + 3)) %>%
    select("geneID", everything()) %>%
    rename_all(~c("geneID", str_c("section", 1:50))) %>%
    
    return()
}

func4 <- function(x, y, z) {
    return(
        (x - 12)^2 + (y - 30)^2 + (z - 26)^2 <= 3 |
          (x - 12)^2 + (y - 20)^2 + (z - 26)^2 <= 3
    )
}

func5 <- function(x, y, z) {
  (x - 12)^2 + (y - 27)^2 + (z - 22)^2 <= 4
}