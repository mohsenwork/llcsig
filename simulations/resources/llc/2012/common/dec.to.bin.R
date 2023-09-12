dec.to.bin <- function( x, ndigits ) {
  #Changed a decimal number into a binary number with ndigits digits.
  Base.b <- array(NA, dim=c(1, ndigits))
  for(i in 1:ndigits){
    Base.b[, ndigits-i+1] <- (x %% 2)
    x <- (x %/% 2)
  }
  Base.b[1, ]
}