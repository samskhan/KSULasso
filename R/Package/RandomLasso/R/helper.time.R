.helper.time.remaining <- function(pb, start.time, current.increment,
                                   end.increment) {
  setTxtProgressBar(pb, current.increment)
  passed <- as.numeric(Sys.time()) - start.time
  remaining <- (passed / (current.increment / end.increment)) - passed
  hr <- floor(remaining / 3600)
  min <- floor(remaining / 60) - (hr * 60)
  sec <- floor(remaining) - (hr * 3600) - (min * 60)
  if (hr < 10 && hr > 0) {hr = paste("0", hr, sep = "")}
  if (min < 10 && min > 0) {min = paste("0", min, sep = "")}
  if (sec < 10 && sec > 0) {sec = paste("0", sec, sep = "")}

  cat("\r", paste("[",hr, ":", min, ":", sec, "] |", sep = ""))
}
