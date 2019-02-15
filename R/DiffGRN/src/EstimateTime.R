PrintEstimateTime <- function(tracked.time, start.time, current.increment,
                              end.increment, save = FALSE) {
  time.snip <- EstimateTime(start.time, current.increment, end.increment)
  cat(paste("Time Remaining:", Seconds2Time(time.snip[1]), "\n"))
  cat(paste("Estimated Completion:", as.POSIXct(time.snip[2], origin = "1970-01-01"), "\n"))
  cat(paste("Percent Complete:", round((100 * time.snip[3]), 2),"%\n"))
  cat(paste("Time Passed:", floor(time.snip[4]), "Seconds", "\n"))
  cat(paste("Average Time:", round(time.snip[5], 1), "Seconds"), "\n")
  cat(paste("Estimated Total Time:", Seconds2Time(time.snip[6]), "\n"))
  if (save) {
    tracked.time[current.increment - 1, 1] = time.snip[1]
    tracked.time[current.increment - 1, 2] = time.snip[2]
    tracked.time[current.increment - 1, 3] = time.snip[3]
    tracked.time[current.increment - 1, 4] = time.snip[4]
    tracked.time[current.increment - 1, 5] = time.snip[5]
    tracked.time[current.increment - 1, 6] = time.snip[6]
    return(tracked.time)
  }
}

EstimateTime <- function(start.time, current.increment, end.increment) {
  start.time <- as.numeric(start.time)
  current.time <- as.numeric(Sys.time())
  percent <- current.increment / end.increment
  passed <- current.time - start.time
  average <- passed / current.increment
  remaining <- (passed / percent) - passed
  complete <- current.time + remaining
  total <- remaining + passed
  values <-  c(remaining, complete, percent, passed, average, total)
  return(values)
}

Seconds2Time <- function(sec) {
  day <- floor(sec / 86400)
  hr <- floor(sec / 3600) - (day * 24)
  min <- floor(sec / 60) - (day * 1440) - (hr * 60)
  sec <- floor(sec) - (day * 86400) - (hr * 3600) - (min * 60)
  return(paste(day, " Days  ", hr, " Hours  ", min, " Minutes  ",
               sec, " Seconds", sep = ""))
}