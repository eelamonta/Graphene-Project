install.packages("ggplot2")
install.packages("pheatmap")
install.packages("dplyr")
install.packages("tidyr")

library(pheatmap)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2) #dev.off() #run this if invalid graphics state error
library(broom)
library(magrittr)

################## SET PARAMETERS ################## 

#set file and directory 
dir <- "Ca_Imaging_Data/220427_CM/Light_Stimulation_CMs/targets/"
fileName = "t59_pva_4.csv"

#set parameters
g=600 #set spike filter size (every read with a max less than this number will be eliminated; defined by 1.5 times baseline with fluo-4)
h=353 #set moving average bin number (average will be taken over every h reads)

#open file
open <- paste0("~/Desktop/",dir)
Data0 <- read.csv("/Users/elamonta/Desktop/Ca_Imaging_Data/230804_CM/230804_t59_analysis/targets2/t59_pva_1.csv")
#Create a data table with only means and times (rowname)
Data1 <- Data0[,grepl("Mean",colnames(Data0))]

#make data frame for analysis stats
folder_stats = data.frame()
folder_stats <- setNames(data.frame(matrix(ncol = 24, nrow = 0)), c("Directory","File","spikefilter_amp", "movavg_binsize", 
                                                                    "num_reads", "num_peaks_avg", "num_peaks_sd", 
                                                                    "TTP_avg", "TTP_sd", "TTPmin_avg", "TTPmin_sd", 
                                                                    "TTPmaxamp_avg", "TTPmaxamp_sd","TTPavgamp_avg", "TTPavgamp_sd",
                                                                    "Freq_Hz_avg", "Freq_Hz_sd","Spearman_avg", "Spearman_sd",
                                                                    "Spearman_4sec_avg", "Spearman_4sec_sd", "Spearman_4sec_min", "Spearman_4sec_max"))
folder_stats <- rbind(folder_stats, data.frame(Directory=dir, File=fileName, spikefilter_amp=g, movavg_binsize=h,
                                               num_reads=NA, num_peaks_avg=NA, num_peaks_sd=NA, 
                                               TTP_avg=NA, TTP_sd=NA, TTPmin_avg=NA, TTPmin_sd=NA, 
                                               TTPmaxamp_avg=NA, TTPmaxamp_sd=NA, TTPavgamp_avg=NA, TTPavgamp_sd=NA,
                                               Freq_Hz_avg=NA, Freq_Hz_sd=NA, Spearman_avg=NA, Spearman_sd=NA,
                                               Spearman_4sec_avg=NA, Spearman_4sec_sd=NA, Spearman_4sec_min=NA, Spearman_4sec_max=NA))

################## DATA EDITING ################## 

#baseline normalization by subtracting the minimum mean value from all other values in the column
Data <- sweep(Data1,2,FUN="-",apply(Data1,2,min))

#convert frame number to seconds (one recording is 30 seconds)
Data0$X <- (30/nrow(Data0))*(Data0$X)

#generate line graph before noise correction
melt_Data <- melt(cbind("X"=Data0$X, Data[,c(1:ncol(Data))]), id.vars = c("X"), variable.name = "sample")
ggplot(melt_Data, aes(x=X, y=value, color=sample)) +  
  geom_line()  +
  xlab("Time (sec)") +
  ylab("Fluorescence") +
  theme(legend.position="none")

#noise correction by moving average
for (i in 1:(ncol(Data)-1)) {
  x = Data[,i]
  #plot(Data[,i], type = "l", lwd = 2, col = "blue", xlab = "x", ylab = "signal")
  avg = length(x)/h
  mov_avg = rep(1/avg, avg)
  signal_movavg = stats::filter(x, mov_avg)
  #plot(signal_movavg, type = "l", lwd = 2, col = "blue", xlab = "x", ylab = "signal")
  Data[,i] <- signal_movavg 
}

#remove any NAs (occur at the front and end due to signal_movavg function)
Data$X <- Data0$X
Data <- Data[rowSums(is.na(Data)) == 0, ]
time<-Data$X

#generate line graph after noise correction
melt_Data <- melt(cbind("X"=Data$X, Data[,c(1:ncol(Data))]), id.vars = c("X"), variable.name = "sample")
ggplot(melt_Data, aes(x=X, y=value, color=sample)) +  
  geom_line()  +
  xlab("Time (sec)") +
  ylab("Fluorescence") +
  theme(legend.position="none")

################## FIND_PEAKS FUNCTION ################## 

#https://github.com/stas-g/findPeaks/blob/master/
#increase m to count more less; decrease m to count more
find_peaks <- function (x, m = 50){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}

################## BASELINE CORRECTION ################## 

#baseline correction to eliminate photobleaching effect
corrected <- data.frame(X=Data$X) #data frame for baseline correction of each read in a file
temp_stats <- data.frame() #data frame used to gather calculated time-to-peak data for each read in a file

#save the column for time in seconds "X"
Data<-Data[,colnames(Data)!="X"]

#for each read in a file
for (name in colnames(Data)){
  readno <- substring(name,5,7)
  read=Data[,name]

  #identify local maxima/minima aka peaks/valleys (increase m to detect fewer peaks)
  peak_mat <- data.frame()
  p1 <- find_peaks(read, m = 10) #finding peaks
  p2 <- c(1, find_peaks(-read, m = 10), length(read)) #finding valleys + start and end points

  if ((length(p1)==0) | (length(p2)==0)) {
    print(paste0("Skipping read ",name," in ", fileName, " because no local maximums or minimums were found..."))
    next() 
  }
  
  #plot to test local maxima/minima assignment
  p <- c(p1, p2) #combining the two
  plot(read, type = 'l', main = paste0('m = ', 5))
  points(p, read[p], col = c(rep('red', length(p1)), rep('blue', length(p2))), pch = 19)

  #calculate the difference between each local maxima and minima (times to peaks)
  for (i in p1) {
   for (j in p2) {
     e <- time[i]
     f <- time[j]
     o <- e-f
     k<- i-j
     m <- read[i]
     n <- read[j]
     peak_mat = rbind(peak_mat, data.frame(read=name, timemax=i, timemin=j, ttp=k, valmax=m, valmin=n, timemax_sec=e, timemin_sec=f, ttp_sec=o))
  }}
  peak_mat <- peak_mat[peak_mat$ttp >= 0, ]

  #keep only the smallest time to peak calculations for every local maxima
  peak_mat2 <- data.frame()
  for (i in p1) {
    peak_sub <- subset(peak_mat, peak_mat$timemax==i)
    l <- min(peak_sub$ttp)
    peak_mat2 = rbind(peak_mat2, peak_sub[peak_sub$ttp == l, ] )
  }
  peak_mat2$valamp <- peak_mat2$valmax - peak_mat2$valmin
  
  #make a data frame of each local minimum + start and end points
  linear <- data.frame(timemin = peak_mat2$timemin, valmin = peak_mat2$valmin)
  linear[nrow(linear)+1,]<-c(length(read), read[length(read)]) #add ending point
  linear <- rbind(c(1, read[1]), linear) #add starting point
  
  #remove duplicate start/end points
  if (linear$timemin[1]==linear$timemin[2]) {
    linear <- linear[2:nrow(linear),]
  }
  if (linear$timemin[nrow(linear)]==linear$timemin[nrow(linear)-1]) {
    linear <- linear[1:nrow(linear)-1,]
  }

  #draw a line between each local minima and subtract that line from the original line so that all local minima are now at 0
  lines <- data.frame("index"=1:nrow(Data), time=time, yline=NA)
  n=1
  for (i in 1:(nrow(linear)-1)) {
    m=n+1
    linear_sub <- linear[n:m, ]
    #plot(linear)
    slopem <- (linear_sub$valmin[1]-linear_sub$valmin[2])/(linear_sub$timemin[1]-linear_sub$timemin[2])
    slopem[is.na(slopem)] <- 0
    b <- linear_sub$valmin[1]-slopem*linear_sub$timemin[1]
    #abline(b, slopem)
    lines$yline[linear_sub$timemin[1]:linear_sub$timemin[2]] <- slopem*linear_sub$timemin[1]:linear_sub$timemin[2]+b
    lines$yline[is.na(lines$yline)] <- 0
    n=n+1
  }

  lines <- na.omit(lines)

  #plot to test if baseline correction worked
  plot(read, type = 'l', main = paste0('m = ', 5))
  points(p, read[p], col = c(rep('red', length(p1)), rep('blue', length(p2))), pch = 19)
  lines(lines$index, lines$yline, col='blue')
  #corrected <- data.frame(corrected, baseline=lines$yline, original=x, new=x-lines$yline) #single test
  
  #add baseline corrected data to final data frame
  corrected <- cbind(corrected, new = read-lines$yline)
}

#rename columns as cells
colnames(corrected)<-gsub("Mean","Cell_",colnames(corrected))
corrected[] <- lapply(corrected, function(x) as.numeric(as.character(x)))

#plot corrected fluorescent profile before baseline normalization
Spikes <- melt(cbind("X"=corrected$X, corrected[,c(1:ncol(corrected))]), id.vars = c("X"), variable.name = "sample")
ggplot(Spikes, aes(x=X, y=value, color=sample)) +  geom_line()  +
  xlab("Time (sec)") +
  ylab("Fluorescence") +
  theme(legend.position="none")

corrected[corrected<0] <- 0

## Baseline normalization by subtracting the minimum mean value from all other values in the column
corrected <- sweep(corrected,2,FUN="-",apply(corrected,2,min))

# Remove columns of inactive cells/background (if max value is < g (defined at top))
colMax <- function(data) sapply(data, max, na.rm = TRUE)
corrected1 <- corrected[, colMax(corrected) >= g]

################## FLUORESCENT PROFILE PLOTS ################## 

#plot corrected fluorescent profile after baseline normalization
Spikes <- melt(cbind("X"=corrected$X, corrected1[,c(1:ncol(corrected1))]), id.vars = c("X"), variable.name = "sample")
ggplot(Spikes, aes(x=X, y=value, color=sample)) +  geom_line()  +
  xlab("Time (sec)") +
  ylab("Fluorescence") +
  theme(legend.position="none")

#sorts samples according to variance, applies new order to original dataset
corrected3 <- corrected1[,colnames(corrected1)!="X"]
MaxSpike <- apply(corrected3,2,var)
Amp <- corrected3[,unlist(MaxSpike) >= g ]
preAmp <- data.frame(MaxSpike)
Amp <- corrected3[,order(-preAmp$MaxSpike),]
Amp$X <- time

#select florescence readings with greatest amplitudes to plot
Spikes <- melt(cbind("X"=Amp$X, Amp[,c(1:5)]), id.vars = c("X"), variable.name = "sample")

#plot fluorescence profile for the 5 reads with the greatest amplitudes
print(ggplot(Spikes, aes(x=X, y=value, color=sample)) +  geom_line()  +
        ggtitle(fileName) +
        xlab("Time (sec)") +
        ylab("Fluorescence") +
        theme(legend.position="none"))

#plot fluorescent reads not overlaid
Spikes <- corrected1[,c(1:ncol(corrected1))]
new_spikes <- data.frame(matrix(ncol = 5, nrow = nrow(Spikes)))
n=1
max=0

for (i in colnames(Spikes)){
  temp_spikes <- data.frame(Spikes[,n])
  temp_spikes[, 1] <- temp_spikes[, 1] + max
  max <- max(temp_spikes)
  new_spikes <- cbind(new_spikes, temp_spikes)
  n=n+1
  #print (n)
}
df_clean <- new_spikes[, colSums(is.na(new_spikes)) == 0]
colnames(df_clean)<-colnames(Spikes)
new_spikes <- melt(cbind("X"=time, df_clean), id.vars = c("X"), variable.name = "sample")

ggplot(new_spikes, aes(x=X, y=value, color=sample)) +  geom_line()  +
  xlab("Time (sec)") +
  ylab("Fluorescence") +
  theme_minimal() +
  theme(legend.position="none") 

################## CALCULATE TIME TO PEAKS ################## 

Data <- corrected1
corrected <- data.frame(X=Data$X) #dataframe for baseline correction each read in a file
temp_stats <- data.frame() #dataframe used to gather calculated time-to-peak data for each read in a file

#save the column for time in seconds "X"
Data<-Data[,colnames(Data)!="X"]

#for each read in a file
for (name in colnames(Data)){
  readno <- substring(name,5,7)
  read=Data[,name]
  
  #identify local maxima/minima aka peaks/valleys
  peak_mat <- data.frame()
  p1 <- find_peaks(read, m = 5) #finding peaks
  p2 <- c(1, find_peaks(-read, m = 5), length(read)) #finding valleys + start and end points
  
  if ((length(p1)==0) | (length(p2)==0)) {
    print(paste0("Skipping read ",name," in ", fileName, " because no local maximums or minimums were found..."))
    next() 
  }
  
  #plot local minima/maxima
  #p <- c(p1, p2) #combining the two
  #plot(read, type = 'l', main = paste0('m = ', 5))
  #points(p, read[p], col = c(rep('red', length(p1)), rep('blue', length(p2))), pch = 19)
  
  #calculate the difference between each local maxima and minima (times to peaks)
  for (i in p1) {
    for (j in p2) {
      e <- time[i]
      f <- time[j]
      o <- e-f
      k<- i-j
      m <- read[i]
      n <- read[j]
      peak_mat = rbind(peak_mat, data.frame(read=name, timemax=i, timemin=j, ttp=k, valmax=m, valmin=n, timemax_sec=e, timemin_sec=f, ttp_sec=o))
    }}
  peak_mat <- peak_mat[peak_mat$ttp >= 0, ]
  
  #keep only the smallest time to peak calculations for every local maxima
  peak_mat2 <- data.frame()
  for (i in p1) {
    peak_sub <- subset(peak_mat, peak_mat$timemax==i)
    l <- min(peak_sub$ttp)
    peak_mat2 = rbind(peak_mat2, peak_sub[peak_sub$ttp == l, ] )
  }
  peak_mat2$valamp <- peak_mat2$valmax - peak_mat2$valmin
  
  #add results to dataframes
  temp_stats <- rbind(temp_stats, data.frame(File=fileName, Read=readno, num_peaks=length(p1), TTPavg=mean(peak_mat2$ttp_sec), 
                                             TTPmin=min(peak_mat2$ttp_sec), TTPavg_amp=mean(peak_mat2$valamp), TTPmax_amp=max(peak_mat2$valamp)))
}

#final TTP data 
TTPdata <- temp_stats[!is.na(as.numeric(as.character(temp_stats$Read))),]

################## WHOLE READ STATISTICS ################## 

#calculate frequency in Hz and Spearman Coefficient for each read during the entire 30sec
stats <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("Read", "Time_Start", "Spearman", "Freq_Hz"))
cor_mat <-cor(corrected1, method = "spearman")

for (i in colnames(corrected1)){
  name=i
  readno <- substring(i,5,7)
  read=corrected1[[i]]
  sub_corrected1 <- data.frame (read)
  sub_corrected1$X <- time
  
  #add spearman score for each read
  Avg <- mean(cor_mat[i,][cor_mat[i,]!=1])

  #frequency calculations for each read during the entire 30sec
  list_peaks<-find_peaks(read, m = 5)
  list_num_peaks<-length(list_peaks)
  Freq = list_num_peaks/30
  
  #add to data frame
  stats = rbind(stats, data.frame(Read=name, Time_Start="Whole", Spearman = Avg, Freq_Hz = Freq))
  
}

################## WHOLE FILE STATISTICS AND PHEATMAP ################## 

#create pheatmap for the Spearman Coefficients for each read in the file
corrected2<-corrected1[,colnames(corrected1)!="X"]
#if you only want to include a certain number of reads in the heatmap
#corrected2 <- corrected2[,c(1:15)]
cor_mat2 <-cor(corrected2, method = "spearman")
print(pheatmap(cor_mat2, main= fileName))

#calculates spearman rank correlation coefficient and average frequency for all reads during the entire 30sec read
Avg <- mean(cor_mat[cor_mat!=1])
stddev <- sd(cor_mat[cor_mat!=1])
total_avg_freq <- mean(stats$Freq_Hz)
total_sd_freq <- sd(stats$Freq_Hz)

stats <- rbind(stats, data.frame(Read = "All", Time_Start = "Whole", Spearman = Avg, Freq_Hz = total_avg_freq))

################## INCREMENTAL READ STATISTICS ################## 

#Calculates spearman coefficient for every 4 second interval in the file
#Light stimulation occurs only in an unknown 4 second interval
#Assuming the light stimulation increases synchrony, it should increase the Spearman Coefficient duringa the 4 sec period

n=0
df4 <- data.frame()

#we start at 0 and stop at 26, because 26-30 is our last window to process
for (second in c(0:26)) {
  sub_corrected <- (corrected1[which(time >= second & time < second+4),])
  sub_corrected <- sub_corrected[, colSums(sub_corrected) > 0]
  
  #calculate spearman score
  cor_mat <-cor(sub_corrected[,1:ncol(sub_corrected)], method = "spearman")
  Avg = mean(cor_mat)

  #calculate frequency in Hz
  list_peaks<-sapply(sub_corrected,function(x) find_peaks(x))
  list_num_peaks<-lapply(list_peaks,length)
  list_num_peaks<-list_num_peaks[-length(list_num_peaks)]
  Freq_Avg = mean(unlist(list_num_peaks))/4
  
  #add to data frame
  stats = rbind(stats, data.frame(Read="All", Time_Start=second, Spearman=Avg, Freq_Hz = Freq_Avg))
  df4 = rbind(df4, data.frame(Read="All", Time_Start=second, Spearman=Avg, Freq_Hz = Freq_Avg))
  
}

################## WRITE DATA FILES (*CSV) ################## 
write.csv(TTPdata, file=paste0("~/Desktop/",dir,"/TTP.csv")) 
write.csv(stats, file=paste0("~/Desktop/",dir,"/STATS.csv")) 

