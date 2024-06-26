library(MASS)
library(distr)
library(mvnmle)
library(mvtnorm)
library(mhsmm)
library(gplots)
library(rasterVis)
library(lattice) 
library(shapes)
library(openxlsx)

DSName = "FixHRF_pu01_au05_NoiseSD0.1"
type = "hsmm"
# 设置工作目录
workDir = paste0("E:\\yjj\\scnu_work\\ishsmm\\data\\HMM_HSMM_mhsmm\\", DSName)
setwd(workDir)

savedir = paste0("E:\\yjj\\scnu_work\\ishsmm\\result\\result1\\", DSName)
if (!dir.exists(savedir)){
  dir.create(savedir)
} 



#setwd("E:\\study\\RTest\\hsmm\\csvData")
# 获取文件名数组
filenames <- Sys.glob(file.path(getwd(),'*.csv'))
#filenames <- Sys.glob(file.path(getwd(),'*.csv'))

#View the subject file names to ensure they have read in correctly.
filenames
# 将数据都读进来，都放在datareadIN这list里面
DataReadIn <- lapply(filenames,function(i){
  table<-read.csv(i, header = FALSE, row.names=NULL)
  table<-scale(table)
  return(table)
})
#对于datareadIN这个list的每一项都做标准差，并且按行
#do.call(函数名, 参数列表) rbind()行合并，lbind 列合并
DataFinal<-do.call(rbind, DataReadIn)
DataFinal<-as.data.frame(DataFinal)
#27489*28 标准化了

# 获取tc的总行数
TimeSerLen <- lapply(filenames,function(i){
  table<-read.csv(i, header = FALSE, row.names=NULL)
  Nrows<-nrow(table) # 可以查询多少行
  return(Nrows)
})
Lengths<-do.call(rbind, TimeSerLen)

#初始化隐藏半马尔可夫模型（HSMM）所需的模型参数

#1.想要拟合的状态个数 
number_states <-  4
# 扫描的总次数，代表人数，roi个数，即特征个数
number_runs <- nrow(Lengths)
number_regions <- ncol(DataFinal)

# 2. 首先将数据分成几个section，每个状态初始化均值（mean）和协方差（covariance）矩阵
output_array = list()
seglength <- ceiling(min(Lengths[,1])/number_states)

for(i in 1:number_states)#number_states 将每个人的（seglength*(i-1)，seglength*i）放到一起
{
  output_array[[i]] <- matrix(NA,nrow=1, ncol=number_regions)
  colnames(output_array[[i]]) <- colnames(DataFinal)
  for(j in 1:number_runs){
    if (i == 1){
      output_array[[i]] <- rbind(output_array[[i]],DataFinal[Lengths[1,1]*(j-1) + c(1:35),])
    }
    if (i == 2){
      output_array[[i]] <- rbind(output_array[[i]],DataFinal[Lengths[1,1]*(j-1) + c(36:58, 127:148),])
    }
    if (i == 3){
      output_array[[i]] <- rbind(output_array[[i]],DataFinal[Lengths[1,1]*(j-1) + c(59:98),])
    }
    if (i == 4){
      output_array[[i]] <- rbind(output_array[[i]],DataFinal[Lengths[1,1]*(j-1) + c(99:126),])
    }
    # output_array[[i]] <- rbind(output_array[[i]],DataFinal[Lengths[1,1]*(j-1) + ((1 + seglength*(i-1)):( seglength*i)),])
    
  }
  output_array[[i]] <- output_array[[i]][-1,]
}

#Now find mean and covariance for each section.
MeanEst <- lapply(output_array,function(i){
  means <- mlest(i)$muhat #取对象的均值
  return(means)
})
CovEst <- lapply(output_array,function(i){
  covs <- mlest(i)$sigmahat#取对象的协方差
  return(covs)
})


# 初始化状态转移概率矩阵 和 初始状态概率
initial <- rep(1/number_states, number_states)
#对角线为0，其他为1/3
TransMatrix <- matrix(NA, nrow=number_states, ncol=number_states)
diag(TransMatrix) <- 0
TransMatrix[is.na(TransMatrix)] <- (1/(number_states-1))


##This saves the initial mean and covariances to a list, which is what will be needed to fit the model.
b <- list(mu = MeanEst, sigma = CovEst)

##This puts the data into a list to fit the model.
Datanew <- list(x = DataFinal, N = Lengths[,1])
class(Datanew) <- "hmm.data"



# d <- list(shape = sample(1:5, size=number_states), scale = sample(10:50, size=number_states), type = "gamma")
d <- list(lambda=c(10,20,30,40),shift=c(8,16,16,10),type='poisson') # 10202010 0.6(8,16,16,10)

model <- hsmmspec(init = initial, trans = TransMatrix, parms.emis = b, dens.emis = dmvnorm.hsmm, sojourn = d)
testrun <- hsmmfit(Datanew, model, mstep=mstep.mvnorm, lock.transition=FALSE, maxit=200, lock.d=FALSE, graphical=TRUE) #50 , 100, 20




cov1<-as.matrix(testrun$model$parms.emission$sigma[[1]])
cov2<-as.matrix(testrun$model$parms.emission$sigma[[2]])
cov3<-as.matrix(testrun$model$parms.emission$sigma[[3]])
cov4<-as.matrix(testrun$model$parms.emission$sigma[[4]])

#Convert the covariance matrices to correlation matrices.
Cor1<-cov2cor(cov1)
Cor2<-cov2cor(cov2)
Cor3<-cov2cor(cov3)
Cor4<-cov2cor(cov4)

#Set the color theme for the state plots we are about to make.
my.theme <- BuRdTheme()
my.at <- seq(-1,1,length.out=length(my.theme$regions$col)-1)
my.ckey <- list(at=my.at, col=my.theme$regions$col)

#After running this code, the states should display as images in R.
plot1 <- levelplot(as.matrix(Cor1[1:ncol(Cor1),ncol(Cor1):1]), par.settings=my.theme, at=my.at, colorkey=my.ckey, scales=list(x=list(rot=90)), main=list(label = "State 1", font = 2), xlab="ROI", ylab="ROI")
plot2 <- levelplot(as.matrix(Cor2[1:ncol(Cor2),ncol(Cor2):1]), par.settings=my.theme, at=my.at, colorkey=my.ckey, scales=list(x=list(rot=90)), main=list(label = "State 2", font = 2), xlab="ROI", ylab="ROI")
plot3 <- levelplot(as.matrix(Cor3[1:ncol(Cor3),ncol(Cor3):1]), par.settings=my.theme, at=my.at, colorkey=my.ckey, scales=list(x=list(rot=90)), main=list(label = "State 3", font = 2), xlab="ROI", ylab="ROI")
plot4 <- levelplot(as.matrix(Cor4[1:ncol(Cor4),ncol(Cor4):1]), par.settings=my.theme, at=my.at, colorkey=my.ckey, scales=list(x=list(rot=90)), main=list(label = "State 4", font = 2), xlab="ROI", ylab="ROI")


plot1
plot2
plot3
plot4
# 保存图片为tiff
plotList <- list(plot1, plot2, plot3, plot4)
for (i in 1:4){
  tiff(paste0(savedir, "\\", type, "_state", as.character(i), ".tif"), width = 800, height = 700, units = "px", res = 100)
  print(plotList[[i]])
  # plot(plotList[[i]])
  dev.off()
}

# 保存状态转移为xlsx
stateTransition = testrun[["yhat"]]
data_to_save <- data.frame(stateTransition)
write.xlsx(data_to_save, file = paste0(savedir, "\\", type, "_stateTransition.xlsx"), rowNames = FALSE, colNames = FALSE)

# 保存4个状态的矩阵为xlsx
for (i in 1:4){
  write.xlsx(data.frame(cov2cor(as.matrix(testrun$model$parms.emission$sigma[[i]]))), file = paste0(savedir, "\\", type, "_state", as.character(i), '.xlsx'), rowNames = FALSE, colNames = FALSE)
}