#加载必要的包
library(ggplot2)
library(reshape2)
library(lattice)
library(hexbin)
library(gridExtra)
library(xtable)
library(splines)
library(survival)
library(grid)
library(lpSolve)
library(pander)
library(Amelia)
library(mice)
library(arm)
#I 数据预处理
d2 <- imputed_dataset
#删除年份这一栏
d2$year <- NULL
d2$daytime<- NULL
#修改force那一栏，转变为01变量
d2$force2 <- d2$force
d2$force2 <- ifelse(d2$force != 0, 1, 0)
# 将 force2 列的值赋给 force 列
d2$force <- d2$force2
# 删除 force2 列
d2$force2 <- NULL
#对cs变量的处理
d2[, (ncol(d2) - 9):ncol(d2)] <- ifelse(d2[, (ncol(d2) - 9):ncol(d2)] == "Y", 1, 0)
# 2. 新建一列数值型变量 cs，数值为对最后十列的求和
d2$cs <- rowSums(d2[, (ncol(d2) - 9):ncol(d2)])
#只保留第1到11，+22列，By setting drop = FALSE, you ensure that the result is still a data frame.
d2 <- d2[, c(1:(ncol(d2) - 11), ncol(d2)), drop = FALSE]

#残差图不一样的问题
d3<-d2
d3$typeofid2<- ifelse(d3$typeofid2 == "R", 0, 1)
# 将 typeofid2 转换为因子型变量
d3$typeofid2 <- as.factor(d3$typeofid2)
d3$force <- as.factor(d3$force)
summary(d3)
#II 开始建立模型
#一、# 建立空模型,this procedure is not to be misused for variable selection, but to supplement the explorative stuides 
form <- force~race2+gender+age2+inout2+ac_incid+ac_time+offunif2+typeofid2+othpers2+cs
form2<- force~race2+gender+ns(age2,df=2)+inout2+ac_incid+ac_time+offunif2+typeofid2+othpers2+cs
form3 <- force~race2+gender+age2+inout2+ac_incid*ac_time+offunif2+typeofid2+othpers2+cs
model1<-glm(form,family = binomial(link = "logit"),data=d3)
model2<-glm(form2,family = binomial(link = "logit"),data=d3)
model3<-glm(form3,family = binomial(link = "logit"),data=d3)

null_model <- glm(force ~ 1 , family = binomial(link = "logit"), data = d2)
#对logistic回归采用似然比检验
oneTermModels<-add1(null_model,form,test="LRT")
#二、终于开始建模
model1<-glm(form,family = binomial(link = "logit"),data=d3) #model1是上面那12个变量的累加模型
#modell<-glm(formm,family = binomial(link = "logit"),data=d3) 
#残差诊断
# 用随机方法抽
# nrow(model1$model) gets the number of rows in the model's dataset
set.seed(123)
subset_size <- 1000000

# Create random indices for the subset
subset_indices <- sample(1:nrow(model3$model), size = min(subset_size, nrow(model3$model)))

# Create a subset model based on the random subset
subset_model3 <- update(model3, subset = subset_indices)
subset_model2 <- update(model2, subset = subset_indices)
subset_model1 <- update(model1, subset = subset_indices)
#确定某一年的残差图,不随机抽
#subset_size <- 160851
#subset_indices <- seq_len(min(subset_size, nrow(model1$model)))
#subset_model3 <- update(model3, subset = subset_indices)

# Create a binned residual plot
binnedplot(fitted(subset_model3, type = "response"), 
           residuals(subset_model3,type="response"), 
           nclass = NULL, 
           xlab = "Expected Values", 
           ylab = "Average residual", 
           main = "Binned residual plot (interaction", 
           cex.pts = 0.8, 
           col.pts = 1,
           col.int = "gray")
binnedplot(fitted(subset_model1), 
           residuals(subset_model1,type="response"), 
           nclass = NULL, 
           xlab = "Expected Values", 
           ylab = "Average residual", 
           main = "Binned residual plot (nothing ", 
           cex.pts = 0.8, 
           col.pts = 1,
           col.int = "gray")
binnedplot(fitted(subset_model2, type = "response"), 
           residuals(subset_model2,type="response"), 
           nclass = NULL, 
           xlab = "Expected Values", 
           ylab = "Average residual", 
           main = "Binned residual plot age spline", 
           cex.pts = 0.8, 
           col.pts = 1,
           col.int = "gray")


binnedplot(x=subset_model1$model$age2, 
           residuals(subset_model3), 
           nclass = NULL, 
           xlab = "Expected Values", 
           ylab = "Average residual", 
           main = "Binned residual plot (Subset of First 50,0000 Observations)", 
           cex.pts = 0.8, 
           col.pts = 1, 
           col.int = "gray")






summary(model1)
drop1(model1,test="LRT")

#三、I开始画残差图和模型诊断

forceDiag<-transform(
  d2,
  .fitted=predict(model1,type = "response"),
  .deviance=residuals(model1),
  .pearson=residuals(model1, type = "pearson")
)

##画原始残差图和拟合值
p1<-qplot(.fitted, .deviance,data=forceDiag,geom="hex")+
  geom_smooth(size=1)+
  xlab("fitted values")+ylab("deviance residuals")
p2<-qplot(.fitted,.pearson,data=forceDiag,geom="hex")+
  geom_smooth(size=1)+
  xlab("fitted values")+ylab("pearson residuals")
p3<-qplot(.fitted,sqrt(abs(.pearson)),data=forceDiag,geom="hex")+
  geom_smooth(size=1)+
  xlab("fitted values")+ylab(expression(sqrt("|Deviance Residuals|")))
grid.arrange(p1,p2,p3,ncol=3)

#画针对logistic回归的残差bin图
# Assuming 'model1' is your logistic regression model
# You may need to replace 'model1' with the actual name of your model

# Subset the first 10,000 observations
# Assuming model1 is your original model
# nrow(model1$model) gets the number of rows in the model's dataset

#确定某一年的残差图
subset_size <- 100000
subset_indices <- seq_len(min(subset_size, nrow(model1$model)))



# Create random indices for the subset
#subset_indices <- sample(1:nrow(model1$model), size = min(subset_size, nrow(model1$model)))

model1<-glm(form,family = binomial(link = "logit"),data=d3)
model3<-glm(form3,family = binomial(link = "logit"),data=d3)

# Create a binned residual plot
# Create a binned residual plot
binnedplot(fitted(subset_model1), 
           residuals(subset_model1, type = "response"), 
           nclass = NULL, 
           xlab = "Expected Values", 
           ylab = "Average residual", 
           main = "Binned residual plot (Subset of First 50,0000 Observations)", 
           cex.pts = 0.8, 
           col.pts = 1,
           col.int = "gray")


# II 画残差bin图和连续型变量之间的关系
binnedplot(x=subset_model1$model$age2, 
           residuals(subset_model1, type = "response"), 
           nclass = NULL, 
           xlab = "Expected Values", 
           ylab = "Average residual", 
           main = "Binned r
           esidual plot (Subset of First 50,0000 Observations)", 
           cex.pts = 0.8, 
           col.pts = 1, 
           col.int = "gray")

binnedplot(x=subset_model1$model$cs, 
           residuals(subset_model1, type = "response"), 
           nclass = NULL, 
           xlab = "Expected Values", 
           ylab = "Average residual", 
           main = "Binned residual plot (Subset of First 50,0000 Observations)", 
           cex.pts = 0.8, 
           col.pts = 1, 
           col.int = "gray")


#画残差图和连续型变量之间的关系，决定是否进行线性变换
q1<- qplot(age2, .deviance,data=forceDiag,geom="hex")+
  geom_smooth(method="lm",size=1)+
  xlab("age2")+ylab("deviance residuals")
print(q1)


q11<- qplot(age2, .deviance,data=forceDiag3,geom="hex")+
  geom_smooth(method="lm",size=1)+
  xlab("age2")+ylab("deviance residuals")
print(q11)

q2<- qplot(cs, .deviance,data=forceDiag,geom="hex")+
  geom_smooth(method="lm",size=1)+
  xlab("cs")+ylab("deviance residuals")
grid.arrange(q1,q2,ncol=2)

#画残差图和因子型变量的关系(似乎用不太上)
q1 <- ggplot(forceDiag, aes(x = gender, y = .deviance)) +
  geom_boxplot(fill = "skyblue", color = "blue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Residuals vs. gender",
       x = "gender",
       y = "Residuals") +
  theme_minimal()


#计算因子型变量分别对应的残差数值

forceDiag%>%
  group_by(gender)%>%
  summarise(mean_resid=mean(.deviance))

forceDiag%>%
  group_by(race2)%>%
  summarise(mean_resid=mean(.deviance))

forceDiag%>%
  group_by(daytime)%>%
  summarise(mean_resid=mean(.deviance))


forceDiag%>%
  group_by(inout2)%>%
  summarise(mean_resid=mean(.deviance))

forceDiag%>%
  group_by(ac_incid)%>%
  summarise(mean_resid=mean(.deviance))

forceDiag%>%
  group_by(ac_time)%>%
  summarise(mean_resid=mean(.deviance))

forceDiag%>%
  group_by(offunif2)%>%
  summarise(mean_resid=mean(.deviance))

forceDiag%>%
  group_by(typeofid2)%>%
  summarise(mean_resid=mean(.deviance))


forceDiag%>%
  group_by(othpers2)%>%
  summarise(mean_resid=mean(.deviance))


forceDiag%>%
  group_by(cs)%>%
  summarise(mean_resid=mean(.deviance))


#model fit检验-2 # 绘制 ROC 曲线 AUC=0.63还可以 
predicted_prob <- predict(model3, type = "response")
 # 创建 ROC 曲线对象
roc_curve <- roc(response = d3$force, predictor = predicted_prob)

## Plot ROC curve
plot(roc_curve, col = "blue", main = "ROC Curve", lwd = 2)

# 添加对角线
abline(a = 0, b = 1, col = "red", lwd = 2, lty = 2)
## 添加AUC值
text(0.8, 0.2, paste("AUC =", round(auc(roc_curve), 2)), col = "blue", cex = 1.2)



##一个新的模型
#
form3 <- force~race2+gender+age2+inout2+ac_incid*ac_time+offunif2+typeofid2+othpers2+cs
model3<-glm(form3,family = binomial(link = "logit"),data=d2)
model4<-glm(form4,family = binomi(link="logit"),data=d3)

forceDiag3<-transform(
  d2,
  .fitted=predict(model3,type = "response"),
  .deviance=residuals(model3),
  .pearson=residuals(model3, type = "pearson")
)

forceDiag4<-transform(
  d2,
  .fitted=predict(model4,type = "response"),
  .deviance=residuals(model4),
  .pearson=residuals(model4, type = "pearson")
)


forceDiag3%>%
  group_by(ac_incid)%>%
  summarise(mean_resid=mean(.deviance))

forceDiag3%>%
  group_by(ac_time)%>%
  summarise(mean_resid=mean(.deviance))



#共线性、变量选择和模型构建、交互项问题、模型共线性问题不大
library("car")
vif_values <- vif(model1)

#III   bootstrp 方法
#不用bootstrap算置信区间（两种方法，书上p225）
confint.default(model1,"age2")
confint(model1,"age2")
#nonparametric bootstrap
B<- 5
n<- nrow(d2)
beta<- numeric(B)
for(b in 1:B){
  i<- sample(n,n,replace=TRUE)
  bootGlm<-glm(form,family = binomial(link = "logit"), data = d2[i,])
  beta[b]<-coefficients(bootGlm)["age2"]
}
#parametric bootsrtap
parbeta<-numeric(B)
d2Samp<-d2
for(b in 1:B){
  d2$force<-simulate(model1)[,1]
  bootGlm<-glm(form,family = binomial(link = "logit"), data = d2Samp)
  parbeta[b]<-coefficients(bootGlm)["age2"]
}

sebeta<-sd(beta)
separbeta<-sd(parbeta)
sebeta
separbeta
#standard error based on ordinary analytic approximations
coef(summary(model1))["age2",2]
betahat<-coefficients(model1)["age2"]
#使用非参数方法的置信区间，解释直接抄书上p226-227
betahat+1.96*sebeta*c(-1,1)
#使用参数方法的置信区间
betahat+1.96*separbeta*c(-1,1)



# 查看转换后的数据集




# Create a binned residual plot
arm::binnedplot(fitted(subset_model1,type = "response"), 
                residuals(subset_model1), 
                nclass = NULL,             
                xlab = "Expected Values", 
                ylab = "Average residual", 
                main = "Binned residual plot (Subset of First 100000 Observations)", 
                cex.pts = 0.8, 
                col.pts = 1, 
                col.int = "gray")

