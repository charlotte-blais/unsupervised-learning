#Unsupervised Learning Examples 

library(ISLR)
data(USArrests)

# Hierarchical Clustering with Euclidean Dissimilarity Measure and Complete linkage: 
dat.scaled = scale(USArrests)
euclidean.dist = hclust(dist(dat.scaled), method="complete")
plot(euclidean.dist,main="Hierarchical Clustering with Scaled Features")
# Cut dendrogram to get 3 clusters
euc.dist.cut = cutree(euclidean.dist, 3)
# Now try without scaling: 
euclidean.dist.nsc = hclust(dist(USArrests), method="complete")
euc.dist.cut.nsc = cutree(euclidean.dist.nsc, 3)
# Compare results
# We can see, as we would expect, that scaling has an impact on the clustering
table(euc.dist.cut.nsc, euc.dist.cut)
summary(USArrests)
# UrbanPop measures the percent of the state's population that lives in an urban area and 
# the other features represent the number of instances of a crime per 100,000 people. Since these scales
# measure different quantities we will want to scale the features to have standard deviation = 1 so that population 
# doesn't dominate the analysis. Another reason to scale the standard deviations to 1 is that assault occurs much more 
# frequently than rape or murder and therefore will likely have a larger effect on the clustering. But, this may not be desirable
# because the gravity of the three offenses may not be considered equal. For these reasons we will scale the mean to 0 and 
# standard deviation to 1 in this case. 



# Short Example on Simulated Data using Principal Component Analysis (PCA)
# Generate simulated data to perform PCA and K-means clustering on 
x = matrix(rnorm(60*50), ncol=50)
# Make a mean shift to differentiate the 3 classes:
x[1:20,5] = x[1:20,5]+5
x[21:40,c(6,8)]=x[21:40,c(6,8)]-8
# even though the variables are drawn 
pca = prcomp(x)
# Plot the first two principal components (and as expected, because this is simulated data, we see separation): 
# Note: the pca object holds the principal component score vectors in pca$x. Equivalently we could matrix multiply the X matrix by 
# pr.out$rotation (the loading vectors) to get the coordinates of principal component scores (which are the coordinates of the data in the rotated 
# coordinate system so we can see the plot of PC1 vs. PC2) (plotted here)
par(mfrow=c(1,2))
biplot(pca, scale=0)
plot(pca$x)

# Let's see what K-means clustering looks like with k=2. We will try with 25 initial configurations.
kmthree=kmeans(x,2,nstart=25)
# We try to minimize the total within-cluster sum of squares for k-means. This is computed as..
# result: 3192.161 (this is high)
kmthree$tot.withinss
# Let's try with k=4
kmfour=kmeans(x,4,nstart=25)
# result: 2786.061 (this is high)
kmfour$tot.withinss

# PCA is a way of reducing dimensionality in a dataset. K-means tends to do poorly in high dimensional datasets because observations can be far away from each other, 
# especially when there aren't that many observations in comparison to the number of predictors as in this situation.  
# A good technique can be to lower the dimensionality of the problem using PCA and then perform K-means on the principal components.
# Let's try performing K-means on the first two principal component score vectors to see if this helps
kmpc = kmeans(pca$x[,1:2],3,nstart=25)
# result: 118.5042 (much better than before). 
kmpc$tot.withinss
# We can look at the proportion of variance explained by each principal component and cumulatively:
# Generally principal components do best when the first few principal components explain the majority of the variance
# In this case we can see from the plots below that the first principal component alone explains nearly 40% of the data - a good sign for PCA. 
pcvar = pca$sdev^2
pve = pcvar/sum(pcvar)
# note: type="b" in plot indicates that points should be joined by lines.
par(mfrow=c(1,2))
plot(pve, xlab="Principal Component", ylab="Proportion of Variance Explained", ylim=c(0,1), type="b")
plot(cumsum(pve),xlab="Principal Component", ylab="Cumulative PVE", ylim=c(0,1),type="b")

# Let's try scaling x before we do k-means clustering
xscale <-scale(x, scale=TRUE)
kmscaled = kmeans(xscale,3,nstart=25)
# As we can see scaling does not improve the total within-cluster variation much probably because the same variable was increased by the same amount 
# for all instances of a particular class and only one class had more than one variable increased. So euclidean distances within a cluster likely either 
# did not change or changed marginally. Also, all observations were drawn from the same distribution on the same scale originally anyway.
#result: 2638.427
kmscaled$tot.withinss
# Compare to the PCA results (some overlap)
biplot(pca, col=kmscaled$cluster)

# Comparing hierarchical clustering results with different types of linkage used 
# Use correlation-based distance
# Add path in lieu of ellipsis 
genes = as.matrix(read.csv(".../GeneExpression.csv",header=F))
genes.dist = as.dist(1-cor(genes))
par(mfrow=c(1,3))
# We can see that in this example (as usually happens) single linkage results in a more trailing dendrogram
# and complete and average linkage dendrograms are more balanced 
plot(hclust(genes.dist, method="complete"),main="Complete Linkage")
plot(hclust(genes.dist, method="single"),main="Single Linkage")
plot(hclust(genes.dist, method="average"),main="Average Linkage")

