##空间转录组基于图像颜色RGB聚类####


img <- brain@images$slice1@image
imgDm <- dim(img)
# Assign RGB channels to data frame
imgRGB <- data.frame(
  x = rep(1:imgDm[2], each = imgDm[1]),
  y = rep(imgDm[1]:1, imgDm[2]),
  R = as.vector(img[,,1]),
  G = as.vector(img[,,2]),
  B = as.vector(img[,,3])
)
#kmean
kClusters <- 5
kMeans <- kmeans(imgRGB[, c("R", "G", "B")], centers = kClusters)
kColours <- rgb(kMeans$centers[kMeans$cluster,])
unicolor <- unique(kColours)
print(unicolor)

# plot cluster result
p5 <- list()
for (i in 1:length(unicolor)){
  binary_col <- kColours
  #set other clusters as background color.
  binary_col[kColours %in% unicolor[-i]] <- "#FAF9F9"
  p5[[i]] <- ggplot(data = imgRGB, aes(x = x, y = y)) + 
    geom_point(colour = binary_col) +
    labs(title = paste("k-Means Clustering: ",i,'; Color: ',unicolor[i])) +
    xlab("x") +
    ylab("y") ## + plotTheme()
  }

print(p5[[1]])


coordinates <- GetTissueCoordinates(object = brain)
coordinates$imagerow <- 600 - coordinates$imagerow
ggplot(coordinates,aes(imagecol,imagerow)) + geom_point()




imgRGB$kmean5 <- kMeans$cluster
posr <- round(coordinates,0)
colnames(posr) <- c('y','x')
posr$spot <- rownames(posr)
image_meta <- inner_join(imgRGB, posr, by = c('x','y'))
rownames(image_meta) <- image_meta$spot
image_meta <- image_meta[names(Idents(brain)),]
image_meta$seurat_clusters <- as.vector(Idents(brain))
brain@meta.data <- cbind(brain@meta.data, image_meta)
SpatialDimPlot(brain, label = TRUE, group.by = 'kmean5', label.size = 5)
