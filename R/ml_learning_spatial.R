##################################################################
######supervised learning & protein localization prediction#######
##################################################################

# x : data frame where rows are proteins and columns are fractions
# model: string specifying which classification model to use
# organelle:organelle-specific protein distribution of desired organelle (e.g., Mitochondria)




ml_learning_spatial <- 
    function(x, model,organelle) {
        
    library(gplots)
    library(ggplot2)
    library(data.table)
    library(dplyr)
    library(tidyr)
    library(reshape2)
    library(caret)
    library(magrittr)
        
    df <- x
    #preparer the training data 
    x_train <- 
        df %>%
        as.data.frame(.) %>%
        set_rownames(.$Gene) %>%
        dplyr::select(-1)%>%
        filter(Marker != "unknown") %>%
        mutate_at(vars(Marker), 
            list(factor))
    
    x_train$Marker <- 
        factor(x_train$Marker, 
            labels = make.names(levels(x_train$Marker)))
    
    
    # define training control using k-fold Cross Validation
    train_control <- 
        trainControl(method="cv",
            classProbs=TRUE,
            number=10, verbose = T,savePred=T)
    
    #train the classifier (e.g, random forest)
    for (ml in model) {
        message(ml)
    set.seed(100)
    rfModel <- 
        train(Marker~., data=x_train,
            trControl=train_control,
            method=ml,
            preProcess = c("center", "scale"),
            metric= "Accuracy")
    }
    
    
    #Predict on the whole data 
    predictdf <-
        df %>%
        as.data.frame(.) %>%
        set_rownames(.$Gene) %>%
        dplyr::select(-c(1,8))
    
    set.seed(101)
    predictdf  <-
        predict(rfModel,predictdf, type = "prob") 
    
##################################################################
#################### visulaize the result
##################################################################    
        
    #organelle-specific score distributions
    datBOX <- 
        predictdf %>%
        gather("Compartment", "probScore", 1:9) %>%
        group_by(Compartment) %>%
        filter(probScore > 0) %>%
        mutate_at(vars(Compartment), 
            list(factor))
    
    
    print(ggplot(datBOX, aes(Compartment, probScore)) +
        geom_boxplot(aes(fill=Compartment)) +
        geom_hline(yintercept=0.5, linetype="dashed", #0.50 cutoff
            color = "red", size=1) +
        theme_bw()) 
    
    
    #visualize the result on PCA
    finaldf <- #predicted scores for each compartment
        predictdf %>%
        tibble::rownames_to_column('Gene') %>%
        gather("Predicted_Compartment", "probScore", 2:10) %>%
        group_by(Gene) %>%
        dplyr::slice(which.max(probScore))  %>%
        left_join(.,df) %>%
        dplyr::select(-Marker) %>%
        ungroup()
    
    
    pca_res <- #keep the numeric data
        prcomp(finaldf[,-c(1,2)], scale. = TRUE)
    print(autoplot(pca_res,
        data = finaldf[,-c(1)], 
        colour = "Predicted_Compartment") +
        theme_bw())
    
    
    
    #organelle-specific protein distributions along the gradient fractions
    
    #can select any compartment(e.g., Mitochondria)
    for(cm in organelle){
        message(cm)
        
    orgDist <- 
        finaldf %>%
        filter(Predicted_Compartment == cm) %>%
        dplyr::select(c(1,4:9)) %>%
        gather("condition", "intensity", 2:7)
    
    print(ggplot(data=orgDist,
        aes(x=condition, y=intensity, group = Gene)) +
        geom_line(color = "red") +
        theme_bw())
        }
return(predictdf)
    }


