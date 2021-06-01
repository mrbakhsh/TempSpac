##################################################################
######supervised learning & protein localization prediction#######
##################################################################

# x : data frame where rows are proteins and columns are fractions
# model: string specifying which classification model to use




ml_learning_spatial <- 
    function(x, model) {
        
        library(gplots)
        library(ggplot2)
        library(data.table)
        library(dplyr)
        library(tidyr)
        library(reshape2)
        library(caret)
        library(magrittr)
        
        
        df <- x
        
        if(!is.data.frame(df)){
            stop("Input data must be data.frame")
        }
        
        if(all(colnames(df) != "Marker") == TRUE){
            stop("Marker is absent from the data.frame")
        }
        
        if(all(colnames(df) != "Protein.Name") == TRUE){
            stop("Protein names is absent from the data.frame")
        }
        #preparer the training data 
        x_train <- 
            df %>%
            as.data.frame(.) %>%
            set_rownames(.$Protein.Name) %>%
            dplyr::select(-1)%>%
            filter(Marker != "unknown") %>% #drop the unknown
            mutate_at(vars(Marker), 
                list(factor))
        
        x_train$Marker <- 
            factor(x_train$Marker, 
                labels = make.names(levels(x_train$Marker)))
        
        
        # define training control using k-fold Cross Validation
        train_control <- 
            trainControl(method="repeatedcv",#repeated cross-validation
                number=10, #number of resampling iterations
                repeats = 10, #sets of folds to for repeated cross-validation
                classProbs=TRUE, 
                verbose = T,savePred=T)
        
        #train the classifier (e.g, random forest)
        for (ml in model) {
            message(ml)
            set.seed(100)
            fit <- 
                train(Marker~., data=x_train,
                    trControl=train_control,
                    method=ml,
                    preProcess = c("center", "scale"), # necessary task
                    metric= "Accuracy"
                )
        }
        
        
        #Predict on the whole data 
        predictdf <-
            df %>%
            as.data.frame(.) %>%
            set_rownames(.$Protein.Name) %>%
            dplyr::select(-c(1,8))
        
        set.seed(101)
        predictdf  <-
            predict(fit,predictdf, type = "prob") 
        output <- 
            predictdf %>%
            tibble::rownames_to_column("Protein.Name") %>%
            gather('Organelle', "ProbScore", 2:ncol(.)) %>%
            group_by(Protein.Name) %>%
            dplyr::slice(which.max(ProbScore))
        
        ##################################################################
        #################### visulaize the result
        ##################################################################    
        
        #organelle-specific score distributions
        dat <- 
            predictdf %>%
            gather("Compartment", "probScore", 1:9) %>%
            group_by(Compartment) %>%
            filter(probScore > 0) %>%
            mutate_at(vars(Compartment), 
                list(factor))
        
        
        print(ggplot(dat, aes(Compartment, probScore)) +
                geom_boxplot(aes(fill=Compartment)) +
                geom_hline(yintercept=0.5, linetype="dashed", #0.50 cutoff
                    color = "red", size=1) +
                theme_bw() +
                theme(panel.grid = element_blank())  +
                theme(text = element_text(colour = "black",size=12)) +
                theme(
                    axis.text = element_text(colour = "black", size = 12),
                    axis.ticks.length = unit(0.25, "cm")))
        
        
        #visualize the result on PCA
        finaldf <- #predicted scores for each compartment
            predictdf %>%
            tibble::rownames_to_column('Protein.Name') %>%
            gather("Predicted_Compartment", "probScore", 2:10) %>%
            group_by(Protein.Name) %>%
            dplyr::slice(which.max(probScore))  %>%
            left_join(.,df) %>%
            dplyr::select(-Marker) %>%
            ungroup()
        
        pca_res <- #keep the numeric data
            prcomp(finaldf[,-c(1,2)], scale. = TRUE)
        print(autoplot(pca_res,
            data = finaldf[,-c(1)], size = "probScore",
            colour = "Predicted_Compartment") +
                scale_size_continuous(range = c(0,1)) +
                theme_bw() +
                theme(panel.grid = element_blank())  +
                theme(text = element_text(colour = "black",size=12)) +
                theme(
                    axis.text = element_text(colour = "black", size = 12),
                    axis.ticks.length = unit(0.25, "cm"))) 
       
        #organelle distribution plots
        orgDist <- 
            finaldf %>%
            dplyr::select(c(1,4:9,2)) %>%
            gather("condition", "intensity", 2:7) 
        
        print(ggplot(data=orgDist,
            aes(x=condition, y=intensity, group = Protein.Name, color =Predicted_Compartment )) +
                geom_line() +
                facet_wrap(~ Predicted_Compartment,ncol=2) +
                theme_bw()+
                theme(panel.grid = element_blank())  +
                theme(text = element_text(colour = "black",size=12)) +
                theme(
                    axis.text = element_text(colour = "black", size = 8),
                    axis.ticks.length = unit(0.25, "cm"))) 
        
        return(output)
    }

