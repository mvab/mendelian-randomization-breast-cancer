# install.packages("devtools")
devtools::install_github("malcolmbarrett/ggdag")

library(ggdag)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
theme_set(theme_dag_blank())


dagify(y ~ x) %>% 
  ggdag()

#You also sometimes see edges that look bi-directed, like this:
  dagify(y ~~ x) %>% 
  ggdag()
  
  #  canonicalize the DAG: Add the latent variable in to the graph
  dagify(y ~~ x) %>% 
    ggdag_canonical()   
  
#A DAG is also acyclic, which means that there are no feedback loops;
# a variable can’t be its own descendant. The above are all DAGs because they are acyclic, but this is not:
    
dagify(y ~ x,
       x ~ a,
       a ~ y) %>% 
    ggdag()   


#Relationships between variables
#Let’s say we’re looking at the relationship between smoking and cardiac arrest.
#We might assume that smoking causes changes in cholesterol, which causes cardiac arrest:
  
smoking_ca_dag <- dagify(cardiacarrest ~ cholesterol,
                           cholesterol ~ smoking + weight,
                           smoking ~ unhealthy,
                           weight ~ unhealthy,
                           labels = c("cardiacarrest" = "Cardiac\n Arrest", 
                                      "smoking" = "Smoking",
                                      "cholesterol" = "Cholesterol",
                                      "unhealthy" = "Unhealthy\n Lifestyle",
                                      "weight" = "Weight"),
                           latent = "unhealthy",
                           exposure = "smoking",
                           outcome = "cardiacarrest")

ggdag(smoking_ca_dag, text = FALSE, use_labels = "label")


earlyBMI_BC <- dagify(igf ~ earlyBMI,
                      BC ~ igf,
                      BC ~ earlyBMI,
                      labels = c(
                        'igf' = "IGF1",
                        'earlyBMI' = "Childhood BMI",
                        'BC' = "Breast cancer"),
                      exposure = 'earlyBMI',
                      outcome = 'BC',
                      latent = 'igf')
                        
ggdag(earlyBMI_BC, text = FALSE, use_labels = "label", stylized=TRUE, 
      edge_type="diagonal")  +
  ggtitle("hello")

pal<-rev(unname(yarrr::piratepal("info2")))[3]
scales::show_col(pal)
ggdag_drelationship(earlyBMI_BC, controlling_for = "igf", 
                 text = FALSE, use_labels = "label",
                 stylized=TRUE, 
                 edge_type="diagonal")  +
  scale_color_manual(values = pal, na.value = "#F08892FF")+
  scale_fill_manual(values = pal, na.value = "#F08892FF")+
  theme(legend.position = "none")





earlyBMI_BC_MVMR <- dagify(igf ~ earlyBMI,
                      BC ~ igf,
                      BC ~ earlyBMI,
                      igf~SNPs, 
                      earlyBMI ~ SNPs,
                      labels = c(
                        'igf' = "IGF1",
                        'earlyBMI' = "Childhood BMI",
                        'BC' = "Breast cancer", 
                        "SNPs" = "IVs"),
                      exposure = 'SNPs',
                      outcome = 'BC')
ggdag(earlyBMI_BC_MVMR, text = FALSE, use_labels = "label", stylized=TRUE, 
      edge_type="link")  +
  ggtitle("hello")

ggdag_drelationship(earlyBMI_BC_MVMR, controlling_for = "igf", 
                    text = FALSE, use_labels = "label",
                    stylized=TRUE, 
                    edge_type="diagonal")  +
  #scale_color_brewer(palette = "Dark2")+
  theme(legend.position = "none")



#### other way

# create DAG object
earlyBMI_BC3 <- dagify(igf ~ earlyBMI,
                      BC ~ igf + earlyBMI)


# tidy the dag object and supply to ggplot
set.seed(100)
earlyBMI_BC3 %>%
  tidy_dagitty() %>%
 # x and y are the order on the plane in those directions
  mutate(x = c(0, 0, 1, 2)) %>%
  mutate(y = c(0, 0, 2, 0)) %>%
  # these are order "coordinates" where you want the `name` to point
  mutate(xend = c(2, 1, 2, NA)) %>%
  mutate(yend = c(0, 2, 0, NA)) %>%

  dag_label(labels = c(
    'igf' = "IGF1",
    'earlyBMI' = "Childhood BMI",
    'BC' = "Breast cancer")
  ) %>%
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_edges(
    edge_colour = "pink",
    edge_width = .8
  ) +
  geom_dag_node(
    color = "#F08892FF",
    alpha = 0.8
  ) +
  #geom_dag_text(color = "black") +
  geom_dag_label_repel(aes(label = label),
                       col = "white",
                       label.size = .4,
                       fill = "#20a486ff",
                       alpha = 0.8,
                       show.legend = FALSE,
                       nudge_x = .7,
                       nudge_y = .3
  ) +
  labs(
    title = " Directed Acyclic Graph",
    subtitle = " Two Variables of Interest with a Mediator"
  ) +
  xlim(c(-1.5, 3.5)) +
  ylim(c(-.33, 2.2))+
  geom_rect(
    xmin = -.5,
    xmax = 3.25,
    ymin = -.25,
    ymax = .65,
    alpha = .04,
    fill = "grey"
  ) 






