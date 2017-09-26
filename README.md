# Fuzzy C Mean Clustering on IRIS Dataset

This is the [Fuzzy C Mean Clustering algorithm](https://home.deib.polimi.it/matteucc/Clustering/tutorial_html/cmeans.html) implemented in C, and used over [IRIS dataset](https://archive.ics.uci.edu/ml/datasets/iris).


Clustering accuracy achieved: 92.67% (approx).


Instead of initializing the membership values (Uij) I initialized the cluster centroids based on the initialization method described [here](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.451.7988&rep=rep1&type=pdf) (Even though this paper is discusses a centroid initialization method for K-means algorithm, not Fuzzy-C-Means).


I then calculated the membership values from the cluster centroids following the standard Fuzzy-C-Means equations. The rest of the program follows the standard norms of Fuzzy-C-Means. 


