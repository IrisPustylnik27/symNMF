import sys
import numpy as np
from sklearn.metrics import silhouette_score, adjusted_rand_score

import symnmf

def main():
    if len(sys.argv) != 3:
        print("An Error Has Occurred")
        sys.exit(1)
    
    k = int(sys.argv[1])
    file_name = sys.argv[2]

    try:
        X = np.loadtxt(file_name, delimiter=',')
    except Exception:
        print("An Error Has Occurred")
        sys.exit(1)

    kmeans_labels = get_kmeans_labels(X, k)
    H_matrix = symnmf.fit(X, k) 
    symnmf_labels = np.argmax(H_matrix, axis=1)
    
    if len(np.unique(symnmf_labels)) < 2 or len(np.unique(kmeans_labels)) < 2:
        print("An Error Has Occurred")
        sys.exit(1)

    sil_kmeans = silhouette_score(X, kmeans_labels)
    sil_symnmf = silhouette_score(X, symnmf_labels)
    ari_score = adjusted_rand_score(symnmf_labels, kmeans_labels)

    print(f"nmf: {sil_symnmf:.4f}")
    print(f"kmeans: {sil_kmeans:.4f}")
    print(f"ari: {ari_score:.4f}")


#copied the Kmeans from HW1 and modified it to work with numpy arrays instead of reading from stdin
class Cluster(object):

    def __init__(self, centroid):
        self.centroid = centroid
        self.points = []

    def update_centroid(self):
        if not self.points:
            return
        total = Vector([0] * len(self.points[0]))
        for point in self.points:
            total += point
        self.centroid = total / len(self.points)

class Vector(object):

    def __init__(self, values):
        self.values = values

    def __add__(self, other):
        return Vector([a + b for a, b in zip(self.values, other.values)])

    def __truediv__(self, scalar):
        return Vector([a / scalar for a in self.values])

    def __len__(self):
        return len(self.values)
    

def euclidean_distance(vec1, vec2):
    return sum((a - b) ** 2 for a, b in zip(vec1.values, vec2.values)) ** 0.5

def get_kmeans_labels(data_matrix, k):
    epsilon = 1e-4
    maxIter = 300
    counter = 0
    isConverged = False
    clusters = []
    vectors = [Vector(row.tolist()) for row in data_matrix]
    
    # Initialize centroids 
    for i in range(k):
        centroid = Vector(vectors[i].values.copy())
        clusters.append(Cluster(centroid))
    
    while not isConverged and maxIter > counter:
        for cluster in clusters:
            cluster.points = []

        for vector in vectors:
            closest_cluster = min(clusters, key=lambda c: euclidean_distance(vector, c.centroid))
            closest_cluster.points.append(vector)

        isConverged = True

        for cluster in clusters:
            old_centroid = cluster.centroid
            cluster.update_centroid()
            deltaCentroid = euclidean_distance(old_centroid, cluster.centroid)
            if (deltaCentroid >= epsilon):
                isConverged = False
        counter += 1
    
    labels = []
    for vector in vectors:
        distances = [euclidean_distance(vector, c.centroid) for c in clusters]
        closest_idx = distances.index(min(distances))
        labels.append(closest_idx)
        
    return labels

if __name__ == "__main__":
    main()
    