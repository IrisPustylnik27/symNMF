import sys

import numpy as np
from sklearn.metrics import adjusted_rand_score, silhouette_score

import symnmf


def main():
    """Run both clustering methods and print their silhouette and ARI scores."""
    if len(sys.argv) != 3:
        exception_handler()

    try:
        k = int(sys.argv[1])
        file_name = sys.argv[2]
    except Exception:
        exception_handler()

    try:
        X = np.loadtxt(file_name, delimiter=",", ndmin=2)
    except Exception:
        exception_handler()

    kmeans_labels = get_kmeans_labels(X, k)

    try:
        H_matrix = symnmf.fit(X, k)
    except Exception:
        exception_handler()
    symnmf_labels = np.argmax(H_matrix, axis=1)

    if len(np.unique(symnmf_labels)) < 2 or len(np.unique(kmeans_labels)) < 2:
        exception_handler()
    try:
        sil_kmeans = silhouette_score(X, kmeans_labels)
        sil_symnmf = silhouette_score(X, symnmf_labels)
        ari_score = adjusted_rand_score(symnmf_labels, kmeans_labels)
    except Exception:
        exception_handler()

    print(f"nmf: {sil_symnmf:.4f}")
    print(f"kmeans: {sil_kmeans:.4f}")
    print(f"ari: {ari_score:.4f}")


def exception_handler():
    """Print the standard error message and terminate with a failure status."""
    print("An Error Has Occurred")
    sys.exit(1)


# Copied the K-means implementation from HW1 and modified it to work
# with NumPy arrays instead of reading from standard input.
class Cluster(object):
    """Store a centroid and the points currently assigned to its cluster."""

    def __init__(self, centroid):
        """Initialize a cluster with a centroid and no assigned points."""
        self.centroid = centroid
        self.points = []

    def update_centroid(self):
        """Replace the centroid with the mean of its assigned points, if any."""
        if not self.points:
            return
        total = Vector([0] * len(self.points[0]))
        for point in self.points:
            total += point
        self.centroid = total / len(self.points)


class Vector(object):
    """Represent a numeric vector used by the K-means implementation."""

    def __init__(self, values):
        """Initialize a vector from its coordinate values."""
        self.values = values

    def __add__(self, other):
        """Return the coordinate-wise sum of this vector and another vector."""
        return Vector([a + b for a, b in zip(self.values, other.values)])

    def __truediv__(self, scalar):
        """Return this vector divided coordinate-wise by a scalar."""
        return Vector([a / scalar for a in self.values])

    def __len__(self):
        """Return the number of coordinates in the vector."""
        return len(self.values)


def euclidean_distance(vec1, vec2):
    """Return the Euclidean distance between two vectors."""
    return sum((a - b) ** 2 for a, b in zip(vec1.values, vec2.values)) ** 0.5


def get_kmeans_labels(data_matrix, k):
    """Cluster the data with HW1 K-means and return one label per row."""
    epsilon = 1e-4
    max_iter = 300
    counter = 0
    is_converged = False
    clusters = []
    vectors = [Vector(row.tolist()) for row in data_matrix]

    for i in range(k):
        centroid = Vector(vectors[i].values.copy())
        clusters.append(Cluster(centroid))

    while not is_converged and max_iter > counter:
        for cluster in clusters:
            cluster.points = []

        for vector in vectors:
            closest_cluster = min(
                clusters,
                key=lambda c: euclidean_distance(vector, c.centroid),
            )
            closest_cluster.points.append(vector)

        is_converged = True
        for cluster in clusters:
            old_centroid = cluster.centroid
            cluster.update_centroid()
            delta_centroid = euclidean_distance(old_centroid, cluster.centroid)
            if delta_centroid >= epsilon:
                is_converged = False
        counter += 1

    labels = []
    for vector in vectors:
        distances = [euclidean_distance(vector, c.centroid) for c in clusters]
        closest_idx = distances.index(min(distances))
        labels.append(closest_idx)

    return labels


if __name__ == "__main__":
    main()
