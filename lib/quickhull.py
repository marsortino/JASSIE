import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
import numpy as np

def convex_hull(points):
    """
    Computes the convex hull of a given list of points.
    Returns index of the points belonging to the convex hull.

    Parameters:
    --
    pointlist :class:`~list`: list of blobs coordinates.
    """
    hull = ConvexHull(points)
    simplex = np.unique(hull.simplices)
    return simplex

def convex_hull_vertices(points):
    """
    Computes the convex hull of a given list of points.
    Returns the vertices of the points belonging to the convex hull.

    Parameters:
    --
    pointlist :class:`~list`: list of blobs coordinates.
    """
    hull = ConvexHull(points)
    return hull.points

def plot_convex_hull(points):
    """
    Computes the convex hull of a given list of points and plots the result.

    Parameters:
    --
    pointlist :class:`~list`: list of blobs coordinates.
    """
    hull = ConvexHull(points)
    status = 0
    for index in hull.vertices:
        if status == 0:
            point_list = np.array(hull.points[index])
            status += 1
        else:
            point_list = np.vstack([point_list, hull.points[index]])

    # Plot the original points
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(points[:, 0], points[:, 1], points[:, 2], c='blue', marker='o', label='Original Points')

    #Plot the convex hull
    for simplex in hull.simplices:
        simplex = np.append(simplex, simplex[0])  # Close the loop
        ax.plot(points[simplex, 0], points[simplex, 1], points[simplex, 2], 'r-')
    
    simplex = np.unique(hull.simplices)

    ax.set_xlabel('X-axis')
    ax.set_ylabel('Y-axis')
    ax.set_zlabel('Z-axis')
    ax.legend()
    plt.title('Convex Hull Plot')
    plt.show()
    return point_list, simplex

