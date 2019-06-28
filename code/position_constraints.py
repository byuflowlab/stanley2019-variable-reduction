import numpy as np
from scipy.spatial import ConvexHull
import sys
sys.dont_write_bytecode = True

def turbineSpacingSquared(turbineX, turbineY):
    """
    Calculates inter-turbine spacing for all turbine pairs
    """
    nTurbines = len(turbineX)
    separation_squared = np.zeros((nTurbines-1)*nTurbines/2)

    k = 0
    for i in range(0, nTurbines):
        for j in range(i+1, nTurbines):
            separation_squared[k] = (turbineX[j]-turbineX[i])**2+(turbineY[j]-turbineY[i])**2
            k += 1
    return separation_squared


def SpacingConstraint(turbineX, turbineY, rotorDiameter, minSpacing=2.0):
    """
    inter turbine spacing constraint
    """
    nTurbines = len(rotorDiameter)
    separation_squared = turbineSpacingSquared(turbineX, turbineY)
    spacing_con = np.zeros(int((nTurbines-1.)*nTurbines/2.))

    k = 0
    for i in range(0, nTurbines):
        for j in range(i+1, nTurbines):
            spacing_con[k] = separation_squared[k] - (0.5*minSpacing*rotorDiameter[i]+0.5*minSpacing*rotorDiameter[j])**2
            k += 1
    return spacing_con


def circularBoundary(turbineX, turbineY, circle_radius, circle_center=np.array([0.,0.])):
    """
    circular boundary constraint
    """
    nTurbines = len(turbineX)
    circle_boundary_constraint = np.zeros(nTurbines)
    for i in range(nTurbines):
        R = np.sqrt((turbineX[i]-circle_center[0])**2+(turbineY[i]-circle_center[1])**2)
        circle_boundary_constraint[i] = circle_radius-R
    return circle_boundary_constraint


def arbitraryBoundary(turbineX, turbineY, boundaryVertices, boundaryNormals):

    if type(turbineX) == np.float64:
        nTurbines = 1
    else:
        nTurbines = len(turbineX)
    locations = np.zeros([nTurbines, 2])
    if nTurbines == 1:
        locations[0] = np.array([turbineX, turbineY])
    else:
        for i in range(0, nTurbines):
            locations[i] = np.array([turbineX[i], turbineY[i]])

    # calculate distance from each point to each face
    boundaryDistances = calculate_distance(locations, boundaryVertices, boundaryNormals)

    return boundaryDistances


def calculate_boundary(vertices):

    # find the points that actually comprise a convex hull
    hull = ConvexHull(list(vertices))

    # keep only vertices that actually comprise a convex hull and arrange in CCW order
    vertices = vertices[hull.vertices]

    # get the real number of vertices
    nVertices = vertices.shape[0]

    # initialize normals array
    unit_normals = np.zeros([nVertices, 2])

    # determine if point is inside or outside of each face, and distance from each face
    for j in range(0, nVertices):

        # calculate the unit normal vector of the current face (taking points CCW)
        if j < nVertices - 1:  # all but the set of point that close the shape
            normal = np.array([vertices[j+1, 1]-vertices[j, 1],
                               -(vertices[j+1, 0]-vertices[j, 0])])
            unit_normals[j] = normal/np.linalg.norm(normal)
        else:   # the set of points that close the shape
            normal = np.array([vertices[0, 1]-vertices[j, 1],
                               -(vertices[0, 0]-vertices[j, 0])])
            unit_normals[j] = normal/np.linalg.norm(normal)

    return vertices, unit_normals



def calculate_distance(points, vertices, unit_normals, return_bool=False):

    """
    :param points: points that you want to calculate the distance from to the faces of the convex hull
    :param vertices: vertices of the convex hull CCW in order s.t. vertices[i] -> first point of face for
           unit_normals[i]
    :param unit_normals: unit normal vector for each face CCW where vertices[i] is first point of face
    :param return_bool: set to True to return an array of bools where True means the corresponding point
           is inside the hull
    :return face_distace: signed perpendicular distance from each point to each face; + is inside
    :return [inside]: (optional) an array of zeros and ones where 1.0 means the corresponding point is inside the hull
    """

    # print points.shape, vertices.shape, unit_normals.shape

    nPoints = points.shape[0]
    nVertices = vertices.shape[0]

    # initialize array to hold distances from each point to each face
    face_distance = np.zeros([nPoints, nVertices])

    if not return_bool:
        # loop through points and find distance to each face
        for i in range(0, nPoints):

            # determine if point is inside or outside of each face, and distance from each face
            for j in range(0, nVertices):

                # define the vector from the point of interest to the first point of the face
                pa = np.array([vertices[j, 0]-points[i, 0], vertices[j, 1]-points[i, 1]])

                # find perpendicular distance from point to current surface (vector projection)
                d_vec = np.vdot(pa, unit_normals[j])*unit_normals[j]

                # calculate the sign of perpendicular distance from point to current face (+ is inside, - is outside)
                face_distance[i, j] = np.vdot(d_vec, unit_normals[j])

        return face_distance

    else:
	# initialize array to hold boolean indicating whether a point is inside the hull or not
        inside = np.zeros(nPoints)

        # loop through points and find distance to each face
        for i in range(0, nPoints):

            # determine if point is inside or outside of each face, and distance from each face
            for j in range(0, nVertices):

                # define the vector from the point of interest to the first point of the face
                pa = np.array([vertices[j, 0]-points[i, 0], vertices[j, 1]-points[i, 1]])

                # find perpendicular distance from point to current surface (vector projection)
                d_vec = np.vdot(pa, unit_normals[j])*unit_normals[j]

                # calculate the sign of perpendicular distance from point to current face (+ is inside, - is outside)
                face_distance[i, j] = np.vdot(d_vec, unit_normals[j])

            # check if the point is inside the convex hull by checking the sign of the distance
            if np.all(face_distance[i] >= 0):
                inside[i] = 1.0

        return face_distance, inside



if __name__=="__main__":

    # turbineX = np.array([0.,0.,0.])
    # turbineY = np.array([0.,500.,1000.])
    # rotorDiameter = np.array([100.,100.,100.])
    #
    # print SpacingConstraint(turbineX, turbineY, rotorDiameter,minSpacing=2.)
    # print circularBoundary(turbineX, turbineY, 10000.)

    import matplotlib.pyplot as plt
    locations = np.loadtxt('layout_amalia.txt')
    turbineX = locations[:, 0]
    turbineY = locations[:, 1]
    plt.plot(turbineX, turbineY, 'o')
    bounds = arbitraryBoundary(turbineX, turbineY)
    plt.show()
