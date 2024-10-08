import numpy as np
import subprocess
import os
import time
import sys

class raytracing:
    """
    Calculates the intersection length of lines originating from a specified source point and passing 
    through the centers of a list of blobs. It determines whether these lines intersect with other blobs and,
    if so, provides the length of the intersection between the line and the other blobs. 
    """

    def __init__(self, blobslist, observer):
        """
        Parameters:
        --
        blobslist :class:`~list` list of objects.block; 
        observer :class:`~numpy.array` array of the observer in cartesian coordinates. They should be converted as u.cm beforehand (although it needs to be passed stripped of .unity)
        """
        self.ray_path(blobslist, observer)

    def distance_between_line_and_point(self, p, q, r):
        """
        Returns the distance between a line, defined passing through two points, and a point. 
        Computes the direction vector D of the line:
            D = p-q   
        Computes the vector from point1 to point3.
            PR = r-p
        It then computes the projection of the vector on the direction vector:
            projection = PR x D/ D x D * D
        where 'x' stands for dot product 
        Finally it normalize the distance of PR-projection
        """

        P = np.array(p)
        Q = np.array(q)
        R = np.array(r)
        
        D = self.normal_vec(P,Q)
        
        # Vector from point P to point R
        PR = R - P
        
        # Calculate the projection of PR onto D
        projection = np.dot(PR, D) / np.dot(D, D) * D
        
        # Calculate the distance between point R and the projected point
        distance = np.linalg.norm(PR - projection)
        
        return distance

    def chord(self, radius, distance):
        """
        Computes the chord defined by the intersection between the line and the sphere 
        """
        return 2*np.sqrt(2*radius*distance-distance**2)

    def interacting_region(self, source, point_on_plane, normal):
        """
        Computes the plane passing through the center of the i-th blob and orthogonal to the line passing
        through the observer and the i-th blob center.
        """

        distance = np.dot(source - point_on_plane, normal)
        if distance == 0:
            return False
        sgn_distance = distance/np.abs(distance)
        return sgn_distance

    def normal_vec(self, point1, point2):
        """
        Computes the directional vector passing through two points.
        """
        P = np.array(point1)
        Q = np.array(point2)  
        return Q-P



    def ray_path(self, blobslist, obs):
        """
        Computes the path of a photon.
        """

        ray_flux = []
        for blob_i in blobslist:
            counter = 0
            raypath = 0*blob_i.radius.unit
            normal = self.normal_vec(obs.coords, blob_i.center)
            obs_sgn = self.interacting_region(obs.coords, blob_i.center, normal)

            for blob_j in blobslist:
                if blob_i != blob_j:
                    blob_j_sgn_1 = self.interacting_region(blob_j.center, blob_i.center, normal)
                    blob_j_sgn_2 = self.interacting_region(blob_j.center, obs.coords, -normal)
                    if (obs_sgn == blob_j_sgn_1) & (obs_sgn == blob_j_sgn_2):
                        value = self.distance_between_line_and_point(blob_i.center, obs.coords, blob_j.center)
                        if value == 0:

                            counter += 1
                            raypath += blob_j.radius
                        elif (value > 0) & (value < blob_j.radius.value):

                            counter +=1
                            raypath += self.chord(blob_j.radius.value, value)*blob_i.radius.unit

            # It now stores the value found in each blob with the proper tag. (source, observer)
            blob_i.rayline(raypath, obs.name)
        return ray_flux

def C_raytracing(target_list, blocklist, filename):
    """
    Computes raytracing using a C algorithm.
    
    target :class:`~tuple` of objects.target
    blobslist :class:`~list` of blocks
    filename :class:`~str` name of the simulation input file.
    """

    radius_unit = blocklist[0].radius.unit

    #output_file_name =  'tmp_' + str(target.name) + '_raytrace_output.txt'
    script_dir = os.path.dirname(__file__)
    tmp_file_path = 'tmp_files/tmp_list_value_' + filename + '.txt'
    tmp_file_path = os.path.join(script_dir, tmp_file_path)

    c_script_path = "c_routines/exec/raytrace"
    abs_file_path = os.path.join(script_dir, c_script_path)
    for target in target_list:
        result = subprocess.run([abs_file_path, str(target.coords[0]), str(target.coords[1]), str(target.coords[2]), tmp_file_path], capture_output=True, text = True, check=True)
        #print(result.stdout)

        output_path = os.path.join(script_dir, 'tmp_files/tmp_raytrace_o_' + filename + '.txt')
        values = open(output_path, 'r').readlines()
        if target.name == 'source':
            for value, block in zip(values, blocklist):
                block.rayline(float((value.split('\n')[0]))*radius_unit, target.name)

        #Debugging purpose
        #hist_indices = []
        #total_encounters = 0
        if target.name == 'observer':
            output_path_index = os.path.join(script_dir, 'tmp_files/tmp_indices_' + filename + '.txt')
            values_index = open(output_path_index, 'r').readlines()
            for value, line, block in zip(values, values_index, blocklist):
                array = [int(x) for x in line.strip().split(",")]
                # Debugging purpose
                # if array[0] != -1:
                #     hist_indices.append(len(array))
                #     for number in array:
                #         total_encounters += number
                block.rayline(float((value.split('\n')[0]))*radius_unit, target.name, order_array = array)

            # Deleting tmp file
            output_path_index = 'rm ' + output_path_index
            subprocess.run([output_path_index], capture_output= True, shell=True)

            # debugging purpose

            # import matplotlib.pyplot as plt
            # n_bin = round(np.power(len(hist_indices),1/2))
            # plt.hist(hist_indices, bins = n_bin)
            # print(hist_indices)
            # plt.xlabel('Valore')
            # plt.ylabel('Frequenza')
            # plt.title('total_'+ str(total_encounters)+'_'+filename)
            # plt.savefig('hist_'+filename+'.png')
            # plt.close()

    # Deleting tmp raytrace input file
    tmp_file_path = 'rm ' + tmp_file_path
    subprocess.run([tmp_file_path], capture_output= True, shell=True)


    # Deleting tmp raytrace.c output file.
    output_path = 'rm ' + output_path
    subprocess.run([output_path], capture_output= True, shell=True)


