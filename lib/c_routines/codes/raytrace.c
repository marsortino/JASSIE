#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

typedef struct{
	double x;
	double y;
	double z;
} vec;
	      

typedef struct{
	vec center;
	double r;
} block;

typedef struct{
	int index;
	double distance;
} block_distance;

double distance_between_line_point(vec p, vec q, vec r);
vec normal_vec(vec source, vec block_i);
double interacting_region(vec source, vec point_on_plane, vec normal);
double chord(double radius, double distance);
double raytracing(int i, int max, block* blocklist, vec source, int* psphere_counter, int* psphere_max_counter, 	block_distance *distance_sorter);
void preset_array(int sphere_counter, int sphere_index, double distance, block_distance *distance_array, double* pdistance_max);
void update_indices(int **indices, int row, int number_of_spheres, block_distance *distance_array);
char* get_file_name(char* filepath);
double compute_distance(block block_j, block block_i);
void sorting_distance(block_distance* distance_array, int sphere_index, long double distance, int max_index);
void check_duplicates(block_distance* distance_sorted, int size);


int main(int argc, char *argv[]){
	if (argc != 5){
		fprintf(stderr, "Usage: %s <observer.x observer.y observer.z tmp_file_path>\n", argv[0]);
		return 1;
	}
	//printf("Starting raytrace procedures...\n");

	clock_t start_t, end_t;
	double total_t;
	start_t = clock();
	vec source;
	int i = 0, j, max, sphere_counter, sphere_counter_max = 0, *psphere_counter, *psphere_counter_max;
	FILE *list_raytrace, *list_index;
	block *blocklist;
	long double *raytrace;
	block_distance *distance_array;

	source.x = strtod(argv[1], NULL);
	source.y = strtod(argv[2], NULL);
	source.z = strtod(argv[3], NULL);

	char file_path[500];
	strcpy(file_path, argv[4]);

	list_raytrace = fopen(file_path, "r");

	fscanf(list_raytrace, "%i", &max);	
	blocklist = malloc(max*sizeof(block));
	while(!feof(list_raytrace)){

		fscanf(list_raytrace, "%le %le %le %le", &blocklist[i].center.x, &blocklist[i].center.y, &blocklist[i].center.z, &blocklist[i].r); 
		i++;
	}

	fclose(list_raytrace);

	raytrace = malloc(max*sizeof(long double));
	char *filename = get_file_name(file_path);
	int str_lenght =  strlen(filename);
	char *outputname_rt = (char *)malloc((str_lenght+29+1) * sizeof(char)); // +15 because tmp_raytrace_o_ are 14 chars.
	strcpy(outputname_rt, "lib/tmp_files/tmp_raytrace_o_");
	strcat(outputname_rt, filename);
	char *outputname_indx = (char *)malloc((str_lenght+26+1) *sizeof(char));
	strcpy(outputname_indx, "lib/tmp_files/tmp_indices_");
	strcat(outputname_indx, filename);

	distance_array = malloc(max*sizeof(block_distance));

	list_raytrace = fopen(outputname_rt, "w");
	list_index = fopen(outputname_indx, "w");
	psphere_counter = &sphere_counter;
	psphere_counter_max = &sphere_counter_max;

	for(i=0; i<max; i++){
		sphere_counter = 0;
		raytrace[i] += raytracing(i, max, blocklist, source, psphere_counter, psphere_counter_max, distance_array);
		for(j = 0; j<sphere_counter-1; j++){
			fprintf(list_index, "%d, ", distance_array[j].index);
		}
		if(sphere_counter == 0){
				fprintf(list_index, "%d\n", -1);
		}
		else fprintf(list_index, "%d\n", distance_array[sphere_counter-1].index);
		
		fprintf(list_raytrace, "%0.15Le\n", raytrace[i]);
	}
	free(distance_array);
	end_t = clock();

	fclose(list_raytrace);
	fclose(list_index);
	free(blocklist);
	free(raytrace);

	total_t = (double)(end_t-start_t)/CLOCKS_PER_SEC;
	// printf("Raytrace total time: %f\n", total_t);	
	return 0;
}

double raytracing(int i, int max, block* blocklist, vec source, int* psphere_counter, int* psphere_max_counter, 	block_distance *distance_sorter){
	/**/
	int j;
	long double raypath = 0, value;
	int obs_sign, block_j_sgn1, block_j_sgn2;
	vec normal;
	double distance, distance_max = 0, *pdistance_max = &distance_max;

	normal = normal_vec(source, blocklist[i].center);
	obs_sign = interacting_region(source, blocklist[i].center, normal);


	for(j=0; j < max; j++){

		if (i!=j){
			block_j_sgn1 = interacting_region(blocklist[j].center, blocklist[i].center, normal);	
			block_j_sgn2 = interacting_region(blocklist[j].center, source, normal);

			if ((obs_sign == block_j_sgn1) && (obs_sign == -block_j_sgn2)){
				value = distance_between_line_point(blocklist[i].center, source, blocklist[j].center);

				if (value == 0){
					raypath+=blocklist[j].r;
					distance = compute_distance(blocklist[j], blocklist[i]);
					preset_array(*psphere_counter, j , distance, distance_sorter, pdistance_max);
					(*psphere_counter)++;

				}
				else if ((value > 0) && (value < blocklist[j].r)){
					raypath += chord(blocklist[j].r, value);
					distance = compute_distance(blocklist[j], blocklist[i]);
					preset_array(*psphere_counter, j , distance, distance_sorter, pdistance_max);
					(*psphere_counter)++;

				}
			}					
		}
	}

	if ((*psphere_counter) > (*psphere_max_counter)) (*psphere_max_counter) = (*psphere_counter);		

	return raypath;
}  /* raytracing */

vec normal_vec(vec source, vec block_i){
	/* Computes the directional vector passing through two points. */
	vec normal;

	normal.x = block_i.x-source.x;
	normal.y = block_i.y-source.y;
	normal.z = block_i.z-source.z;
	
	return normal;
}	
	
double interacting_region(vec source, vec point_on_plane, vec normal){
	/*
	 * Computes the plane passing throught the center of the i-th blob and 
	 * orthogonal to the line passing through the observer and the i-th blob center */

	double distance;
	vec vec1;

	vec1.x = point_on_plane.x - source.x;
	vec1.y = point_on_plane.y - source.y;
	vec1.z = point_on_plane.z - source.z;
	/* Computing the dot product between vec1 and normal */
	distance = vec1.x*normal.x+vec1.y*normal.y+vec1.z*normal.z;
	distance = distance/fabsl(distance);
	return distance;
}   /* interacting_region */

double chord(double radius, double distance){
	/* Computes the chord of the blob determined by the intersection via the blob and the rayline*/
	return 2*sqrt(2*radius*distance-distance*distance);
}  /* chord */

double distance_between_line_point(vec p, vec q, vec r){
	/*
	 Returns the distance between a line, defined passing through two points, and a point.
        Computes the direction vector D of the line:
            D = p-q
        Computes the vector from point1 to point3.
            PR = r-p
        It then computes the projection of the vector on the direction vector:
            projection = PR x D/ D x D * D
        where 'x' stands for dot product
        Finally it normalize the distance of PR-projection
		*/
	double distance, dot1, dot2;
	vec D, PR, projection;

	D = normal_vec(p, q);

	PR.x = r.x-p.x;
	PR.y = r.y-p.y;
	PR.z = r.z-p.z;
	
	dot1 = PR.x*D.x+PR.y*D.y+PR.z*D.z;
	dot2 = D.x*D.x+D.y*D.y+D.z*D.z;
	
	projection.x = dot1/dot2*D.x;
	projection.y = dot1/dot2*D.y;
	projection.z = dot1/dot2*D.z;
	
	distance = sqrt(pow(PR.x-projection.x,2)+pow(PR.y-projection.y,2)+pow(PR.z-projection.z,2));
	return distance;
}


void preset_array(int sphere_counter, int sphere_index, double distance, block_distance *distance_array, double* pdistance_max){
	/*
	Updates the indices matrix, sorting them by distance from the blob i-th.
	*/

	if(sphere_counter < 2){
		if(sphere_counter == 0){
 		*pdistance_max = distance;
		// Initialize the sorting distance array:
		distance_array[sphere_counter].index = sphere_index;
		distance_array[sphere_counter].distance = distance;
		}
		else{
			if (distance > *pdistance_max){
				*pdistance_max = distance;
				// The sorting array is setted
				distance_array[sphere_counter].index = sphere_index;
				distance_array[sphere_counter].distance = distance;
			} 
			else{
				// Switching the two index position and setting the array:
				distance_array[sphere_counter].index = distance_array[0].index;
				distance_array[sphere_counter].distance = distance_array[0].distance;
				distance_array[0].index = sphere_index;
				distance_array[0].distance = distance;
			}
		}
	// The sorting array is now setted, with 0 having the lower distance from the blob i-th.
		return;
	}

	else{
		sorting_distance(distance_array, sphere_index, distance, sphere_counter);
		return;
	}

}


char* get_file_name(char* filepath){
	/*
	Detects the correct file name and define the new output name.
	*/
	char *filename;

    // Find the last occurrence of '/'
    filename = strrchr(filepath, '/');

    // If a '/' was found, increment the pointer to point to the filename
    if (filename != NULL) {
        filename++; // Move past the '/'
    } else {
        // If no '/' was found, the entire path is the filename
        filename = filepath;
    }
    // Print the filename
	filename = filename + 15; 
	return filename;
}

double compute_distance(block block_j, block block_i){
	/*
	Computes the distance between two blocks.
	*/
	double x,y,z, distance;
	x =  block_j.center.x - block_i.center.x;
	y =  block_j.center.y - block_i.center.y;
	z =  block_j.center.z - block_i.center.z;
	distance = sqrt(pow(x, 2)+pow(y,2)+pow(z,2));

	return distance;
}

void sorting_distance(block_distance* distance_array, int sphere_index, long double distance, int max_index){
	/* 
	Sorts the given array by distance and returns the place in
	*/
	int flag = 0, i=0, j;
		// Check the condition for the first element in the array:
		if (distance_array[0].distance > distance){
			for (j=max_index+1; j>0; j--){
				distance_array[j].distance = distance_array[j-1].distance;
				distance_array[j].index = distance_array[j-1].index;
			}
			distance_array[0].distance = distance;
			distance_array[0].index = sphere_index;
			flag ++;
		}
		// Check the condition for the last element in the array;
		if (distance_array[max_index-1].distance < distance){
			distance_array[max_index].distance = distance;
			distance_array[max_index].index = sphere_index;
			flag++;
		}
	while(flag == 0){
		if (distance_array[i].distance <= distance && distance <= distance_array[i+1].distance){
			for(j=max_index+1; j>i; j--){
				distance_array[j+1].distance = distance_array[j].distance;
				distance_array[j+1].index = distance_array[j].index;
			}
			distance_array[i+1].distance = distance;
			distance_array[i+1].index = sphere_index;
			flag++;
			return;
		}
		else{
			i++;
			if (i>max_index){
				printf("Something went wrong in sorting distance. I'm not supposed to be here...\n");
				flag++;
			}
		}
	}
	return;

}

void check_duplicates(block_distance* distance_sorted, int size){
	/* Only for debugging purposes*/
	int i, j;
	for (i = 0; i < size; i++) {
        for (j = i + 1; j < size; j++) {
            if (distance_sorted[i].index == distance_sorted[j].index) { // check if current element is equal to another
                break;
			}
		}
	}
}