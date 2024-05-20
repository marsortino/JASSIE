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

double distance_between_line_point(vec p, vec q, vec r);
vec normal_vec(vec source, vec block_i);
double interacting_region(vec source, vec point_on_plane, vec normal);
long double chord(double radius, double distance);
double raytracing(int block_num, int max, block* blocklist, vec source, int **indices, int* psphere_counter, int* psphere_counter_max);
void update_indices(int row, int col, int sphere_index, int **indices);
char* file_name(char* filepath);

int main(int argc, char *argv[]){
	
	if (argc != 5){
		fprintf(stderr, "Usage: %s <observer.x observer.y observer.z tmp_file_path>\n", argv[0]);
		return 1;
	}
	printf("Starting...\n");

	clock_t start_t, end_t;
	double total_t;
	start_t = clock();
	vec source;
	int i = 0, j, max, sphere_counter, sphere_counter_max = 0, *psphere_counter, *psphere_counter_max;
	FILE *list;
	block *blocklist;
	long double *raytrace;


	source.x = strtod(argv[1], NULL);
	source.y = strtod(argv[2], NULL);
	source.z = strtod(argv[3], NULL);

	char file_path[500];

	strcpy(file_path, argv[4]);


	puts(file_path);
	list = fopen(file_path, "r");

	//list = fopen("../../tmp_files/block_list_value.txt", "r");
	fscanf(list, "%i", &max);	

	blocklist = malloc(max*sizeof(block));

	while(!feof(list)){

		fscanf(list, "%le %le %le %le", &blocklist[i].center.x, &blocklist[i].center.y, &blocklist[i].center.z, &blocklist[i].r); 
		i++;
	}

	fclose(list);

	raytrace = malloc(max*sizeof(long double));
	*file_path = file_name(file_path);
	// We define a matrix n^2 made by -1. Each row represents a blob. The code computes which blob are in between the i-th blob and point P.
	// When this happens, the code updates the value in that row with the index of the blob that lies in between.

	int **indices = (int **)malloc(max * sizeof(int *));
	for(i = 0; i < max; i++){
		indices[i] = (int *)malloc(max * sizeof(int));
		memset(indices[i], -1, max * sizeof(int));
	}


	list = fopen("lib/tmp_files/tmp_raytrace_output.txt", "w");

	psphere_counter = &sphere_counter;
	psphere_counter_max = &sphere_counter_max;

	for(i=0; i<max; i++){
		sphere_counter = 0;
		raytrace[i] += raytracing(i, max, blocklist, source, indices, psphere_counter, psphere_counter_max);

		fprintf(list, "%0.15Le\n", raytrace[i]);

		//printf("block [%d] raytrace value: %Le\n", i, raytrace[i]);
	}
	end_t = clock();

	fclose(list);
	free(blocklist);
	free(raytrace);
	free(indices);
	total_t = (double)(end_t-start_t)/CLOCKS_PER_SEC;
	printf("Raytrace total time: %f\n", total_t);	
	return 0;
}

double raytracing(int i, int max, block* blocklist, vec source, int **indices, int* psphere_counter, int* psphere_max_counter){
	/**/
	int j;
	long double raypath = 0, value;
	int obs_sign, block_j_sgn1, block_j_sgn2;
	vec normal;


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
					update_indices(i,  *psphere_counter, j, indices);
					(*psphere_counter)++;
				}
				else if ((value > 0) && (value < blocklist[j].r)){
					raypath += chord(blocklist[j].r, value);
					update_indices(i,  *psphere_counter, j, indices);
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

long double chord(double radius, double distance){
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


void update_indices(int row, int col, int sphere_index, int **indices){
	/*
	Updates the indices matrix. 
	*/
	indices[row][col] = sphere_index;
}

char* file_name(char* filepath){
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
    printf("Filename: %s\n", filename);

    return filename;

}