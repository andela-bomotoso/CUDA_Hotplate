/* A program to solve the hotplate problem  using GPU

  Author: Bukola Grace Omotoso
  MNumber: M01424979
  ID: bgo2e
  Last Modified: 11/27/2018
  
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#include<cuda_runtime.h>
#include<device_launch_parameters.h>
#include<cuda.h>



#define CHECK(call) \
{ \
 const cudaError_t error = call; \
 if (error != cudaSuccess) \
 { \
 printf("Error: %s:%d, ", __FILE__, __LINE__); \
 printf("code:%d, reason: %s\n", error, cudaGetErrorString(error)); \
 exit(1); \
 } \
}


float** buildHotplate(int rows, int columns) {
float** hotplate;
hotplate = (float**) malloc(rows*sizeof(float*));
for (int i = 0; i < rows; i++)
   hotplate[i] = (float*) malloc(columns*sizeof(float));
   return hotplate;
}

float* flattenArray(float** arrayToFlatten, int num_rows, int num_cols){
	float* flattenedArray = (float*) malloc(num_rows*num_cols*sizeof(float));
	int counter = 0;
	
	for (int row = 0; row < num_rows; row++){
		
	{
		for(int col = 0; col < num_cols; col++){
			flattenedArray[counter] = arrayToFlatten[row][col];
			counter++;
		}
	}
	
}
		return flattenedArray;
}

 void initializeHotPlate(int num_rows, int num_cols, float** hotplate, float** hotplateClone, int top_temp, int left_temp, int right_temp, int bottom_temp)	{
 	int num_outer_grid = (2 * num_rows) + (2 * (num_cols - 2));
	float outer_grid_sum = (top_temp * (num_cols - 2)) + (left_temp * (num_rows - 1)) + (bottom_temp * num_cols) + (right_temp * (num_rows - 1));
    float initial_inner_val = outer_grid_sum / num_outer_grid;

    for (int row = 0; row < num_rows; row++) {
            for (int column = 0; column < num_cols; column++) {

                //top values override the top row except the edges
                 if ((row == 0) & (column != 0 & column != num_cols - 1)) {
                	hotplate[row][column] = top_temp;
                    hotplateClone[row][column] = top_temp;
                } 
                else if (column == 0 && (row != (num_rows-1))) {
                    hotplate[row][column] = left_temp;
                    hotplateClone[row][column] = left_temp;
                }
                else if (column == (num_cols - 1) && (row != (num_rows-1))) {
                    hotplate[row][column] = right_temp;
                    hotplateClone[row][column] = right_temp;
                }
                else if(row == (num_rows -1 )){
                    hotplate[row][column] = bottom_temp;
                    hotplateClone[row][column] = bottom_temp;
                }
                if ((row != 0) && (row != num_rows - 1) && (column != 0) && (column != num_cols - 1))
                    hotplate[row][column] = initial_inner_val;
            }
        }

 }

 void swapHotplate(float *a, float *b) {

    float *tmp = a;
    a = b;
    b = tmp;
}


__global__ 
void generateHeat(int num_rows, int num_cols, float* hotplate, float* hotplateClone, float* d_maximums, float epsilon) {
        float max_difference = 0;
        float previous_val;
        float current_val;
        float diff;
        
         int row = blockIdx.x * blockDim.x + threadIdx.x;
         
            if (row > 0 && row < (num_rows-1)){
                for (int col = 1; col < (num_cols - 1); col++) {
					
				int idx = (row * num_cols) + col;
				float top = hotplate[idx - num_cols];				
				float bottom = hotplate[idx + num_cols];
				float left =  hotplate[idx - 1];
				float right =  hotplate[idx + 1];

                    previous_val = hotplate[idx];
                    current_val = (top + bottom + left + right) / 4.0;
                    diff = fabsf(previous_val - current_val);
                    if (diff > max_difference){
                        max_difference = diff;
                    }
                    hotplateClone[idx] = current_val;
                }
                
                d_maximums[row] = max_difference;
			}
}


/*Get the maximum values from all threads*/
float max_max_diff(float arr[], int n)
{
    int i;
        float max = arr[0];
    for (i = 1; i < n; i++)
        if (arr[i] > max)
            max = arr[i];
    return max;
}



int main(int argc, char const *argv[])
{		int num_rows = atoi(argv[1]);
		int num_cols = atoi(argv[2]);
		int top_temp = atoi(argv[3]);
		int left_temp = atoi(argv[4]);
		int right_temp = atoi(argv[5]);
		int bottom_temp = atoi(argv[6]);
		float epsilon = atof(argv[7]);
		
		float* flattenedhotplate;
		float* flattenedhotplateClone;
		float* maximums;

		
		int gridsize = 8;
		int block = 0;
		int block1 = (num_rows/gridsize);
		if (gridsize > num_rows)
			block = 1;
		else if ((block1 * gridsize) < num_rows){
			block = block1 + 1;
		}else{
			block = block1;
		}
		size_t nBytes = num_rows*num_cols * sizeof(float);
		double max_difference = epsilon + 1;
		int counter = 0;

        float** hotplate =  buildHotplate(num_rows, num_cols);
        float** hotplateClone = buildHotplate(num_rows, num_cols);
        
        
         
        initializeHotPlate(num_rows, num_cols, hotplate, hotplateClone, top_temp, left_temp, right_temp, bottom_temp);    
         
         flattenedhotplate = (float*) malloc(num_cols*num_rows*sizeof(float));
         flattenedhotplateClone = (float*) malloc(num_cols*num_rows*sizeof(float));
         maximums = (float*)malloc(num_rows*sizeof(float));
         
     
         flattenedhotplate = flattenArray(hotplate, num_rows, num_cols);
         flattenedhotplateClone = flattenArray(hotplateClone, num_rows, num_cols);

         
         float *d_hotplate; 
         float *d_hotplateClone;
         float *d_maximums;
	

		CHECK(cudaMalloc((float**)&d_hotplate, nBytes));
		CHECK(cudaMalloc((float**)&d_hotplateClone, nBytes));
		CHECK(cudaMalloc((float**)&d_maximums, num_rows*sizeof(float)));
		
		
		CHECK(cudaMemcpy( d_hotplate, flattenedhotplate, nBytes, cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_hotplateClone, flattenedhotplateClone, nBytes, cudaMemcpyHostToDevice));
		
		printf("%10s%10s\n", "Iteration", "Epsilon");

		 while(max_difference > epsilon){
        
		generateHeat<<<block,gridsize>>>(num_rows, num_cols, d_hotplate, d_hotplateClone, d_maximums, epsilon);
	
		cudaDeviceSynchronize();
		CHECK(cudaMemcpy(maximums, d_maximums, num_rows*sizeof(float), cudaMemcpyDeviceToHost));
		max_difference = max_max_diff(maximums, num_rows-1);
        
            float *T = d_hotplate;
            d_hotplate = d_hotplateClone;
            d_hotplateClone  =  T;
            
            if (counter > 0 && (counter & (counter - 1)) == 0)
                printf("%6d%15.6f\n", counter, max_difference);
            if (max_difference < epsilon) {
                printf("%6d%15.6f\n", counter, max_difference);
                break;
            }
            counter++;
	}
	
	cudaFree(d_hotplate);
	cudaFree(d_hotplateClone);
	cudaFree(d_maximums);
	
		return 0;

}

