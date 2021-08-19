#include <stdio.h>
#include <stdlib.h>
int main() {
    int rc = 2;   // rc is row  number
    int n = 4;    // n is cloumn number
    int i = 0, j = 0;
    int **ptr2 = NULL;
    int cnt = 0;
    
    ptr2 = (int **)malloc(rc * sizeof(int));
    for (i = 0; i < rc; i++) ptr2[i] = (int *)malloc(n * sizeof(int));
    for (i = 0; i < rc; i++) {
      for(j = 0; j < n; j++) ptr2[i][j] = cnt++;
    }
    for (i = 0; i < rc; i++) {
      for(j = 0; j < n; j++) printf("ptr2[%d][%d] = %d",i,j,ptr2[i][j]);
      printf("\n");
    }
    for (i = 0; i < rc; i++) free(ptr2[i]);
    free(ptr2);
    return 0;

} 