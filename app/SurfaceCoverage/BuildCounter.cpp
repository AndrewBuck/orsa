#include <iostream>
#include <stdio.h>

int main() {
  
  unsigned int counter = 0;
  
  FILE * fp = fopen("BuildCounter.dat","r");
  
  if (fp) {
    fscanf(fp,"%u",&counter);
    fclose(fp);
  }
  
  ++counter;
  
  fp = fopen("BuildCounter.dat","w");
  
  if (fp) {
    fprintf(fp,"%u\n",counter);
    fclose(fp);
  } else {
    printf("cannot increase build counter\n");
  }
  
  fp = fopen("BuildCounter.h","w");
  if (fp) {
    fprintf(fp,"#define __BUILD_COUNTER__ \"%u\"\n",counter);
    fclose(fp);
  } else {
    printf("cannot write header file\n");
  }
  
  printf("build %i\n",counter);
  
  return 0;
}
