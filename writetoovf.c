#ifndef WRITETOOVF_C
#define WRITETOOVF_C

#include <stdio.h>

int M,N,K;


void writetoovf(double ***m1,double ***m2,double ***m3,double dx,double dy,double dz,int iteration) {
    char filename[50]; // 假设文件名长度不超过49个字符
    sprintf(filename, "output_%d.ovf", iteration); // 生成带有迭代次数的文件名

    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        printf("Error opening file!\n");
        return;
    }

    fprintf(file, "# OOMMF: rectangular mesh v1.0\n");
    fprintf(file, "# Segment count: 1\n");
    fprintf(file, "# Begin: Segment\n");
    fprintf(file, "# Begin: Header\n");
    fprintf(file, "# Title: M-field\n");
    fprintf(file, "# meshtype: rectangular\n");
    fprintf(file, "# meshunit: m\n");
    fprintf(file, "# xstepsize: %e\n", dx);
    fprintf(file, "# ystepsize: %e\n", dy);
    fprintf(file, "# zstepsize: %e\n", dz);
    fprintf(file, "# xnodes: %d\n", M);
    fprintf(file, "# ynodes: %d\n", N);
    fprintf(file, "# znodes: %d\n", K);
    fprintf(file, "# End: Header\n");
    fprintf(file, "# Begin: Data Text\n");
    for (int k = 1; k <= K; k++) {
        for (int j = 1; j <= N; j++) {
            for (int i = 1; i <= M; i++) {
                fprintf(file, "%f %f %f\n", m1[i][j][k], m2[i][j][k], m3[i][j][k]);
            }
        }
    }
    fprintf(file, "# End: Data Text\n");
    fprintf(file, "# End: Segment\n");

    fclose(file);
}

#endif	