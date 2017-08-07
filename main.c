#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>

#define ROOT 0
#define MESSAGE_SIZE 50
#define END_OF_FILE "FIM"
#define SEED 666

#define ITERACOES 100
#define QTD_CENTROIDES 16

#define ARQUIVO "../coordinates.csv"
#define PTS_CENTROIDES "coordinates_centroides.csv"

#define QTD_PARAMETROS 6
#define LATITUDE 0
#define LONGITUDE 1
#define LATITUDE_ACUMULADA 2
#define LONGITUDE_ACUMULADA 3
#define QTD_PONTOS 4
#define SSE 5

double distancia(double latitude_centroide, double longitude_centroide, double latitude_ponto, double longitude_ponto) {
    double latitude = pow(latitude_centroide - latitude_ponto, 2);
    double longitude = pow(longitude_centroide - longitude_ponto, 2);
    return pow(latitude + longitude, 1.0 / 2.0);
}

void master(int myID, int qtdProcessos, MPI_Request *req, MPI_Status *status) {
    char mensagem[MESSAGE_SIZE];
    FILE *fp = fopen(ARQUIVO, "r");
    FILE *pt_centroides = fopen(PTS_CENTROIDES, "r");
    int slave;
    double centroides[QTD_CENTROIDES][QTD_PARAMETROS];
    double aux[QTD_CENTROIDES][QTD_PARAMETROS];
    //    srand(SEED);
    for (int i = 0; i < QTD_CENTROIDES; i++) {//melhorar random?
        fscanf(pt_centroides, " %lF %c %lF", &centroides[i][LATITUDE], mensagem, &centroides[i][LONGITUDE]);
        centroides[i][LATITUDE_ACUMULADA] = 0;
        centroides[i][LONGITUDE_ACUMULADA] = 0;
        centroides[i][QTD_PONTOS] = 0;
    }
    for (int i = 0; i < ITERACOES; i++) {
        MPI_Bcast(centroides, QTD_CENTROIDES*QTD_PARAMETROS, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
        fscanf(fp, " %s", mensagem); //descarte
        while (fscanf(fp, " %s", mensagem) != EOF) {//EOF = -1
            MPI_Recv(&slave, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, status);
            MPI_Isend(mensagem, MESSAGE_SIZE, MPI_CHAR, slave, 0, MPI_COMM_WORLD, req);
        }
        strcpy(mensagem, END_OF_FILE);
        for (int i = qtdProcessos - 1; i > 0; i--) {
            MPI_Recv(&slave, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, status);
            MPI_Isend(mensagem, MESSAGE_SIZE, MPI_CHAR, slave, 0, MPI_COMM_WORLD, req);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        for (int j = qtdProcessos - 1; j > 0; j--) {
            MPI_Recv(aux, QTD_CENTROIDES*QTD_PARAMETROS, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, status);
            for (int i = 0; i < QTD_CENTROIDES; i++) {
                centroides[i][SSE] += aux[i][SSE];
                centroides[i][LATITUDE_ACUMULADA] += aux[i][LATITUDE_ACUMULADA];
                centroides[i][LONGITUDE_ACUMULADA] += aux[i][LONGITUDE_ACUMULADA];
                centroides[i][QTD_PONTOS] += aux[i][QTD_PONTOS];
            }
        }
        for (int j = 0; j < QTD_CENTROIDES; j++) {
            if (centroides[j][QTD_PONTOS] > 0) {
                centroides[j][LATITUDE] = centroides[j][LATITUDE_ACUMULADA] / centroides[j][QTD_PONTOS];
                centroides[j][LONGITUDE] = centroides[j][LONGITUDE_ACUMULADA] / centroides[j][QTD_PONTOS];
                //                centroides[j][SSE] = centroides[j][SSE] / centroides[j][QTD_PONTOS];
            }
            if (i < ITERACOES - 1) {
                centroides[j][SSE] = 0;
                centroides[j][LATITUDE_ACUMULADA] = 0;
                centroides[j][LONGITUDE_ACUMULADA] = 0;
                centroides[j][QTD_PONTOS] = 0;
            }
        }
        rewind(fp);
    }
    double sseGlobal = .0;
    for (int i = 0; i < QTD_CENTROIDES; i++) {
        sseGlobal += centroides[i][SSE];
        printf("Latitude:%.2lF,Longitude:%.2lF\n", centroides[i][LATITUDE], centroides[i][LONGITUDE]);
    }
    printf("Erro:%lF\n\n", sseGlobal);

}

void slave(int myID, int qtdProcessos, MPI_Request *req, MPI_Status * status) {
    char mensagem[MESSAGE_SIZE];
    double latitude, longitude;
    double centroides[QTD_CENTROIDES][QTD_PARAMETROS];
    for (int i = 0; i < ITERACOES; i++) {
        MPI_Bcast(centroides, QTD_CENTROIDES*QTD_PARAMETROS, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
        while (strcmp(mensagem, END_OF_FILE)) {
            MPI_Send(&myID, 1, MPI_INT, ROOT, 0, MPI_COMM_WORLD);
            MPI_Recv(mensagem, MESSAGE_SIZE, MPI_CHAR, ROOT, 0, MPI_COMM_WORLD, status);
            if (strcmp(mensagem, END_OF_FILE)) {
                latitude = atof(strtok(mensagem, ","));
                longitude = atof(strtok(NULL, "\n"));
                int aux = 0;
                double menorDistancia = distancia(latitude, longitude, centroides[0][LATITUDE], centroides[0][LONGITUDE]);
                for (int i = 1; i < QTD_CENTROIDES; i++) {
                    if (menorDistancia > distancia(latitude, longitude, centroides[i][LATITUDE], centroides[i][LONGITUDE])) {
                        menorDistancia = distancia(latitude, longitude, centroides[i][LATITUDE], centroides[i][LONGITUDE]);
                        aux = i;
                    }
                }
                centroides[aux][SSE] += pow(menorDistancia, 2);
                centroides[aux][LATITUDE_ACUMULADA] += latitude;
                centroides[aux][LONGITUDE_ACUMULADA] += longitude;
                centroides[aux][QTD_PONTOS] += 1;
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Send(centroides, QTD_CENTROIDES*QTD_PARAMETROS, MPI_DOUBLE, ROOT, 0, MPI_COMM_WORLD);
        strcpy(mensagem, "");
    }
}

int main(int argc, char** argv) {

    int myID, qtdProcessos;

    MPI_Request req;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myID);
    MPI_Comm_size(MPI_COMM_WORLD, &qtdProcessos);


    if (myID == ROOT) {
        double time = MPI_Wtime();
        master(myID, qtdProcessos, &req, &status);
        printf("\nQuantidade de Centroides:%d", QTD_CENTROIDES);
        printf("\nQuantidade de iterações:%d", ITERACOES);
        printf("\nQuantidade de processos:%d", qtdProcessos);
        printf("\nTempo em segundos:%lF\n\n", MPI_Wtime() - time);
    } else {
        slave(myID, qtdProcessos, &req, &status);
    }

    MPI_Finalize();
    return (EXIT_SUCCESS);
}

