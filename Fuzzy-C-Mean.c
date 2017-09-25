#include <stdio.h>
#include <stdlib.h>

double patterns[150][5],original[150][5],data;
int c=3; //no. of clusters;
double m=9; //Fuzziness Index;
int n=0; //data points;
double ClusterCenter[3][4];

double EucleadianDistance(double x[][5],double cc[][4],int i,int j)
{
    int k=0;
    double distance = 0;
    for(k=0;k<4;k++)
    {
        distance+=pow((cc[j][k]-x[i][k]),2);

    }
    distance = sqrt(distance);
    return distance;
}

void displaycenters()
{

    int i;

    printf("Cluster 1 Centroid: ");
    for(i=0;i<4;i++)
    {
        printf("%lf ",ClusterCenter[0][i]);

    }
    printf("\nCluster 2 Centroid: ");
    for(i=0;i<4;i++)
    {
        printf("%lf ",ClusterCenter[1][i]);

    }
    printf("\nCluster 3 Centroid: ");
    for(i=0;i<4;i++)
    {
        printf("%lf ",ClusterCenter[2][i]);

    }
    printf("\n\n");
}

void retrieveData()
{
  FILE *fp;
  int i,j;


  fp = fopen("iris.txt", "r");
    fscanf (fp, "%lf", &data);

    j=0;
    i=0;
    n=0;

    while (!feof (fp))
    {

      if(j<4)
      {
          patterns[i][j]=data;
          j++;
      }
      else
      {
          j=0;
          i++;
          n++;
          patterns[i][j]=data;
          j++;
      }
      fscanf (fp, "%lf", &data);
    }

    n++;
    fclose (fp);


    fp = fopen("originaliris.txt", "r");
    fscanf (fp, "%lf", &data);

    j=0;
    i=0;

    while (!feof (fp))
    {

      if(j<5)
      {
          original[i][j]=data;
          j++;
      }
      else
      {
          j=0;
          i++;
          original[i][j]=data;
          j++;
      }
      fscanf (fp, "%lf", &data);
    }

      fclose (fp);
}

void initializeCentroids()
{
   double origin[1][4];
   int i,j,k;
   double d[150][2];
   double temp;
   int count[3];

   for(i=0;i<3;i++)
      count[i]=0;

   for(i=0;i<4;i++)
     origin[0][i]=0;

   for(i=0;i<n;i++)
   {
       d[i][0]= EucleadianDistance(patterns,origin,i,0);
       d[i][1]= i;
   }

   /*for(i=0;i<n;i++)
   {
       printf("distance %d from 0: %lf\n",i,d[i][0]);
   }*/

   for(i=0;i<n-1;i++)
   {
       for(j=i+1;j<n;j++)
       {
          if(d[i][0]>d[j][0])
          {
                for(k=0;k<2;k++){
                 temp = d[i][k];
                 d[i][k]=d[j][k];
                 d[j][k]=temp;
                }
          }
       }
   }

   /*for(i=0;i<n;i++)
   {
       printf("%d %lf\n",(int)d[i][1],d[i][0]);
   }*/

   /*for(i=0;i<n;i++)
   {
       printf("distance %lf from 0: %lf\n",d[i][1],d[i][0]);
   }*/

   for(i=0;i<c;i++)
   {
       for(j=0;j<4;j++)
         ClusterCenter[i][j]=0;
   }

   for(i=0;i<n;i++)
   {
       if(i<(int)n/3)
       {
           for(j=0;j<4;j++){
             ClusterCenter[0][j]+=patterns[(int)d[i][1]][j];
             //printf("%lf ",ClusterCenter[0][j]);
           }
           //printf("\n");
           count[0]++;


       }
       else if(i>=(int)n/3 && i<(int)(n*2)/3)
       {
           for(j=0;j<4;j++)
             ClusterCenter[1][j]+=patterns[(int)d[i][1]][j];
           count[1]++;
       }
       else if(i>=(int)(n*2)/3)
       {
           for(j=0;j<4;j++)
             ClusterCenter[2][j]+=patterns[(int)d[i][1]][j];
           count[2]++;
       }
   }

       for(j=0;j<c;j++)
         for(k=0;k<4;k++)
           ClusterCenter[j][k]=ClusterCenter[j][k]/count[j];

   /*for(i=0;i<n;i++)
   {
       if(i==24)
       {
          for(j=0;j<4;j++)
             ClusterCenter[0][j]=patterns[(int)d[i][1]][j];
       }
       else if(i==74)
       {
           for(j=0;j<4;j++)
             ClusterCenter[1][j]=patterns[(int)d[i][1]][j];
       }
       else if(i==124)
       {
           for(j=0;j<4;j++)
             ClusterCenter[2][j]+=patterns[(int)d[i][1]][j];
       }
   }*/





}

void FCM()
{
   int i,j,k,step=0;
   double vnum[3][4],vden[3][4];
   double oldJ=0,newJ=999999999;
   double u[150][3],oldu[150][3];
   double d[150][3],dk[150];
   double maxudiff = 0, maxu = 0;
   int confusionMatrix[3][3],max,clusterindex[3],match;
   double percentage = 0;
   double minJ,minU[150][3],minCC[3][4];
   double sumU[150];
   FILE *fp;

   retrieveData();
   initializeCentroids();
   printf("Initial Cluster Centers:\n\n");
   displaycenters();

   step=0;
   minJ=9999999999999999;
   while(step==0 || maxudiff>=0.00001)
   {
        step++;
        printf("\nSTEP: %d\n\n",step);

        //initialize centroid numerator and denominator
        for(i=0;i<c;i++)
        {
           for(j=0;j<4;j++)
           {
              vnum[i][j]=0; //centroid numerator
              vden[i][j]=0; //centroid denominator
           }
        }
        printf("Current Centroids:\n\n");
        displaycenters();


         for(i=0;i<n;i++)
         {
                for(j=0;j<c;j++)
                {
                    if(step==1)
                      oldu[i][j]=0;
                    else
                      oldu[i][j]=u[i][j];
                }
        }


        for(i=0;i<n;i++)
        {
            dk[i]=0;
            for(j=0;j<c;j++)
            {
                d[i][j] = EucleadianDistance(patterns,ClusterCenter,i,j);
                if(i==2)
                  printf("d %d %d = %lf\n",i,j,d[i][j]);
                dk[i]+=d[i][j];

            }
            if(i==2)
                  printf("dk %d = %lf\n",i,dk[i]);
            for(j=0;j<c;j++)
            {
                if(d[i][j]==0)
                    u[i][j]=999;
                else
                    u[i][j] = 1/(d[i][j]/dk[i]);
                u[i][j]=pow(u[i][j],(2/(m-1)));
                if(i==2)
                  printf("\nu %d %d = %lf\n",i,j,u[i][j]);
            }
        }

        //Normalize Us
        for(i=0;i<n;i++)
        {
            sumU[i]=0;
            for(j=0;j<c;j++)
            {
                sumU[i]+=u[i][j];
            }

        }
        for(i=0;i<n;i++)
        {
            for(j=0;j<c;j++)
            {
                u[i][j]=u[i][j]/sumU[i];
                if(i==2)
                  printf("\nnormalized u %d %d = %lf\n",i,j,u[i][j]);
            }
        }

        for(i=0;i<n;i++)
        {
            for(j=0;j<c;j++)
            {
                for(k=0;k<4;k++)
                {
                    vnum[j][k]+=pow(u[i][j],m)*patterns[i][k];
                    vden[j][k]+=pow(u[i][j],m);
                }
            }
        }
        for(i=0;i<c;i++)
        {
            for(j=0;j<4;j++)
            {
                ClusterCenter[i][j]=vnum[i][j]/vden[i][j];
            }
        }
        //displaycenters();

        oldJ=newJ;
        newJ=0;

        for(i=0;i<n;i++)
        {
            for(j=0;j<c;j++)
            {
                newJ+=pow(u[i][j],m)*pow(d[i][j],2);
            }
        }

        printf("\nOld Objective Function Value: %lf\n",oldJ);
        printf("New Objective Function Value: %lf\n",newJ);

        if(newJ<minJ)
        {
            minJ=newJ;
            for(i=0;i<n;i++)
            {
                for(j=0;j<c;j++)
                    minU[i][j]=u[i][j];
            }
            for(i=0;i<c;i++)
            {
                for(j=0;j<4;j++)
                    minCC[i][j]=ClusterCenter[i][j];
            }
        }
        maxudiff=0;
        for(i=0;i<n;i++)
        {
            for(j=0;j<c;j++)
            {
                if(i==2)
                   printf("%lf ", fabs(u[i][j]-oldu[i][j]));
                if(fabs(u[i][j]-oldu[i][j])>maxudiff)
                    maxudiff=fabs(u[i][j]-oldu[i][j]);

            }
            if(i==2)
             printf("\n");
        }

        printf("\nMaximum difference between previous and current membership degree: %lf\n",maxudiff);


        if(step==100)
            break;

   }

   for(i=0;i<n;i++)
   {
       for(j=0;j<c;j++)
         u[i][j]=minU[i][j];
   }
   for(i=0;i<c;i++)
   {
       for(j=0;j<4;j++)
          ClusterCenter[i][j]=minCC[i][j];
   }

   printf("\nU 1 0: %lf, U 1 1: %lf, U 1 2: %lf\n",u[1][0],u[1][1],u[1][2]);

   printf("\nFinal Cluster Centroids:\n\n");
   displaycenters();

    for(i=0;i<n;i++)
    {
        maxu=0;
        for(j=0;j<c;j++)
        {
           if(u[i][j]>maxu){
             maxu=u[i][j];
             patterns[i][4]=j;
           }
        }
    }

   for(i=0;i<c;i++)
   {
       for(j=0;j<c;j++)
        confusionMatrix[i][j]=0;
   }

   for(i=0;i<c;i++)
   {
       for(j=0;j<n;j++)
       {
           if(patterns[j][4]==i)
              confusionMatrix[i][(int)original[j][4]]+=1;
       }
   }

   printf("\nCONFUSION MATRIX:\n");
   for(i=0;i<c;i++)
   {
       for(j=0;j<c;j++)
        printf("%d ",confusionMatrix[i][j]);
       printf("\n");
   }

   for(i=0;i<c;i++)
   {
       max=0;
       for(j=0;j<c;j++)
       {
           if(confusionMatrix[i][j]>max)
           {
              max=confusionMatrix[i][j];
              clusterindex[i]=j;
           }
       }

   }

   for(i=0;i<n;i++)
   {
       patterns[i][4] = clusterindex[(int)patterns[i][4]];
   }


    fp = fopen ("irisoutput.txt", "w");

   //printf("\n\nMembers of Cluster 1 (I.setosa):\n");
   fputs("\n\nMembers of Cluster 1 (I.setosa):\n\n",fp);

   //printf("Press enter to continue...\n\n");
   //getch();
    k=0;
    for(i=0;i<n;i++)
    {
        if(patterns[i][4]==0){
        k++;
        fprintf(fp,"%d %s",k,") ");
        for(j=0;j<4;j++){
        //printf("%lf ",patterns[i][j]);
        fprintf(fp,"%lf %s",patterns[i][j]," ");
        }
        //printf("\n");
        fprintf(fp,"\n");

        }
    }

   //printf("\n\nMembers of Cluster 2 (I.versicolor):\n");
   fputs("\n\nMembers of Cluster 2 (I.versicolor):\n\n",fp);
   //printf("Press enter to continue...\n\n");
   //getch();
   k=0;
    for(i=0;i<n;i++)
    {
        if(patterns[i][4]==1){
        k++;
        fprintf(fp,"%d %s",k,") ");
        for(j=0;j<4;j++){
        //printf("%lf ",patterns[i][j]);
        fprintf(fp,"%lf %s",patterns[i][j]," ");
        }
        //printf("\n");
        fprintf(fp,"\n");

        }
    }

    //printf("\n\nMembers of Cluster 3 (I.verginica):\n");
    fputs("\n\nMembers of Cluster 3 (I.verginica):\n\n",fp);
    //printf("Press enter to continue...\n\n");
    //getch();
    k=0;
    for(i=0;i<n;i++)
    {
        if(patterns[i][4]==2){
        k++;
        fprintf(fp,"%d %s",k,") ");
        for(j=0;j<4;j++){
        //printf("%lf ",patterns[i][j]);
        fprintf(fp,"%lf %s",patterns[i][j]," ");
        }
        //printf("\n");
        fprintf(fp,"\n");
        }
    }

    match=0;
    for(i=0;i<n;i++)
    {
        if(patterns[i][4]==original[i][4])
            match+=1;
    }
    percentage =  (((double)match)/n)*100;

    printf("\n\nAccuracy of clustering: %lf Percent\n\n",percentage);
    fprintf(fp,"%s %lf %s","\n\nAccuracy of clustering: ",percentage," Percent");

    fclose (fp);


}

int main()
{
    FCM();
    return 0;
}
