#include <stdio.h>
#include <string.h>

int main(int argc, char const *argv[])
{
    FILE * source;
    FILE * destination;

    char str[] = "";

    int i,j,count,counter,length;

    int flag = 0;
    int recordflag = 0;

    int stepcount;

    stepcount = 0;
    count = 0;
    counter = -1;


    for(i=0;i<argc;i++)
    {
        if(strcmp(argv[i],"-s")==0)
        {    
            strcpy(str,argv[i+1]);
            source = fopen(str,"r");
        }

        if(strcmp(argv[i],"-d")==0)
        {    
            strcpy(str,argv[i+1]);
            source = fopen(str,"w");
        }

        if(strcmp(argv[i],"-c")==0)
        {
            strcpy(str,argv[i+1]);

            length = strlen(str);
            for(j=0;j<length;j++)
                count = count * 10 + ( str[i] - '0' );
        }
    }

    while(fscanf(source,"%s",str) != EOF){
        if(strcmp(str,"step")==0){
            flag = 0;
            recordflag = 0;
            counter++;

            fscanf(source,"%s",str);
            fscanf(source,"%s",str);

            length = strlen(str);
            for(j=0;j<length;j++)
                stepcount = stepcount * 10 + ( str[i] - '0' );


            fprintf(destination,"step = %d\n",stepcount);
        }

        if(counter == 0) flag = 1;
        if(counter >= count) counter = -1;

        if(strcmp(str,"Real")==0)
        {
            fscanf(source,"%s",str);
            recordflag = 1;
            fprintf(destination,"\nReal Part:\n");
        }

        if(strcmp(str,"Imaginary")==0)
        {
            fscanf(source,"%s",str);
            fprintf(destination,"\nImaginary Part:\n");
        }

        if(flag == 1 && recordflag == 1){
            fprintf(destination,"%s\n",str);
        }
    }
    fclose(destination);
    fclose(source);

}