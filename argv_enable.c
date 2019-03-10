#include <stdio.h>
#include <string.h>
#include <ctype.h>

#define MAX_ARGVINT 9999999

int mystr2int(char *strptr);
int mystr2int(char *strptr){
    int i,j;
    int buf = 0;
    int result = 0;
    int strlength;

    int anotherchance;

    if(isdigit(strptr[0])){
        strlength = strlen(strptr);
        for(i = 0; i < strlength; i++){
            buf = strptr[i] - '0';
            for(j = 1; j < strlength-i; j++){
                buf *= 10;
            }
            result += buf;
        }
        return result;
    }
    else{
        printf("Input Error. %s is not an POSITIVE INTEGER.\n Please Re-Enter.", strptr);
        exit(0);
    }
}
