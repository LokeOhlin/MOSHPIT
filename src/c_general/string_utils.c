#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cgeneral.h>
#include <ctype.h>
/*  
 *  Contains own defined methods for manipulating/using strings
 */


// Method for comparing strings of different lengths.
// Treats carriage returns as space such that
//      char string1[9]="string  ";
//      char string2[7]="string";
//      compStr(string1, string2) == 0;
// is true 
int compStr(char *s1, char*s2, size_t sz){
    while (sz != 0) {
        // At end of both strings, equal.
        if ((*s1 == '\0') && (*s2 == '\0')) break;

        // Treat spaces at end and end-string the same.
        if ((*s1 == '\0') && (*s2 == ' ')) { s2++; sz--; continue; }
        if ((*s1 == ' ') && (*s2 == '\0')) { s1++; sz--; continue; }

        // Detect difference otherwise.
        if (*s1 != *s2) return -1;

        s1++; s2++; sz--;
    }
    return 1;
}


void trim(char * s) {
    char * p = s;
    int l = strlen(p);

    while(isspace(p[l - 1])) p[--l] = 0;
    while(* p && isspace(* p)) ++p, --l;

    memmove(s, p, l + 1);
}   
