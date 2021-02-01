#ifndef CGENERAL
#define CGENERAL

#define MAX_STRING_LENGTH 80

#define MAX(a,b)((a>b)?a:b)
#define MIN(a,b)((a>b)?b:a)

// Lists of various types
typedef struct int_list_t{
  char name[MAX_STRING_LENGTH];
  int value;
} int_list_t;

typedef struct real_list_t{
  char name[MAX_STRING_LENGTH];
  double value;
} real_list_t;

typedef struct log_list_t{
  char name[MAX_STRING_LENGTH];
  int value; /*use int_list_t?*/
} log_list_t;

typedef struct str_list_t{
  char name[MAX_STRING_LENGTH];
  char value[MAX_STRING_LENGTH];
} str_list_t;

// Own defined string comparison. 
int compStr(char *s1, char*s2, size_t sz);
// Trim whitespaces and stuff
void trim(char *s);
// Progress bar
void DoProgress( int step, int total );
// Check boundry conditions. Imposes periodic
double checkBnds(double x, double xmin, double xmax);
// Create multiple dirs recursively
int recmkdir(const char *dir);
// Array search methods
int binarySearch(double val, double *array, int size);
int sequentialSearch(double val, double *array, int size);
#endif
