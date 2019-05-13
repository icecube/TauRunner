#ifdef __cplusplus
// extern "C" {
#endif

void initJre();
void stopJre();
void setStderr(char*);
void initMMC(char*, int);
void deleteMMC();
void propagate(char*, double, double, double, double, double, double, double);
void createNext();
void endProp();
int getNum();
void startRead(int);
void endRead();
char *getName(int);
double getDouble(char*, int);
int getType();
void setStart(int);

#ifdef __cplusplus
// }
#endif
