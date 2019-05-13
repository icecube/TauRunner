#include <jni.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#define PATH_SEPARATOR ';'
#else /* UNIX */
#define PATH_SEPARATOR ':'
#endif

static JNIEnv *env;
static JavaVM *jvm;

static jclass Amanda, Particle;
static jobject amanda, particle, aobj, eobj;

static int flag;

void initJre(){

  jint res;
  JavaVMInitArgs vm_args;
  JavaVMOption options[4];

  options[0].optionString = "-Djava.compiler=JIT";          /* NONE disables JIT */
  options[1].optionString = "-Djava.class.path=mmc.jar";    /* user classes */
  options[2].optionString = "-Djava.library.path=.";        /* set native library path */
  options[3].optionString = "-verbose:none";                /* jni prints JNI-related messages */

  vm_args.version = JNI_VERSION_1_2;
  vm_args.options = options;
  vm_args.nOptions = 4;
  vm_args.ignoreUnrecognized = 1;

  res = JNI_CreateJavaVM(&jvm, (void **)&env, &vm_args);
  if (res < 0) {
    fprintf(stderr, "cannot create Java VM\n");
    exit(1);
  }

  fprintf(stderr, "jre started\n");

}


void stopJre(){

  (*jvm)->DestroyJavaVM(jvm);
  fprintf(stderr, "jre stopped\n");

}


void setStderr(char *filename){

  jclass Output;

  Output = (*env)->FindClass(env, "mmc/Output");
  if (Output == 0) {
    fprintf(stderr, "cannot find the mmc Output class\n");
    exit(1);
  }

  jmethodID mid = (*env)->GetStaticMethodID(env, Output, "setStderr", "(Ljava/lang/String;)Z");
  if (mid == 0) {
    fprintf(stderr, "cannot find the setStderr method\n");
    exit(1);
  }

  jstring jstr = (*env)->NewStringUTF(env, filename);
  if (jstr == 0) {
    fprintf(stderr, "Out of memory\n");
    exit(1);
  }

  (*env)->CallStaticBooleanMethod(env, Output, mid, jstr);
  (*env)->DeleteLocalRef(env, jstr);
  fprintf(stderr, "Stderr redirected to the file %s\n", filename);

}


void initMMC(char *options, int iflag){

  jmethodID mid;
  jstring jstr;

  flag=iflag;
  Amanda = (*env)->FindClass(env, flag==1?"tfa/Amanda":"gen/AtmFlux");
  if (Amanda == 0) {
    fprintf(stderr, "cannot find the main class\n");
    exit(1);
  }

  mid = (*env)->GetMethodID(env, Amanda, "<init>", "()V");
  if (mid == 0) {
    fprintf(stderr, "cannot find the constructor\n");
    exit(1);
  }

  jobject amanda_l = (*env)->NewObject(env, Amanda, mid);
  if (amanda_l == 0) {
    fprintf(stderr, "cannot start the constructor\n");
    exit(1);
  }

  amanda = (*env)->NewGlobalRef(env, amanda_l);
  if (amanda == 0) {
    fprintf(stderr, "cannot create global reference to amanda\n");
    exit(1);
  }
  (*env)->DeleteLocalRef(env, amanda_l);

  mid = (*env)->GetMethodID(env, Amanda, "setup", "(Ljava/lang/String;)V");
  if (mid == 0) {
    fprintf(stderr, "cannot find the setup method\n");
    exit(1);
  }

  jstr = (*env)->NewStringUTF(env, options);
  if (jstr == 0) {
    fprintf(stderr, "Out of memory\n");
    exit(1);
  }

  Particle = (*env)->FindClass(env, "mmc/Particle");
  if (Particle == 0) {
    fprintf(stderr, "cannot find the Particle class\n");
    exit(1);
  }

  (*env)->CallVoidMethod(env, amanda, mid, jstr);
  (*env)->DeleteLocalRef(env, jstr);
  fprintf(stderr, "MMC initialized\n");

}


void deleteMMC(){
  (*env)->DeleteGlobalRef(env, amanda);
  fprintf(stderr, "MMC deleted\n");
}


void propagate(char *name, double x, double y, double z, double theta, double phi, double e, double t){

  jmethodID pmid, amid;
  jstring jstr;

  // create new particle

  pmid = (*env)->GetMethodID(env, Particle, "<init>", "(Ljava/lang/String;DDDDDDD)V");
  if (pmid == 0) {
    fprintf(stderr, "cannot find the Particle constructor\n");
    exit(1);
  }

  jstr = (*env)->NewStringUTF(env, name);
  if (jstr == 0) {
    fprintf(stderr, "Out of memory\n");
    exit(1);
  }

  jobject particle_l = (*env)->NewObject(env, Particle, pmid, jstr, x, y, z, theta, phi, e, t);
  if (particle_l == 0) {
    fprintf(stderr, "cannot start the Particle constructor\n");
    exit(1);
  }
  (*env)->DeleteLocalRef(env, jstr);

  particle = (*env)->NewGlobalRef(env, particle_l);
  if (particle == 0) {
    fprintf(stderr, "cannot create global reference to the initial particle\n");
    exit(1);
  }
  (*env)->DeleteLocalRef(env, particle_l);

  // propagate the particle

  amid = (*env)->GetMethodID(env, Amanda, "propagate", "(Lmmc/Particle;)[Lmmc/Particle;");
  if (amid == 0) {
    fprintf(stderr, "cannot find the propagate method\n");
    exit(1);
  }

  jobject aobj_l=(*env)->CallObjectMethod(env, amanda, amid, particle);
  if (aobj_l == 0) {
    fprintf(stderr, "cannot run the propagate method\n");
    exit(1);
  }

  aobj = (*env)->NewGlobalRef(env, aobj_l);
  if (aobj == 0) {
    fprintf(stderr, "cannot create global reference to the particle array\n");
    exit(1);
  }
  (*env)->DeleteLocalRef(env, aobj_l);
}


void createNext(){
  jmethodID amid;
  amid = (*env)->GetMethodID(env, Amanda, "createNext", "()[Lmmc/Particle;");
  if (amid == 0) {
    fprintf(stderr, "cannot find the createNext method\n");
    exit(1);
  }

  jobject aobj_l=(*env)->CallObjectMethod(env, amanda, amid);
  if (aobj_l == 0) {
    fprintf(stderr, "cannot run the createNext method\n");
    exit(1);
  }

  aobj = (*env)->NewGlobalRef(env, aobj_l);
  if (aobj == 0) {
    fprintf(stderr, "cannot create global reference to the particle array\n");
    exit(1);
  }
  (*env)->DeleteLocalRef(env, aobj_l);
}


void endProp(){
  if(flag!=2) (*env)->DeleteGlobalRef(env, particle);
  (*env)->DeleteGlobalRef(env, aobj);
}


int getNum(){
  return (*env)->GetArrayLength(env, (jobjectArray) aobj);
}

void startRead(int i){
  jobject eobj_l = (*env)->GetObjectArrayElement(env, (jobjectArray) aobj, i);
  if (eobj_l == 0) {
    fprintf(stderr, "cannot initialize the particle array element\n");
    exit(1);
  }

  eobj = (*env)->NewGlobalRef(env, eobj_l);
  if (eobj == 0) {
    fprintf(stderr, "cannot create global reference to the particle array element\n");
    exit(1);
  }
  (*env)->DeleteLocalRef(env, eobj_l);
}

void endRead(){
  (*env)->DeleteGlobalRef(env, eobj);
}

const char *getName(int i){
  static char name[256];
  jobject sobj=(*env)->GetObjectField(env, i==0?particle:eobj, (*env)->GetFieldID(env, Particle, "name", "Ljava/lang/String;"));
  const char *aname=(*env)->GetStringUTFChars(env, (jstring) sobj, NULL);
  strncpy(name, aname, 255);
  (*env)->ReleaseStringUTFChars(env, (jstring) sobj, aname);
  (*env)->DeleteLocalRef(env, sobj);
  return name;
}

double getDouble(char *s, int i){
  switch(i){
  case 0: return (*env)->GetDoubleField(env, eobj, (*env)->GetFieldID(env, Particle, s, "D"));
  case 1: return (*env)->GetDoubleField(env, particle, (*env)->GetFieldID(env, Particle, s, "D"));
  case 2: return (*env)->GetDoubleField(env, amanda, (*env)->GetFieldID(env, Amanda, s, "D"));
  default: return 0;
  }
}

int getType(){
  return (*env)->GetIntField(env, eobj, (*env)->GetFieldID(env, Particle, "type", "I"));
}

void setStart(int f){
  if(getDouble("rw", 2)!=0) (*env)->SetBooleanField(env, amanda, (*env)->GetFieldID(env, Amanda, "dw", "Z"), f);
}
