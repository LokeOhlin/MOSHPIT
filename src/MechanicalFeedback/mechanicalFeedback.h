#include <cgeneral.h>
#ifndef FEEDBACK
#define FEEDBACK
#include <mechanicalFeedback.h>

int setFeedbackPars();
int checkFeedbackPars(char *name, char *value);
int initFeedback();
int doFeedbackStep(double dt, double *dt_feedback);
int Feedback_initIO();
int Feedback_output();
#endif
extern int nrealFBPars, nintFBPars;
extern int fb_useWind;
extern double fb_MdotWind, fb_vWind, fb_dEmaxWind;
extern real_list_t *FBDPars;
extern int_list_t  *FBIPars;
