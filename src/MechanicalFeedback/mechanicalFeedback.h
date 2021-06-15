#ifndef FEEDBACK
#define FEEDBACK
#include <mechanicalFeedback.h>

int setFeedbackPars();
int checkFeedbackPars(char *name, char *value);
int initFeedback();
int doFeedbackStep(double dt, double *dt_feedback);
#endif
extern int fb_useWind;
extern double fb_MdotWind, fb_vWind, fb_dEmaxWind;
