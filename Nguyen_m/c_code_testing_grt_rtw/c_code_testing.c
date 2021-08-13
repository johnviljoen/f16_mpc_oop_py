/*
 * c_code_testing.c
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "c_code_testing".
 *
 * Model version              : 1.2
 * Simulink Coder version : 9.5 (R2021a) 14-Nov-2020
 * C source code generated on : Tue Jun 29 20:40:44 2021
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "c_code_testing.h"
#include "c_code_testing_private.h"
#include "c_code_testing_dt.h"

/* Block signals (default storage) */
B_c_code_testing_T c_code_testing_B;

/* Continuous states */
X_c_code_testing_T c_code_testing_X;

/* External inputs (root inport signals with default storage) */
ExtU_c_code_testing_T c_code_testing_U;

/* External outputs (root outports fed by signals with default storage) */
ExtY_c_code_testing_T c_code_testing_Y;

/* Real-time model */
static RT_MODEL_c_code_testing_T c_code_testing_M_;
RT_MODEL_c_code_testing_T *const c_code_testing_M = &c_code_testing_M_;

/*
 * This function updates continuous states using the ODE3 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  /* Solver Matrices */
  static const real_T rt_ODE3_A[3] = {
    1.0/2.0, 3.0/4.0, 1.0
  };

  static const real_T rt_ODE3_B[3][3] = {
    { 1.0/2.0, 0.0, 0.0 },

    { 0.0, 3.0/4.0, 0.0 },

    { 2.0/9.0, 1.0/3.0, 4.0/9.0 }
  };

  time_T t = rtsiGetT(si);
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE3_IntgData *id = (ODE3_IntgData *)rtsiGetSolverData(si);
  real_T *y = id->y;
  real_T *f0 = id->f[0];
  real_T *f1 = id->f[1];
  real_T *f2 = id->f[2];
  real_T hB[3];
  int_T i;
  int_T nXc = 1;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);

  /* Save the state values at time t in y, we'll use x as ynew. */
  (void) memcpy(y, x,
                (uint_T)nXc*sizeof(real_T));

  /* Assumes that rtsiSetT and ModelOutputs are up-to-date */
  /* f0 = f(t,y) */
  rtsiSetdX(si, f0);
  c_code_testing_derivatives();

  /* f(:,2) = feval(odefile, t + hA(1), y + f*hB(:,1), args(:)(*)); */
  hB[0] = h * rt_ODE3_B[0][0];
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[0]);
  rtsiSetdX(si, f1);
  c_code_testing_step();
  c_code_testing_derivatives();

  /* f(:,3) = feval(odefile, t + hA(2), y + f*hB(:,2), args(:)(*)); */
  for (i = 0; i <= 1; i++) {
    hB[i] = h * rt_ODE3_B[1][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[1]);
  rtsiSetdX(si, f2);
  c_code_testing_step();
  c_code_testing_derivatives();

  /* tnew = t + hA(3);
     ynew = y + f*hB(:,3); */
  for (i = 0; i <= 2; i++) {
    hB[i] = h * rt_ODE3_B[2][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1] + f2[i]*hB[2]);
  }

  rtsiSetT(si, tnew);
  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model step function */
void c_code_testing_step(void)
{
  if (rtmIsMajorTimeStep(c_code_testing_M)) {
    /* set solver stop time */
    if (!(c_code_testing_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&c_code_testing_M->solverInfo,
                            ((c_code_testing_M->Timing.clockTickH0 + 1) *
        c_code_testing_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&c_code_testing_M->solverInfo,
                            ((c_code_testing_M->Timing.clockTick0 + 1) *
        c_code_testing_M->Timing.stepSize0 +
        c_code_testing_M->Timing.clockTickH0 *
        c_code_testing_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(c_code_testing_M)) {
    c_code_testing_M->Timing.t[0] = rtsiGetT(&c_code_testing_M->solverInfo);
  }

  /* Outport: '<Root>/Output' incorporates:
   *  Integrator: '<Root>/Integrator'
   */
  c_code_testing_Y.Output = c_code_testing_X.Integrator_CSTATE;

  /* Sum: '<Root>/Sum' incorporates:
   *  Inport: '<Root>/Input'
   *  Integrator: '<Root>/Integrator'
   */
  c_code_testing_B.Sum = c_code_testing_U.Input -
    c_code_testing_X.Integrator_CSTATE;
  if (rtmIsMajorTimeStep(c_code_testing_M)) {
    /* Matfile logging */
    rt_UpdateTXYLogVars(c_code_testing_M->rtwLogInfo,
                        (c_code_testing_M->Timing.t));
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep(c_code_testing_M)) {
    /* External mode */
    rtExtModeUploadCheckTrigger(2);

    {                                  /* Sample time: [0.0s, 0.0s] */
      rtExtModeUpload(0, (real_T)c_code_testing_M->Timing.t[0]);
    }

    if (rtmIsMajorTimeStep(c_code_testing_M)) {/* Sample time: [0.001s, 0.0s] */
      rtExtModeUpload(1, (real_T)(((c_code_testing_M->Timing.clockTick1+
        c_code_testing_M->Timing.clockTickH1* 4294967296.0)) * 0.001));
    }
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep(c_code_testing_M)) {
    /* signal main to stop simulation */
    {                                  /* Sample time: [0.0s, 0.0s] */
      if ((rtmGetTFinal(c_code_testing_M)!=-1) &&
          !((rtmGetTFinal(c_code_testing_M)-
             (((c_code_testing_M->Timing.clockTick1+
                c_code_testing_M->Timing.clockTickH1* 4294967296.0)) * 0.001)) >
            (((c_code_testing_M->Timing.clockTick1+
               c_code_testing_M->Timing.clockTickH1* 4294967296.0)) * 0.001) *
            (DBL_EPSILON))) {
        rtmSetErrorStatus(c_code_testing_M, "Simulation finished");
      }

      if (rtmGetStopRequested(c_code_testing_M)) {
        rtmSetErrorStatus(c_code_testing_M, "Simulation finished");
      }
    }

    rt_ertODEUpdateContinuousStates(&c_code_testing_M->solverInfo);

    /* Update absolute time for base rate */
    /* The "clockTick0" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick0"
     * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick0 and the high bits
     * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++c_code_testing_M->Timing.clockTick0)) {
      ++c_code_testing_M->Timing.clockTickH0;
    }

    c_code_testing_M->Timing.t[0] = rtsiGetSolverStopTime
      (&c_code_testing_M->solverInfo);

    {
      /* Update absolute timer for sample time: [0.001s, 0.0s] */
      /* The "clockTick1" counts the number of times the code of this task has
       * been executed. The resolution of this integer timer is 0.001, which is the step size
       * of the task. Size of "clockTick1" ensures timer will not overflow during the
       * application lifespan selected.
       * Timer of this task consists of two 32 bit unsigned integers.
       * The two integers represent the low bits Timing.clockTick1 and the high bits
       * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
       */
      c_code_testing_M->Timing.clockTick1++;
      if (!c_code_testing_M->Timing.clockTick1) {
        c_code_testing_M->Timing.clockTickH1++;
      }
    }
  }                                    /* end MajorTimeStep */
}

/* Derivatives for root system: '<Root>' */
void c_code_testing_derivatives(void)
{
  XDot_c_code_testing_T *_rtXdot;
  _rtXdot = ((XDot_c_code_testing_T *) c_code_testing_M->derivs);

  /* Derivatives for Integrator: '<Root>/Integrator' */
  _rtXdot->Integrator_CSTATE = c_code_testing_B.Sum;
}

/* Model initialize function */
void c_code_testing_initialize(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* initialize real-time model */
  (void) memset((void *)c_code_testing_M, 0,
                sizeof(RT_MODEL_c_code_testing_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&c_code_testing_M->solverInfo,
                          &c_code_testing_M->Timing.simTimeStep);
    rtsiSetTPtr(&c_code_testing_M->solverInfo, &rtmGetTPtr(c_code_testing_M));
    rtsiSetStepSizePtr(&c_code_testing_M->solverInfo,
                       &c_code_testing_M->Timing.stepSize0);
    rtsiSetdXPtr(&c_code_testing_M->solverInfo, &c_code_testing_M->derivs);
    rtsiSetContStatesPtr(&c_code_testing_M->solverInfo, (real_T **)
                         &c_code_testing_M->contStates);
    rtsiSetNumContStatesPtr(&c_code_testing_M->solverInfo,
      &c_code_testing_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&c_code_testing_M->solverInfo,
      &c_code_testing_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&c_code_testing_M->solverInfo,
      &c_code_testing_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&c_code_testing_M->solverInfo,
      &c_code_testing_M->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&c_code_testing_M->solverInfo, (&rtmGetErrorStatus
      (c_code_testing_M)));
    rtsiSetRTModelPtr(&c_code_testing_M->solverInfo, c_code_testing_M);
  }

  rtsiSetSimTimeStep(&c_code_testing_M->solverInfo, MAJOR_TIME_STEP);
  c_code_testing_M->intgData.y = c_code_testing_M->odeY;
  c_code_testing_M->intgData.f[0] = c_code_testing_M->odeF[0];
  c_code_testing_M->intgData.f[1] = c_code_testing_M->odeF[1];
  c_code_testing_M->intgData.f[2] = c_code_testing_M->odeF[2];
  c_code_testing_M->contStates = ((X_c_code_testing_T *) &c_code_testing_X);
  rtsiSetSolverData(&c_code_testing_M->solverInfo, (void *)
                    &c_code_testing_M->intgData);
  rtsiSetSolverName(&c_code_testing_M->solverInfo,"ode3");
  rtmSetTPtr(c_code_testing_M, &c_code_testing_M->Timing.tArray[0]);
  rtmSetTFinal(c_code_testing_M, 10.0);
  c_code_testing_M->Timing.stepSize0 = 0.001;

  /* Setup for data logging */
  {
    static RTWLogInfo rt_DataLoggingInfo;
    rt_DataLoggingInfo.loggingInterval = (NULL);
    c_code_testing_M->rtwLogInfo = &rt_DataLoggingInfo;
  }

  /* Setup for data logging */
  {
    rtliSetLogXSignalInfo(c_code_testing_M->rtwLogInfo, (NULL));
    rtliSetLogXSignalPtrs(c_code_testing_M->rtwLogInfo, (NULL));
    rtliSetLogT(c_code_testing_M->rtwLogInfo, "tout");
    rtliSetLogX(c_code_testing_M->rtwLogInfo, "");
    rtliSetLogXFinal(c_code_testing_M->rtwLogInfo, "");
    rtliSetLogVarNameModifier(c_code_testing_M->rtwLogInfo, "rt_");
    rtliSetLogFormat(c_code_testing_M->rtwLogInfo, 4);
    rtliSetLogMaxRows(c_code_testing_M->rtwLogInfo, 0);
    rtliSetLogDecimation(c_code_testing_M->rtwLogInfo, 1);
    rtliSetLogY(c_code_testing_M->rtwLogInfo, "");
    rtliSetLogYSignalInfo(c_code_testing_M->rtwLogInfo, (NULL));
    rtliSetLogYSignalPtrs(c_code_testing_M->rtwLogInfo, (NULL));
  }

  /* External mode info */
  c_code_testing_M->Sizes.checksums[0] = (1593011083U);
  c_code_testing_M->Sizes.checksums[1] = (2764943512U);
  c_code_testing_M->Sizes.checksums[2] = (3256060653U);
  c_code_testing_M->Sizes.checksums[3] = (248967043U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[1];
    c_code_testing_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(c_code_testing_M->extModeInfo,
      &c_code_testing_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(c_code_testing_M->extModeInfo,
                        c_code_testing_M->Sizes.checksums);
    rteiSetTPtr(c_code_testing_M->extModeInfo, rtmGetTPtr(c_code_testing_M));
  }

  /* block I/O */
  (void) memset(((void *) &c_code_testing_B), 0,
                sizeof(B_c_code_testing_T));

  /* states (continuous) */
  {
    (void) memset((void *)&c_code_testing_X, 0,
                  sizeof(X_c_code_testing_T));
  }

  /* external inputs */
  c_code_testing_U.Input = 0.0;

  /* external outputs */
  c_code_testing_Y.Output = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    c_code_testing_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 14;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.BTransTable = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.PTransTable = &rtPTransTable;
  }

  /* Matfile logging */
  rt_StartDataLoggingWithStartTime(c_code_testing_M->rtwLogInfo, 0.0,
    rtmGetTFinal(c_code_testing_M), c_code_testing_M->Timing.stepSize0,
    (&rtmGetErrorStatus(c_code_testing_M)));

  /* InitializeConditions for Integrator: '<Root>/Integrator' */
  c_code_testing_X.Integrator_CSTATE = c_code_testing_P.Integrator_IC;
}

/* Model terminate function */
void c_code_testing_terminate(void)
{
  /* (no terminate code required) */
}
