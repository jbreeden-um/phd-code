/* Functions for Safety Constraints */
struct SafetyDataType SafetyData = {0};
/* Product of 1x3 vector with 3x4 matrix */
void VT3xM4(double A[3], double B[3][4], double C[4]){
   C[0] = A[0]*B[0][0] + A[1]*B[1][0] + A[2]*B[2][0];
   C[1] = A[0]*B[0][1] + A[1]*B[1][1] + A[2]*B[2][1];
   C[2] = A[0]*B[0][2] + A[1]*B[1][2] + A[2]*B[2][2];
   C[3] = A[0]*B[0][3] + A[1]*B[1][3] + A[2]*B[2][3];
}
/* Product of 1x3 vector with 3x3 matrix with 3x1 vector (returns a scalar) */
double VxMxV(double A[3], double M[3][3], double B[3]){
   return A[0]*M[0][0]*B[0] + A[1]*M[1][0]*B[0] + A[2]*M[2][0]*B[0]
        + A[0]*M[0][1]*B[1] + A[1]*M[1][1]*B[1] + A[2]*M[2][1]*B[1]
        + A[0]*M[0][2]*B[2] + A[1]*M[1][2]*B[2] + A[2]*M[2][2]*B[2];
}
/* Cross Product Matrix */
void Vx(double V[3], double M[3][3]){
   M[0][0] = 0.0;
   M[0][1] = -V[2];
   M[0][2] = V[1];
   M[1][0] = V[2];
   M[1][1] = 0.0;
   M[1][2] = -V[0];
   M[2][0] = -V[1];
   M[2][1] = V[0];
   M[2][2] = 0.0;
}
/* The following three functions describe generic pointing constraints */
double h_func(struct SCType *S, struct CmdVecType *CV, double theta){
   struct AcType *AC;
   double sB[3];
   AC = &S->AC;
   QxV(AC->qbn, CV->N, sB);
   return VoV(CV->R, sB) - cos(theta);
}
double hdot_func(struct SCType *S, struct CmdVecType *CV){
   struct AcType *AC;
   double sB[3], wx[3][3], sdotN[3], sdotB[3];
   AC = &S->AC;
   QxV(AC->qbn, CV->N, sB); /* Target vector in body frame */
   VxV(CV->wn, CV->N, sdotN); /* Change of target vector w.r.t N expressed in N */
   QxV(AC->qbn, sdotN, sdotB); /* Change of target vector w.r.t. N expressed in B */
   Vx(AC->wbn, wx);
   return VoV(sdotB, CV->R) + VxMxV(sB, wx, CV->R); /* Sum of variation of the target w.r.t. N and body w.r.t. N */
}
double phi_func(struct SCType *S, double u[4], struct CmdVecType *CV){
   struct AcType *AC;
   long i, j;
   double Z11[3][3] = {{5.982580878211763,   0.000000000000000,   0.000000000000000},
                       {-0.000000000000000,   7.940273264504397,   0.000000000000000},
                       {0.000000000000000,   0.000000000000000,  16.336748748482520}};
   double Z12[3][4] = {{0.000000000000000,   0.000000000000000,  -4.884756834288861,   4.884756834288861},
                       {0.000000000000000,   7.486161429636335,  -3.743080714816829,  -3.743080714816829},
                       {16.336748748482524,  -5.445582897368602,  -5.445582897384170,  -5.445582897384170}};
   double wxJw[3];
   double wdot[3] = {0.0, 0.0, 0.0};
   double wx[3][3], wx2[3][3], rx[3][3];
   double sB[3], sdotN[3], sdotB[3], sddotN[3], sddotB[3];

   AC = &S->AC;
   QxV(AC->qbn, CV->N, sB);
   VxV(CV->wn, CV->N, sdotN);
   QxV(AC->qbn, sdotN, sdotB);
   VxV(CV->wn, sdotN, sddotN);
   QxV(AC->qbn, sddotN, sddotB);

   VxV(AC->wbn, S->Hvb, wxJw);
   for (i=0; i<3; i++){
      for (j=0; j<3; j++){
         wdot[i] += Z11[i][j] * (-wxJw[j]);
      }
      for (j=0; j<4; j++){
         wdot[i] += Z12[i][j] * u[j];
      }
   }

   Vx(AC->wbn, wx);
   MxM(wx, wx, wx2);
   Vx(CV->R, rx);
   return VoV(sddotB, CV->R) + 2*VxMxV(sdotB, wx, CV->R) + VxMxV(sB, wx2, CV->R) - VxMxV(sB, rx, wdot);
}
/* Inverts the absSq function to determine required phi value */
double get_phiQ(struct SCType *S, struct CmdVecType *CV, double theta, double DT){
   const double delta2 = 1.103000000000000e-05;
   const double Delta2 = 1.103000000000000e-05;
   const double mu = 0.001670000000000;
   const double M2plus = 1.639470002367842e-04;
   const double M3plus = 0.006200000000000;
   double c_h, c_hdot, phi1, phi2, root1, root2;

   /* The following lines essentially invert an affine function */
   c_h = h_func(S, CV, theta) + hdot_func(S, CV)*DT + 1.0/2*M2plus*DT*DT + 1.0/2*M3plus*DT*DT*DT;
   phi1 = -2/(DT*DT)*(c_h + delta2);
   
   /* The follow lines essentially invert a monotone increasing function composed of two quadratic functions */
   c_hdot = hdot_func(S, CV) + M2plus*DT + 1.0/2*M3plus*DT*DT;
   root1 = pow(1.0/2*DT*DT + c_hdot*DT/mu, 2) - 4*( DT*DT/(2*mu))*(c_h + Delta2 + c_hdot*c_hdot/(2*mu));
   root2 = pow(1.0/2*DT*DT - c_hdot*DT/mu, 2) - 4*(-DT*DT/(2*mu))*(c_h + Delta2 - c_hdot*c_hdot/(2*mu));
   if (root1 > 0 && root2 > 0){
      // printf("Multiple solutions exist on line %d\n", __LINE__);
      double phi2_left, phi2_right, phi_crit;
      phi2_right = (-(1.0/2*DT*DT + c_hdot*DT/mu) + sqrt(root1)) / ( DT*DT/mu);
      phi2_left = (-(1.0/2*DT*DT - c_hdot*DT/mu) + sqrt(root2)) / (-DT*DT/mu);
      phi_crit = -c_hdot/DT; /* value for which the quadratic term is zero; note the slope is always DT*DT/2.0 at phi_crit */
      if (phi2_left >= phi_crit && phi2_right >= phi_crit){
         phi2 = phi2_right; /* The left side parabola is only valid when phi <= phi2_crit, so choose the right side */
      } else if (phi2_left <= phi_crit && phi2_right <= phi_crit){
         phi2 = phi2_left; /* The right side parabola is only valid when phi >= phi2_crit, so choose the left side */
      } else {
         printf("Ambiguous solution for phi2 on line %d\n", __LINE__);
      }
   } else if (root1 >= 0){
      phi2 = (-(1.0/2*DT*DT + c_hdot*DT/mu) + sqrt(root1)) / ( DT*DT/mu); /* Choose the larger root. Denominator is positive */
   } else if (root2 >= 0){
      phi2 = (-(1.0/2*DT*DT - c_hdot*DT/mu) + sqrt(root2)) / (-DT*DT/mu); /* Choose the smaller root. Denominator is negative */
   } else {
      printf("No solutuions exist on line %d\n", __LINE__);
      phi2 = 0;
   }
   return min(phi1, phi2);
}
/* Converts required phi value into A and b constraint matrix */
double get_ConQ(struct SCType *S, struct CmdVecType *CV, double theta, double A[4]){
   struct AcType *AC;
   double b;
   double sB[3], temp3[3];
   double Z12[3][4] = {{0.000000000000000,   0.000000000000000,  -4.884756834288861,   4.884756834288861},
                       {0.000000000000000,   7.486161429636335,  -3.743080714816829,  -3.743080714816829},
                       {16.336748748482524,  -5.445582897368602,  -5.445582897384170,  -5.445582897384170}};
   AC = &S->AC;
   b = get_phiQ(S, CV, theta, AC->DT) - phi_func(S, (double[4]){0.0, 0.0, 0.0, 0.0}, CV);
   QxV(AC->qbn, CV->N, sB);
   VxV(CV->R, sB, temp3); /* -(s^T)*(r^x) = (r x s)^T */
   VT3xM4(temp3, Z12, A);
   return b;
}
/* Returns all 4 constraints as a single matrix */
void GetConstraints(struct SCType *S, double A[4][4], double b[4]){
   struct AcType *AC;
   long i;

   // Ideally stuff like this would be read from a file, but for the moment I'm hardcoding it.
   double Je[3][3] = {{0.167151940000000,  -0.000000000000000,  -0.000000000000000},
                     {0.000000000000000,   0.125940250000000,  -0.000000000000000},
                     {-0.000000000000000,  -0.000000000000000,   0.061211690000000}};
   double Z12[3][4] = {{0.000000000000000,   0.000000000000000,  -4.884756834288861,   4.884756834288861},
                       {0.000000000000000,   7.486161429636335,  -3.743080714816829,  -3.743080714816829},
                       {16.336748748482524,  -5.445582897368602,  -5.445582897384170,  -5.445582897384170}};
   const double E_max = 5.091739267514272e-05;
   const double M2omega_lin = 8.300000000000000e-05;
   const double M1omega = 5.788735868341308e-07;
   const double theta_i = 25*D2R;
   const double theta_l = 30*D2R;
   const double theta_s = 45*D2R;
   
   double A_w[4], A_i[4], A_s[4], A_l[4];
   double b_w, b_i, b_s, b_l;
   double temp3[3];
   /* The pointing constraints are described as command vectors so we can keep track of the motion of the bright objects to be avoided */
   struct CmdVecType CV_s;
   struct CmdVecType CV_l;
   struct CmdVecType CV_i;

   /* Tracker to Sun Constraint */
   CV_s.Frame = FRAME_N;
   CV_s.Mode = CMD_TARGET;
   CV_s.TrgType = TARGET_WORLD;
   CV_s.TrgWorld = SOL;
   CV_s.R[0] = sin(10*D2R);
   CV_s.R[1] = -cos(10*D2R);
   CV_s.R[2] = 0.0;
   for(i=0;i<3;i++) CV_s.W[i] = 0.0;
   FindCmdVecN(S, &CV_s);
   SafetyData.h_s = h_func(S, &CV_s, theta_s);
   SafetyData.hdot_s = hdot_func(S, &CV_s);

   /* Tracker to Moon Constraint */
   CV_l.Frame = FRAME_N;
   CV_l.Mode = CMD_TARGET;
   CV_l.TrgType = TARGET_WORLD;
   CV_l.TrgWorld = LUNA;
   CV_l.R[0] = sin(10*D2R);
   CV_l.R[1] = -cos(10*D2R);
   CV_l.R[2] = 0.0;
   for(i=0;i<3;i++) CV_l.W[i] = 0.0;
   FindCmdVecN(S, &CV_l);
   SafetyData.h_l = h_func(S, &CV_l, theta_l);
   SafetyData.hdot_l = hdot_func(S, &CV_l);

   /* Instrument to Sun Constraint */
   CV_i.Frame = FRAME_N;
   CV_i.Mode = CMD_TARGET;
   CV_i.TrgType = TARGET_WORLD;
   CV_i.TrgWorld = SOL;
   CV_i.R[0] = 0.0;
   CV_i.R[1] = 0.0;
   CV_i.R[2] = 1.0;
   for(i=0;i<3;i++) CV_i.W[i] = 0.0;
   FindCmdVecN(S, &CV_i);
   SafetyData.h_i = h_func(S, &CV_i, theta_i);
   SafetyData.hdot_i = hdot_func(S, &CV_i);

   AC = &S->AC;

   /* Angular Velocity Constraint */
   MTxV(Je, AC->wbn, temp3);
   VT3xM4(temp3, Z12, A_w);
   for (i=0; i<4; i++) A_w[i] *= 2*AC->DT;
   b_w = E_max - VxMxV(AC->wbn, Je, AC->wbn) - M1omega*AC->DT - M2omega_lin*AC->DT*AC->DT/2.0;
   SafetyData.h_w = VxMxV(AC->wbn, Je, AC->wbn) - E_max;

   /* Pointing Constraints */
   b_i = get_ConQ(S, &CV_i, theta_i, A_i);
   b_s = get_ConQ(S, &CV_s, theta_s, A_s);
   b_l = get_ConQ(S, &CV_l, theta_l, A_l);

   memcpy(A[0], A_w, 4*sizeof(double));
   memcpy(A[1], A_i, 4*sizeof(double));
   memcpy(A[2], A_s, 4*sizeof(double));
   memcpy(A[3], A_l, 4*sizeof(double));
   b[0] = b_w;
   b[1] = b_i;
   b[2] = b_s;
   b[3] = b_l;
}

#include "osqp.h"
c_float QP_H_vals[4] = {1.0, 1.0, 1.0, 1.0};
c_int QP_H_rows[4] = {0, 1, 2, 3}; /* what row each H_val is in */
c_int QP_H_ptr[5] = {0, 1, 2, 3, 4}; /* the index of H_val at which each column starts, followed by the number of elements of Hvals */
c_float QP_F_vals[4];
c_float QP_A_vals[20]; /* values of A matrix one column at a time */
c_int QP_A_rows[20] = {0, 1, 2, 3, 4, 0, 1, 2, 3, 5, 0, 1, 2, 3, 6, 0, 1, 2, 3, 7};
c_int QP_A_ptr[5] = {0, 5, 10, 15, 20};
c_float QP_lower[8];
c_float QP_upper[8];
OSQPWorkspace *osqp_work = NULL;
OSQPSettings  *osqp_settings = NULL;
OSQPData      *osqp_data = NULL;

void call_QP(struct SCType *S, double A[4][4], double b[4], double u0[4], double u[4]) {
   long i, j;

   // Workspace structures
   if (osqp_work == NULL){
      QP_lower[0] = -1.0/0.0;
      QP_lower[1] = -1.0/0.0;
      QP_lower[2] = -1.0/0.0;
      QP_lower[3] = -1.0/0.0;
      QP_lower[4] = (c_float)(-S->Whl[0].Tmax);
      QP_lower[5] = (c_float)(-S->Whl[1].Tmax);
      QP_lower[6] = (c_float)(-S->Whl[2].Tmax);
      QP_lower[7] = (c_float)(-S->Whl[3].Tmax);
      QP_upper[4] = (c_float)(S->Whl[0].Tmax);
      QP_upper[5] = (c_float)(S->Whl[1].Tmax);
      QP_upper[6] = (c_float)(S->Whl[2].Tmax);
      QP_upper[7] = (c_float)(S->Whl[3].Tmax);
      QP_A_vals[4] = 1.0;
      QP_A_vals[9] = 1.0;
      QP_A_vals[14] = 1.0;
      QP_A_vals[19] = 1.0;

      osqp_settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));
      osqp_set_default_settings(osqp_settings);
      osqp_settings->alpha = 1.0; // Change alpha parameter
      osqp_settings->max_iter = 20; // deliberately set to a low number
      osqp_settings->verbose = 0;
      osqp_settings->warm_start = 1;
      osqp_settings->eps_abs = 1e-6;
      osqp_settings->eps_rel = 1e-6;
      
      osqp_data = (OSQPData *)c_malloc(sizeof(OSQPData));
      osqp_data->n = 4;
      osqp_data->m = 8;
      osqp_data->P = csc_matrix(4, 4, 4, QP_H_vals, QP_H_rows, QP_H_ptr); // quadratic term
      osqp_data->q = QP_F_vals; // linear term
      osqp_data->A = csc_matrix(8, 4, 20, QP_A_vals, QP_A_rows, QP_A_ptr); // constraint matrix
      osqp_data->l = QP_lower;
      osqp_data->u = QP_upper;

      // Setup workspace
      osqp_setup(&osqp_work, osqp_data, osqp_settings); 
   }
   for (i=0; i<4; i++){
      QP_upper[i] = (c_float)b[i];
      QP_F_vals[i] = -2*u0[i];
      for (j=0; j<4; j++){ // i is the row, j is the column
         QP_A_vals[i+5*j] = (c_float)A[i][j];
      }
   }

   /**
    * Passing OSQP_NULL tells OSQP to update the whole array.
    * We only need to update the first 16 elements (input constraints are constant),
    * but this is simpler than passing an array of 16 indices.
    * */
   osqp_update_A(osqp_work, QP_A_vals, OSQP_NULL, 20);
   osqp_update_lin_cost(osqp_work, QP_F_vals);
   osqp_update_bounds(osqp_work, QP_lower, QP_upper);

   // Solve Problem
   osqp_solve(osqp_work);
   // printf("OSQP Status: %d\n", work->info->status_val);
   // printf("OSQP Status: %s\n", work->info->status);

   for (i=0; i<4; i++) u[i] = (double)osqp_work->solution->x[i];
}

/**********************************************************************/
/*  This simple control law is suitable for rapid prototyping.        */
void PrototypeFSW(struct SCType *S)
{
      struct AcType *AC;
      struct AcPrototypeCtrlType *C;
      struct BodyType *B;
      struct CmdType *Cmd;
      // double alpha[3],Iapp[3];
      double Iapp[3];
      double Hvnb[3],Herr[3],werr[3];
      long Ig,i,j;

      AC = &S->AC;
      C = &AC->PrototypeCtrl;
      Cmd = &AC->Cmd;
            
      if (Cmd->Parm == PARM_AXIS_SPIN) {
         if (C->Init) {
            C->Init = 0;
            C->Kprec = 3.0E-2;
            C->Knute = 1.0;
         }
         
         SpinnerCommand(S);
         
         B = &S->B[0];
         
         MxV(B->CN,Cmd->Hvn,Hvnb);
         
         for(i=0;i<3;i++) {
            Herr[i] = S->Hvb[i] - Hvnb[i];
            werr[i] = AC->wbn[i] - Cmd->wrn[i];
            C->Tcmd[i] = -C->Knute*werr[i];
            if (MAGV(Herr) < 0.5*MAGV(Cmd->Hvn)) {
               C->Tcmd[i] -= C->Kprec*Herr[i];
            }
            AC->IdealTrq[i] = Limit(C->Tcmd[i],-0.1,0.1); 
         }
         
      }
      else {
         if (C->Init) {
            C->Init = 0;
            
            for(Ig=0;Ig<AC->Ng;Ig++) {
               FindAppendageInertia(Ig,S,Iapp);
               for(j=0;j<3;j++) {
                  FindPDGains(Iapp[j],0.05,1.0,
                     &AC->G[Ig].AngRateGain[j],
                     &AC->G[Ig].AngGain[j]);
                  AC->G[Ig].MaxAngRate[j] = 0.5*D2R;
                  AC->G[Ig].MaxTrq[j] = 0.1;
               }
            }

            for(i=0;i<3;i++) {
               FindPDGains(AC->MOI[i][i],0.25,0.7,&AC->ThreeAxisCtrl.Kr[i],&AC->ThreeAxisCtrl.Kp[i]);
            }
         }

         /* Find qrn, wrn and joint angle commands */
         int three_axis = 0;
         if (three_axis) ThreeAxisAttitudeCommand(S);

         double *sN;
         double *rB; // Instrument vector in the body frame is +Z axis
         double orth[3], sB[3];
         double angle;
         struct CmdVecType *PV;
         PV = &S->AC.Cmd.PriVec;
         FindCmdVecN(S, PV);
         sN = PV->N; // Target in the inertial frame
         rB = PV->R;
         QxV(AC->qbn,sN,sB); // Target in the body frame
         VxV(rB, sB, orth); // vector to rotate about in body frame
         UNITV(orth);
         angle = VoV(sB, rB);
         if (fabs(angle) > 1) angle = signum(angle);
         angle = min(acos(angle), 0.2); // rotation angle
         AC->qbr[3] = cos(angle/2); // this is the opposite of "dq" in order to work with 42 standard notations
         for (i=0; i<3; i++) AC->qbr[i] = sin(-angle/2)*orth[i];

         if (fabs(fmod(SimTime, 240.0) - 100.0) < 0.1) printf("New Vector: %lf %lf %lf\n", PV->N[0], PV->N[1], PV->N[2]);

         /* Form attitude error signals */
         if (three_axis) QxQT(AC->qbn,Cmd->qrn,AC->qbr);
         Q2AngleVec(AC->qbr,C->therr);
         for(i=0;i<3;i++) C->werr[i] = AC->wbn[i] - Cmd->wrn[i];

         /* Closed-loop attitude control */
         if (three_axis){
            //VectorRampCoastGlide(C->therr,C->werr,C->wc,C->amax,C->vmax,alpha);
            // for(i=0;i<3;i++) AC->IdealTrq[i] = AC->MOI[i][i]*alpha[i];
            VectorRampCoastGlide(C->therr,C->werr,C->wc,C->amax,C->vmax,C->Tcmd);
         }

         /**
          * For this spacecraft, we are in effect only doing 2 DOF attitude control.
          * The reason is that for this spacecraft, the solar arrays are assumed fixed (no rotation allowed),
          * and the solar arrays are opposite the instrument. Moreover, we assume pre-set data uplink windows.
          * As such, there is no reason for there to be 3 DOF attidue control. 
          * 3 DOF control can easily be achieved by switching three_axis to true,
          * but for this particular mission, it is assumed unnecessary. 
          * Note that switching to 3 DOF control could effect the presence/absence of local minima.
          * */

         double wheels_inv[3][4] = {{-0.000000000000000,   0.000000000000000,   0.612372435431858,  -0.612372435431858},
                                    {0.000000000000442,  -0.707106780881842,   0.353553390440574,   0.353553390440574},
                                    {-0.750000001293378,   0.249999999568814,   0.249999999568904,   0.249999999568904}};
         for (i=0; i<3; i++){
            if (three_axis)
               for(i=0;i<3;i++) C->Tcmd[i] *= AC->MOI[i][i];
            else
               C->Tcmd[i] = -AC->ThreeAxisCtrl.Kr[i]*C->werr[i] - AC->ThreeAxisCtrl.Kp[i]*C->therr[i];
         }
         for (i=0; i<4; i++){
            AC->Whl[i].Tcmd = 0;
            for (j=0; j<3; j++){
               AC->Whl[i].Tcmd += -wheels_inv[j][i]*C->Tcmd[j];
               // AC->IdealTrq[j] = C->Tcmd[j];
            }
         } 

         // printf("%lf %lf %lf %lf\n", POV.q[0], POV.q[1], POV.q[2], POV.q[3]);
         if (SimTime == 0.0){
            POV.q[0] = 0.932115;
            POV.q[1] = 0.025010;
            POV.q[2] = -0.359609;
            POV.q[3] = 0.034901;
         }

         double u0[4], A[4][4], b[4], u[4];
         for (i=0; i<4; i++) u0[i] = AC->Whl[i].Tcmd;
         GetConstraints(S, A, b);
         int use_qp = 0;
         if (use_qp){
            call_QP(S, A, b, u0, u);
            for (i=0; i<4; i++) AC->Whl[i].Tcmd = u[i];
            // printf("%lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf\n", u0[0], u[0], u0[1], u[1], u0[2], u[2], u0[3], u[3]);
         }
      }

}