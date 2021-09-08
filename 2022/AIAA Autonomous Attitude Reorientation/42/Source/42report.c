/**
  * We added three new files to 42's outputs.
  * Each file needs to be represented in the below function in three places, copied below.
  * */

void Report(void)
{
      static FILE *Safetyfile;
      static FILE *AttErrfile;
      static FILE *Whlfile;
	  
      if (First) {
         Safetyfile = FileOpen(InOutPath,"Safety.42","w");
         AttErrfile = FileOpen(InOutPath,"AttErr.42","w");
         Whlfile = FileOpen(InOutPath,"Whl.42","w");
      }

      if (OutFlag) {
         if (SC[0].Exists) {
            fprintf(Safetyfile,"%lf %lf %lf %lf %lf %lf %lf\n",SafetyData.h_i,SafetyData.hdot_i,
               SafetyData.h_s,SafetyData.hdot_s,SafetyData.h_l,SafetyData.hdot_l,SafetyData.h_w);
            fprintf(AttErrfile,"%lf %lf %lf\n",SC[0].AC.PrototypeCtrl.therr[0],
               SC[0].AC.PrototypeCtrl.therr[1],SC[0].AC.PrototypeCtrl.therr[2]);
            for (i=0; i<SC[0].Nw; i++){
               fprintf(Whlfile,"%lf %lf %lf",SC[0].Whl[i].w,SC[0].Whl[i].Trq,SC[0].AC.Whl[i].Tcmd);
               if (i<SC[0].Nw-1) fprintf(Whlfile," ");
               else fprintf(Whlfile,"\n");
            }
         }
      }
}