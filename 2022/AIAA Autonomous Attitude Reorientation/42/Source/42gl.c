/**
  * In order to get the graphics shown in the animation, I needed to be able to display the sun and moon vectors.
  * Rather than making a new function to do this, I instead modified an existing function.
  * I simply commented out the B and H vectors in the TRUTH_VECTORS step of the graphics, and kept the sun vector.
  * I then added the moon vector to the same function.
  */

void DrawNearAuxObjects(void)
{
      for(Isc=0;Isc<Nsc;Isc++) {
         S = &SC[Isc];
         if (S->Exists) {
            glPushMatrix();

            if (ScIsVisible(POV.Host.RefOrb,Isc,PosR) ) {
               glTranslated(PosR[0],PosR[1],PosR[2]);
               if (Isc == POV.Host.SC) {
                  if (CamShow[TRUTH_VECTORS]) {
                     B = &S->B[0];
                     glPushMatrix();
                     glTranslated(B->pn[0],B->pn[1],B->pn[2]);
                     RotateL2R(World[Orb[POV.Host.RefOrb].World].CNH);
                     
                     for(i=0;i<3;i++) Vec[i] = World[SOL].PosH[i]-SC[POV.Host.SC].PosH[i];
                     UNITV(Vec);
                     DrawVector(Vec,"S"," ",SvbColor,AxisLength,
                        1.0,TRUE);
                        
                     for(i=0;i<3;i++) Vec[i] = World[LUNA].PosH[i]-SC[POV.Host.SC].PosH[i];
                     UNITV(Vec);
                     DrawVector(Vec,"L"," ",HvbColor,AxisLength,
                        1.0,TRUE);
                     glPopMatrix();
                  }
               }
            }
            glPopMatrix();
         }
      }

}